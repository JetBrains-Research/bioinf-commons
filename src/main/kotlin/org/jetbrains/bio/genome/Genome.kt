package org.jetbrains.bio.genome

import com.google.common.annotations.VisibleForTesting
import com.google.common.collect.ComparisonChain
import com.google.common.collect.Maps
import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.genome.sequence.TwoBitReader
import org.jetbrains.bio.genome.sequence.TwoBitSequence
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.checkOrRecalculate
import org.jetbrains.bio.util.div
import java.io.IOException
import java.io.UncheckedIOException
import java.lang.ref.WeakReference
import java.nio.file.Path

/**
 * The genome.
 *
 * Currently supported builds are:
 *
 *   mm9, mm10
 *   hg18, hg19, hg38
 *
 * A [Genome] is completely defined by the UCSC build string, however,
 * please use [Genome] instead of the build string for all public methods.
 *
 * If your API is capable of processing a subset of chromosomes, consider
 * using [GenomeQuery] to indicate that.
 *
 * @author Sergei Lebedev
 */

class Genome private constructor(
        /** Build in UCSC nomenclature, e.g. `"mm9"`. */
        val build: String,
        /** Absolute path to the genome chromosome sizes path,
         * if **not exists** will be downloaded automatically. */
        val chromSizesPath: Path = Configuration.genomesPath / build / "$build.chrom.sizes"
) : Comparable<Genome> {

    /** Absolute path to the genome data folder. */
    val dataPath = Configuration.genomesPath / build

    /** Species token, e.g. `"mm"`. */
    val species: String get() = build.takeWhile { !it.isDigit() }

    /**
     * Content of chrom.sizes file. Preserves original order during iteration.
     */
    internal val chromSizesMap: LinkedHashMap<String, Int> by lazy {
        // Fetch chrom sizes path if not exists
        chromSizesPath.checkOrRecalculate { (p) ->
            AnnotationsConfig[build].chromsizesUrl.downloadTo(p)
        }
        LOG.debug("Loading chrom.sizes $chromSizesPath")
        val map = LinkedHashMap<String, Int>()
        CSVFormat.TDF.parse(chromSizesPath.bufferedReader()).use { `in` ->
            `in`.records.forEach {
                try {
                    map[it[0]] = it[1].toInt()
                } catch (t: Throwable) {
                    LOG.error("Failed to parse chrom.sizes file: $chromSizesPath, line: ${it.joinToString("\t")}", t)
                }
            }
            LOG.debug("DONE Loading chrom.sizes $chromSizesPath")
            return@lazy map
        }
    }

    val chromosomes: List<Chromosome> by lazy {
        chromSizesMap.map { (name, length) ->
            Chromosome.getOrCreateChromosome(this, name, length)
        }
    }

    /**
     * A resolver for chromosome names.
     *
     * Currently support the following names:
     *
     *     chr21
     *     CHRX
     *     21 or X
     *
     * @author Sergei Lebedev
     */
    val chromosomeNamesMap: Map<String, Chromosome> by lazy {
        val map = hashMapOf<String, Chromosome>()
        for (chr in chromosomes) {
            val name = chr.name
            map[name.substringAfter("chr")] = chr
            map[name] = chr
            map[name.toLowerCase()] = chr
        }
        map
    }

    val twoBitPath: Path
        get() {
            val path = dataPath / "$build.2bit"
            path.checkOrRecalculate { output ->
                AnnotationsConfig[build].sequenceUrl.downloadTo(output.path)
            }
            return path
        }

    val transcripts: Collection<Transcript> by lazy { Transcripts.all(this).values() }

    val genes: List<Gene> by lazy { groupTranscripts(transcripts) }

    val gtfAnnotationsPath: Path
        get() {
            val annotationsPath = dataPath / "$build.annotations.gtf"
            annotationsPath.checkOrRecalculate("GTF annotations") { (path) ->
                Ensembl.getGTF(this).bufferedReader().use {
                    Ensembl.convertGTF(this, it, path)
                }
            }
            return annotationsPath
        }

    override fun compareTo(other: Genome) =
            ComparisonChain.start().compare(build, other.build).compare(chromSizesPath, other.chromSizesPath).result()

    fun presentableName(): String {
        val config = AnnotationsConfig[build]
        return "${config.species}: $build${when {
            config.alias != null -> " (${config.alias})"
            else -> ""
        }}"
    }

    companion object {
        const val TEST_ORGANISM_BUILD = "to1"

        @VisibleForTesting
        internal val LOG = Logger.getLogger(Genome::class.java)

        /** Use cache to avoid extra chrom.sizes loading. */
        private val CACHE = Maps.newConcurrentMap<Pair<String, Path?>, Genome>()

        operator fun get(build: String, chromSizesPath: Path? = null) =
                CACHE.computeIfAbsent(build to chromSizesPath) {
                    if (chromSizesPath != null) Genome(build, chromSizesPath) else Genome(build)
                }!!
    }
}

/**
 * The Chromosome.
 *
 * @author Sergei Lebedev
 */
@Suppress("DataClassPrivateConstructor") // Constructor is used in [invoke] call
data class Chromosome private constructor(
        /** Reference to the owner genome. */
        val genome: Genome,
        /** Unique chromosome name usually prefixed by `"chr"`, e.g. `"chr19"`. */
        val name: String,
        /** Length defined in chrom.sizes file. */
        val length: Int) : Comparable<Chromosome> {

    /**
     * Weak reference for sequence caching
     */
    @Transient
    private var sequenceRef: WeakReference<TwoBitSequence> = WeakReference<TwoBitSequence>(null)

    val sequence: TwoBitSequence
        @Synchronized get() {
            var s = sequenceRef.get()
            if (s == null) {
                s = try {
                    val twoBitLength = TwoBitReader.length(genome.twoBitPath, name)
                    check(length == twoBitLength) {
                        "Chromosome $name length differs in chrom.sizes($length) and 2bit file($twoBitLength)"
                    }
                    TwoBitReader.read(genome.twoBitPath, name)
                } catch (e: IOException) {
                    throw UncheckedIOException(
                            "Error loading $name from ${genome.twoBitPath}", e)
                }

                sequenceRef = WeakReference(s)
            }

            return s
        }

    val range: Range get() = Range(0, length)

    val centromere: Range?
        get() {
            val centromeres = gaps.asSequence()
                    .filter(Gap::isCentromere).map(Gap::location)
                    .toList()
            return when (centromeres.size) {
                0 -> null  // "no centromeres found for this chr/species, e.g. chrM or D. melanogaster
                1 -> centromeres.first().toRange()
                else -> {
                    // Standards are for breaking! hg38 uses a completely
                    // different approach to centromere annotations. Each
                    // centromere is split into multiple consequent chunks.
                    // Thus the "whole" centromere can be obtained as a
                    // bounding range.
                    val startOffset = centromeres.map { it.startOffset }.min()!!
                    val endOffset = centromeres.map { it.endOffset }.max()!!
                    Range(startOffset, endOffset)
                }
            }
        }

    val transcripts: List<Transcript> get() = Transcripts.all(genome)[this]

    val genes: List<Gene> by lazy { groupTranscripts(transcripts) }

    val repeats: List<Repeat> get() = Repeats.all(genome)[this]

    val gaps: List<Gap> get() = Gaps.all(genome)[this]

    val cytoBands: List<CytoBand> get() = CytoBands.all(genome)[this]

    val cpgIslands: List<CpGIsland> get() = CpGIslands.all(genome)[this]

    override fun compareTo(other: Chromosome) = ComparisonChain.start()
            .compare(genome, other.genome)
            .compare(name, other.name)
            .result()

    override fun toString() = "${genome.build}:$name"

    companion object {
        /** Use cache to avoid extra sequence loading. */
        private val CACHE = Maps.newConcurrentMap<Pair<Genome, String>, Chromosome>()

        /** Constructs a chromosome from a human-readable name, e.g. "chrM". */
        operator fun invoke(genome: Genome, name: String): Chromosome {
            check(name in genome.chromosomeNamesMap) {
                "unknown chromosome $name for genome ${genome.presentableName()}"
            }
            val canonicalName = genome.chromosomeNamesMap[name]!!.name
            val length = genome.chromSizesMap[canonicalName]!!
            return getOrCreateChromosome(genome, canonicalName, length)
        }

        internal fun getOrCreateChromosome(genome: Genome, canonicalName: String, length: Int): Chromosome {
            return CACHE.computeIfAbsent(genome to canonicalName) { Chromosome(genome, canonicalName, length) }
        }

        internal val ADAPTER = object : TypeAdapter<Chromosome>() {
            override fun read(`in`: JsonReader) = with(`in`) {
                val token = nextString()
                val build = token.substringBefore(":")
                val name = token.substringAfter(":")
                invoke(Genome[build], name)
            }

            override fun write(out: JsonWriter, chromosome: Chromosome) {
                out.value("${chromosome.genome.build}:${chromosome.name}")
            }
        }.nullSafe()
    }
}


/**
 * DNA strand.
 *
 * @author Roman Chernyatchik
 */
enum class Strand(val char: Char) {
    PLUS('+'), MINUS('-');

    fun asFilter() = when (this) {
        PLUS -> StrandFilter.PLUS
        MINUS -> StrandFilter.MINUS
    }

    fun isPlus() = this === PLUS

    fun isMinus() = this === MINUS

    fun opposite() = choose(MINUS, PLUS)

    @Suppress("NOTHING_TO_INLINE")  // Otherwise the arguments are boxed.
    inline fun <T> choose(ifPlus: T, ifMinus: T) = if (isPlus()) ifPlus else ifMinus

    /**
     * Useful in the case when the arguments evaluation can be expensive,
     * so we don't want to waste resources on both.
     * However, when [T] is a primitive type, the value will be boxed.
     * Choose [choose] wisely.
     */
    inline fun <T> choose(ifPlus: () -> T, ifMinus: () -> T) =
            if (isPlus()) ifPlus() else ifMinus()

    override fun toString() = char.toString()


    companion object {
        internal val LOG = Logger.getLogger(Strand::class.java)
    }

}

fun Int.toStrand() = if (this > 0) Strand.PLUS else Strand.MINUS

fun String.toStrand() = single().toStrand()

fun Char.toStrand() = when (this) {
    '+' -> Strand.PLUS
    '-' -> Strand.MINUS
    '.' -> Strand.PLUS // stranded
    else -> throw IllegalStateException("Unknown strand: $this")
}

enum class StrandFilter(private val char: Char) {
    BOTH('='), PLUS('+'), MINUS('-');

    fun accepts(strand: Strand) = when (this) {
        BOTH -> true
        PLUS -> strand.isPlus()
        MINUS -> strand.isMinus()
    }

    override fun toString() = char.toString()
}