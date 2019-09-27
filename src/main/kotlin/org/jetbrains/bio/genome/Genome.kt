package org.jetbrains.bio.genome

import com.google.common.annotations.VisibleForTesting
import com.google.common.collect.Maps
import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.genome.sequence.TwoBitReader
import org.jetbrains.bio.genome.sequence.TwoBitSequence
import org.jetbrains.bio.genome.sequence.TwoBitWriter
import org.jetbrains.bio.util.*
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

        /** Annotations download urls and genome settings descriptor **/
        val annotationsConfig: GenomeAnnotationsConfig?,

        /** Absolute path to the genome data folder. */
        val dataPath: Path?,

        /** Absolute path to the genome chromosome sizes path,
         * if **not exists** will be downloaded automatically. */
        chromSizesPath: Path?,

        /**
         * Content of chrom.sizes file. Preserves original order during iteration.
         */
        chromSizesMapLambda: () -> LinkedHashMap<String, Int>,

        val cpgIslandsPath: Path?,
        val cytobandsPath: Path?,
        val repeatsPath: Path?,
        gapsPath: Path?,
        private val twoBitPath: Path?,
        private val genesGTFPath: Path?,
        genesDescriptionsPath: Path?
) {
    val chromSizesPath by lazy { ensureNotNull(chromSizesPath, "Chromosomes Sizes") }
    val gapsPath by lazy { ensureNotNull(gapsPath, "Gaps") }
    fun twoBitPath(downloadIfMissed: Boolean = true) =
            ensureNotNull(twoBitPath, "Genome *.2bit Sequence").also { twoBitPath ->
                if (downloadIfMissed) {
                    twoBitPath.checkOrRecalculate { output ->
                        val config = annotationsConfig
                        requireNotNull(config) {
                            "Cannot save Genome Sequence to $twoBitPath. Annotations information for $build isn't available."
                        }
                        val sUrl = config.sequenceUrl
                        when {
                            sUrl.endsWith(".2bit") -> sUrl.downloadTo(output.path)

                            else -> {
                                val suffix = listOf(".fa", ".fa.gz", ".fasta", ".fasta.gz").firstOrNull() { sUrl.endsWith(it) }
                                requireNotNull(suffix) { "Unsupported sequence type: $sUrl" }

                                val faPath = "${output.path}$suffix".toPath()
                                sUrl.downloadTo(faPath)
                                TwoBitWriter.convert(faPath, output.path)
                                faPath.delete()
                            }
                        }
                    }
                }
            }

    /**
     * Required for Star, RSEM
     */
    val genesNormalizedGtfPath: Path  by lazy {
        val gtfPath = genesGtfPath(false)
        val baseName = when {
            gtfPath.name.endsWith(".gtf.gz") -> {
                gtfPath.name.subSequence(0, gtfPath.name.length - ".gtf.gz".length)
            }
            gtfPath.extension == "gtf" -> gtfPath.stem
            else -> error("Unsupported GTF path: $gtfPath")
        }

        val normGtfPath = gtfPath.parent / "$baseName.norm.gtf"

        normGtfPath.checkOrRecalculate("GTF annotations") { (path) ->
            genesGtfPath(true).bufferedReader().use { reader ->
                Ensembl.convertGTF(this, reader, path)
            }
        }
        normGtfPath
    }

    val genesDescriptionsPath: Path by lazy { ensureNotNull(genesDescriptionsPath, "Gene Description") }

    /**
     * Ensure *.gtf file exists and download it if necessary
     */
    fun genesGtfPath(downloadIfMissed: Boolean = true) =
            ensureNotNull(genesGTFPath, "Genes GTF Annotations").also { genesGTFPath ->
                if (downloadIfMissed) {
                    genesGTFPath.checkOrRecalculate("Genes") { output ->
                        val config = annotationsConfig
                        requireNotNull(config) {
                            "Cannot save genes GTF to $genesGTFPath. Annotations information isn't available for $build."
                        }
                        config.gtfUrl.downloadTo(output.path)
                    }
                }
            }

    /** Species token, e.g. `"mm"`. */
    val species: String get() = build.takeWhile { !it.isDigit() }

    internal val chromSizesMap: LinkedHashMap<String, Int> by lazy(chromSizesMapLambda)

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
            val altName = if (name.startsWith("chr")) {
                name.substring(3)
            } else {
                "chr$name"
            }
            map[altName] = chr
            map[altName.toLowerCase()] = chr
            map[name] = chr
            map[name.toLowerCase()] = chr
        }
        val chrAltName2CanonicalMapping = annotationsConfig?.chrAltName2CanonicalMapping ?: emptyMap()
        chrAltName2CanonicalMapping.forEach {(altName,canonicalName) ->
            val chr = map[canonicalName]
            requireNotNull(chr) {
                "Unknown chromosome '$canonicalName' in genome[$build, chrom.sizes: $chromSizesPath]"
            }
            map[altName] = chr
        }
        map
    }

    val transcripts: Collection<Transcript> by lazy { Transcripts.all(this).values() }

    val genes: List<Gene> by lazy { groupTranscripts(transcripts) }

    fun presentableName(): String {
        val additional = StringBuffer()
        annotationsConfig?.alias?.let {
            if (it.isNotEmpty()) {
                additional.append(it)
            }
        }
        annotationsConfig?.description?.let {
            if (it.isNotEmpty()) {
                if (additional.isNotEmpty()) {
                    additional.append("; ")
                }
                additional.append(it)
            }
        }
        val additionalStr = if (additional.isEmpty()) "" else " [$additional]"
        return "${annotationsConfig?.species ?: "Unknown Species"}: $build$additionalStr"
    }

    private fun <T: Any> ensureNotNull(value:T?, tag: String): T {
        requireNotNull(value) {
            "$tag information was requested but it isn't available in '$build' genome." +
                    " This could happen if you are using customized genome which doesn't configure path to" +
                    " the information which is accessed later in your program."

        }
        return value
    }

    companion object {
        const val TEST_ORGANISM_BUILD = "to1"

        @VisibleForTesting
        internal val LOG = Logger.getLogger(Genome::class.java)

        /** Use cache to avoid extra chrom.sizes loading. */
        private val CACHE = Maps.newConcurrentMap<String, Genome>()

        /**
         * Get or init default genome, which downloads all missing files to [Configuration.genomesPath]
         * See [Genome] constructor
         */
        operator fun get(build: String) =
                getOrAdd(build, false) {
                    val dataPath = Configuration.genomesPath / build
                    val to = build == TEST_ORGANISM_BUILD

                    val annCfg: GenomeAnnotationsConfig = when {
                        to -> GenomeAnnotationsConfig(
                            "Test Organism", null, "test", "<n/a>", emptyMap(), false,
                            "<n/a>", "<n/a>", "<n/a>", "<n/a>", "<n/a>",
                            null, "<n/a>", null)

                        else -> {
                            if (!AnnotationsConfig.initialized) {
                                // Init with default settings only if not initialized by user before us
                                AnnotationsConfig.init(Configuration.genomesPath / "annotations.yaml")
                            }
                            AnnotationsConfig[build]
                        }
                    }
                    val genesDescriptionsPath: Path? = if (to) null else dataPath / "description.tsv"
                    val genesGTFPath: Path = dataPath / (when {
                        to -> "genes.gtf.gz"
                        else -> annCfg.gtfUrl.substringAfterLast('/')
                    })

                    Genome(
                        build,
                        annotationsConfig = annCfg,
                        // we don't expect test genome to save smth in test data dir
                        dataPath = if (to) null else dataPath,

                        chromSizesPath = dataPath / "$build.chrom.sizes",
                        chromSizesMapLambda = { chromSizesFromPath(build, dataPath / "$build.chrom.sizes", annCfg) },
                        cpgIslandsPath = annCfg.cpgIslandsUrl?.let { dataPath / CpGIslands.ISLANDS_FILE_NAME },
                        cytobandsPath = annCfg.cytobandsUrl?.let { dataPath / CytoBands.FILE_NAME },
                        repeatsPath = dataPath / Repeats.FILE_NAME,
                        gapsPath = dataPath / Gaps.FILE_NAME,
                        twoBitPath = dataPath / "$build.2bit",
                        genesGTFPath = genesGTFPath,
                        genesDescriptionsPath = genesDescriptionsPath
                    )
                }

        /**
         * Get or init customized genome
         *
         * See [Genome] constructor
         */
        operator fun get(
                chromSizesPath: Path,
                build: String = buildNameFrom(chromSizesPath),
                annotationsConfig: GenomeAnnotationsConfig? = null,
                dataPath: Path? = null,
                cpgIslandsPath: Path? = null,
                cytobandsPath: Path? = null,
                repeatsPath: Path? = null,
                gapsPath: Path? = null,
                twoBitPath: Path? = null,
                genesGTFPath: Path? = null,
                genesDescriptionsPath: Path? = null
        ) = getOrAdd(build, true) {
            val chromSizesDir = chromSizesPath.parent
            Genome(
                build,
                annotationsConfig = annotationsConfig,
                dataPath = dataPath,
                chromSizesPath = chromSizesPath,
                chromSizesMapLambda = { chromSizesFromPath(build, chromSizesPath, annotationsConfig) },
                cpgIslandsPath = cpgIslandsPath,
                cytobandsPath = cytobandsPath,
                repeatsPath = repeatsPath,
                gapsPath = gapsPath,
                twoBitPath = twoBitPath ?: (chromSizesDir / "$build.2bit").let { if (it.exists) it else null },
                genesGTFPath = genesGTFPath,
                genesDescriptionsPath = genesDescriptionsPath

            )
        }

        /**
         * A custom genome with given chromosome sizes -- no paths required
         */
        operator fun get(
                build: String,
                chromSizesMap: LinkedHashMap<String, Int>
        ) = getOrAdd(build, false) { Genome(
            build = build,
            annotationsConfig = null,
            dataPath = null,
            chromSizesPath = null,
            chromSizesMapLambda = { chromSizesMap },
            cpgIslandsPath = null,
            cytobandsPath = null,
            repeatsPath = null,
            gapsPath = null,
            twoBitPath = null,
            genesGTFPath = null,
            genesDescriptionsPath = null
        ) }

        private fun getOrAdd(build: String, customized: Boolean, genomeProvider: () -> Genome): Genome =
                if (!customized) {
                    // If genome getter not customized:
                    // * genome not initialized => return genome with default paths
                    // * already initialized => return cached genome even if cached genome is a customized one

                    // It is useful for CLI tools & genome string parsing like in DataConfig experiments
                    // We could init "hg19" with custom paths (on first access), then DataConfig static
                    // deserialization could load our custom genome for "hg19" instead of genome with default
                    // paths.
                    //
                    // P.S: We could remove this code and leave impl like for customized genome (see the other branch)
                    CACHE.computeIfAbsent(build) {
                        genomeProvider()
                    }
                } else {
                    // Create genome, it is cheap. Need further to compare with cached version
                    val newGenome = genomeProvider()

                    val cachedGenome = CACHE.computeIfAbsent(build) {
                        newGenome
                    }

                    // Assume 'build' is genome id, otherwise it could be misleading and we could get different
                    // [Genome] instances where we do not expect this.

                    // Compare that genome is same as cached one:
                    require(cachedGenome.chromSizesPath == newGenome.chromSizesPath) {
                        "Cannot load genome for ${newGenome.chromSizesPath}: genome '$build' already initialized" +
                                " with ${cachedGenome.chromSizesPath}"
                    }

                    // XXX: maybe implement equals (e.g. convert genome to dataclass) & require(genome == newInstance) ?

                    cachedGenome
                }

        private fun buildNameFrom(chromSizesPath: Path): String {
            val fileName = chromSizesPath.fileName.toString()

            if (!fileName.endsWith(".chrom.sizes")) {
                val build = fileName.substringBefore(".")
                LOG.warn("Unexpected chrom sizes file name: $fileName, expected <build>.chrom.sizes. " +
                        "Detected build: $build")
                return build
            }
            val build = fileName.substringBeforeLast(".chrom.sizes")
            LOG.debug("Chrom sizes name: $fileName. Detected build: $build")
            return build
        }

        private fun chromSizesFromPath(
                build: String, chromSizesPath: Path, annotationsConfig: GenomeAnnotationsConfig?
        ): LinkedHashMap<String, Int> {

            // Fetch chrom sizes path if not exists
            chromSizesPath.checkOrRecalculate { (p) ->
                requireNotNull(annotationsConfig) {
                    "Cannot save chromosomes sizes to $chromSizesPath. Annotations information isn't available" +
                            " for $build."
                }

                annotationsConfig.chromsizesUrl.downloadTo(p)
            }

            LOG.debug("Loading chrom.sizes $chromSizesPath")
            val map = LinkedHashMap<String, Int>()
            CSVFormat.TDF.parse(chromSizesPath.bufferedReader()).use { parser ->
                parser.records.forEach {
                    try {
                        map[it[0]] = it[1].toInt()
                    } catch (t: Throwable) {
                        LOG.error(
                            "Failed to parse chrom.sizes file: $chromSizesPath," +
                                    " line: ${it.joinToString("\t")}", t
                        )
                    }
                }
                LOG.debug("DONE Loading chrom.sizes $chromSizesPath")
                return map
            }
        }
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
        val length: Int) {

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
                    val twoBitPath = genome.twoBitPath()
                    val twoBitLength = TwoBitReader.length(twoBitPath, name)
                    check(length == twoBitLength) {
                        "Chromosome $name length differs in chrom.sizes($length) and 2bit file($twoBitLength)"
                    }
                    TwoBitReader.read(twoBitPath, name)
                } catch (e: IOException) {
                    throw UncheckedIOException(
                        "Error loading $name from ${genome.twoBitPath(false)}", e)
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