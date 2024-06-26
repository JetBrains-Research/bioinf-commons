package org.jetbrains.bio.genome

import com.google.common.collect.Maps
import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.format.TwoBitReader
import org.jetbrains.bio.genome.format.TwoBitSequence
import org.jetbrains.bio.genome.format.TwoBitWriter
import org.jetbrains.bio.genome.format.writeAsFasta
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.IOException
import java.io.UncheckedIOException
import java.lang.ref.WeakReference
import java.nio.channels.ClosedByInterruptException
import java.nio.file.Path

/**
 * The genome.
 *
 * A [Genome] can be defined by [GenomeAnnotationsConfig] or manually.
 *
 * If your API is capable of processing a subset of chromosomes, consider
 * using [GenomeQuery] to indicate that.
 *
 * @author Sergei Lebedev
 */

class Genome private constructor(
        /** Unique build identifier, e.g. `"mm9"`. */
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
        val gapsPath: Path?,
        private val twoBitPath: Path?,
        private val fastaPath: Path?,
        private val genesGTFPath: Path?,
        genesDescriptionsPath: Path?,
        /**
         * Additional optional chromosome mapping for [chromosomeNamesMap], e.g. when genome is configured
         * from *.chrom.sizes file and not loaded from annotations.yaml.
         * E.g. it helps open locations created by UCSC hg19 genome in GRCh37 based genomes, where 'chrM' should be
         * converted to 'MT'
         */
        val chrAltName2CanonicalMapping: Map<String, String>
) {
    val chromSizesPath by lazy { ensureNotNull(chromSizesPath, "Chromosomes Sizes") }
    fun twoBitPath(downloadIfMissing: Boolean = true) =
            ensureNotNull(twoBitPath, "Genome *.2bit Sequence").also { twoBitPath ->
                if (downloadIfMissing) {
                    twoBitPath.checkOrRecalculate { output ->
                        val config = annotationsConfig
                        requireNotNull(config) {
                            "Cannot save Genome Sequence to $twoBitPath. Annotations information for $build isn't available."
                        }
                        val sUrl = config.sequenceUrl
                        when {
                            sUrl.endsWith(".2bit") -> sUrl.downloadTo(output.path)

                            else -> {
                                val suffix =
                                        listOf(".fa", ".fa.gz", ".fasta", ".fasta.gz").firstOrNull { sUrl.endsWith(it) }
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


    fun fastaPath(generateIfNeeded: Boolean = true) =
            ensureNotNull(fastaPath, "Genome *.fa Sequence").also { fastaPath ->
                if (generateIfNeeded) {
                    fastaPath.checkOrRecalculate { output ->
                        writeAsFasta(output.path)
                    }
                }
            }

    /**
     * Required for Star, RSEM
     */
    val genesNormalizedGtfPath: Path by lazy {
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
    fun genesGtfPath(downloadIfMissing: Boolean = true) =
            ensureNotNull(genesGTFPath, "Genes GTF Annotations").also { genesGTFPath ->
                if (downloadIfMissing) {
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

    val chromSizesMap: LinkedHashMap<String, Int> by lazy(chromSizesMapLambda)

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
     *     + mapping from genome [Genome.chrAltName2CanonicalMapping] field
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
            map[altName.lowercase()] = chr
            map[name] = chr
            map[name.lowercase()] = chr
        }
        chrAltName2CanonicalMapping.forEach { (altName, canonicalName) ->
            val chr = map[canonicalName]
            requireNotNull(chr) {
                "Unknown chromosome '$canonicalName' in genome[$build, chrom.sizes: $chromSizesPath], avaialbe chromosomes:" +
                        "${chromosomes.map { it.name }}"
            }
            map[altName] = chr
        }
        map
    }

    /**
     * Conventional to alternative names mapping
     *
     * Currently support the following names:
     *
     *     chr21
     *     CHRX
     *     21 or X
     *     + mapping from genome [Genome.chrAltName2CanonicalMapping] field
     * @author Roman Cherniatchik
     */
    val chromosomeNamesToAltNamesMap: Map<String, List<String>> by lazy {
        val map = hashMapOf<String, MutableList<String>>()
        chromosomeNamesMap.entries.forEach { (name, chr) ->
            val chrName = chr.name
            if (chrName != name) {
                val names = map.getOrPut(chrName) { arrayListOf() }
                names.add(name)
            }
        }
        map
    }
    val transcripts: Collection<Transcript> by lazy { Transcripts.all(this).values() }

    val genes: List<Gene> by lazy { groupTranscripts(transcripts) }

    fun presentableName(): String {
        val additional = StringBuffer()

        // we omit `build` and `ucsc_alias` from `names` since they are mentioned separately
        annotationsConfig?.names
                ?.filter { it != build && it != annotationsConfig.ucscAlias }
                ?.joinToString()
                ?.let { additional.append(it) }

        annotationsConfig?.ucscAlias?.let {
            if (it == build) return@let
            if (additional.isNotEmpty()) {
                additional.append(", ")
            }
            additional.append("UCSC alias: $it")
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

    private fun <T : Any> ensureNotNull(value: T?, tag: String): T {
        requireNotNull(value) {
            "$tag information was requested but it isn't available in '$build' genome." +
                    " This could happen if you are using customized genome which doesn't configure path to" +
                    " the information which is accessed later in your program."

        }
        return value
    }

    companion object {
        const val TEST_ORGANISM_BUILD = "to1"

        internal val LOG = LoggerFactory.getLogger(Genome::class.java)

        /** Use cache to avoid extra chrom.sizes loading. */
        private val CACHE = Maps.newConcurrentMap<String, Genome>()

        /**
         * Get or init default genome, which downloads all missing files to [Configuration.genomesPath]
         * See [Genome] constructor
         */
        operator fun get(build: String) = getCustomised(
                build, build,
                customized = false,
                forceUpdateCache = false,
                annConfigModifier = null
        )

        fun getCustomised(
                build: String,
                parentBuild: String,
                customized: Boolean,
                forceUpdateCache: Boolean = false,
                annConfigModifier: ((GenomeAnnotationsConfig) -> GenomeAnnotationsConfig)?
        ) =
                getOrAdd(build, customized, forceUpdateCache) {

                    val annCfg: GenomeAnnotationsConfig = when (parentBuild) {
                        // For tests
                        TEST_ORGANISM_BUILD -> GenomeAnnotationsConfig(
                                "Test Organism", TEST_ORGANISM_BUILD, listOf("to1"),
                                "test", "<n/a>", emptyMap(), false,
                            "<n/a>", "<n/a>", "<n/a>", "<n/a>", "<n/a>",
                                null, "<n/a>", null
                        )

                        else -> {
                            if (!AnnotationsConfigLoader.initialized) {
                                // Init with default settings only if not initialized by user before us
                                AnnotationsConfigLoader.init(Configuration.genomesPath / "annotations.yaml")
                            }
                            AnnotationsConfigLoader[parentBuild]
                        }
                    }
                    val dataPath = Configuration.genomesPath / parentBuild
                    val genesDescriptionsPath: Path = dataPath / "description.tsv"
                    val genesGTFPath: Path = dataPath / (when (parentBuild) {
                        TEST_ORGANISM_BUILD -> "genes.gtf.gz"
                        else -> annCfg.gtfUrl.substringAfterLast('/')
                    })

                    val annCfgUpdated = if (annConfigModifier == null) annCfg else annConfigModifier(annCfg)
                    Genome(
                            build,
                            annotationsConfig = annCfgUpdated,
                            // we don't expect test genome to save smth in test data dir
                            dataPath = if (parentBuild == TEST_ORGANISM_BUILD) null else dataPath,

                            chromSizesPath = dataPath / "$parentBuild.chrom.sizes",
                            chromSizesMapLambda = {
                                chromSizesFromPath(
                                        parentBuild,
                                        dataPath / "$parentBuild.chrom.sizes",
                                        annCfg
                                )
                            },
                            cpgIslandsPath = annCfgUpdated.cpgIslandsUrl?.let { dataPath / CpGIslands.ISLANDS_FILE_NAME },
                            cytobandsPath = annCfgUpdated.cytobandsUrl?.let { dataPath / CytoBands.FILE_NAME },
                            repeatsPath = dataPath / Repeats.FILE_NAME,
                            gapsPath = annCfgUpdated.gapsUrl?.let { dataPath / Gaps.FILE_NAME },
                            twoBitPath = dataPath / "$parentBuild.2bit",
                            fastaPath = dataPath / "$parentBuild.fa",
                            genesGTFPath = genesGTFPath,
                            genesDescriptionsPath = genesDescriptionsPath,
                            chrAltName2CanonicalMapping = annCfgUpdated.chrAltName2CanonicalMapping
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
                fastaPath: Path? = null,
                genesGTFPath: Path? = null,
                genesDescriptionsPath: Path? = null,
                chrAltName2CanonicalMapping: Map<String, String>? = null
        ) = getOrAdd(build, customized = true, forceUpdateCache = false) {
            val chromSizesDir = chromSizesPath.parent

            val chrMapping = HashMap<String, String>()
            annotationsConfig?.chrAltName2CanonicalMapping?.let { chrMapping.putAll(it) }
            chrAltName2CanonicalMapping?.let { chrMapping.putAll(it) }

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
                    fastaPath = fastaPath ?: (chromSizesDir / "$build.fa").let { if (it.exists) it else null },
                    genesGTFPath = genesGTFPath,
                    genesDescriptionsPath = genesDescriptionsPath,
                    chrAltName2CanonicalMapping = chrMapping
            )
        }

        /**
         * A custom genome with given chromosome sizes -- no paths required
         */
        operator fun get(
                build: String,
                chromSizesMap: LinkedHashMap<String, Int>
        ) = getOrAdd(build, customized = false, forceUpdateCache = false) {
            Genome(
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
                    fastaPath = null,
                    genesGTFPath = null,
                    genesDescriptionsPath = null,
                    chrAltName2CanonicalMapping = emptyMap()
            )
        }

        private fun getOrAdd(
                build: String,
                customized: Boolean,
                forceUpdateCache: Boolean,
                genomeProvider: () -> Genome
        ): Genome =
                if (!customized && !forceUpdateCache) {
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

                    if (forceUpdateCache) {
                        CACHE[build] = newGenome
                    }

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
                LOG.warn(
                        "Unexpected chrom sizes file name: $fileName, expected <build>.chrom.sizes. " +
                                "Detected build: $build"
                )
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
        /** Unique chromosome, could "19" or "MT" or be prefixed by `"chr"` (UCSC format), e.g. `"chr19", "chrM"`. */
        val name: String,
        /** Length defined in chrom.sizes file. */
        val length: Int
) {

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
                } catch (e: ClosedByInterruptException) {
                    // Thread requested sequence was interrupted, e.g. in JBR browser while rendering
                    throw e
                } catch (e: IOException) {
                    val errMsg = "Error loading $name from ${genome.twoBitPath(false)}"
                    LOG.error(errMsg, e) // Workaround: Cause of 'UncheckedIOException' not seen in stdout console, so add it to message
                    throw UncheckedIOException("Error loading $name from ${genome.twoBitPath(false)}, caused by ${e.javaClass.canonicalName}: msg=${e.message}", e)
                }

                sequenceRef = WeakReference(s)
            }

            return s!!
        }

    val range: Range get() = Range(0, length)

    val chromosomeRange: ChromosomeRange get() = ChromosomeRange(0, length, this)

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
                    val startOffset = centromeres.minOf { it.startOffset }
                    val endOffset = centromeres.maxOf { it.endOffset }
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
        internal val LOG = LoggerFactory.getLogger(Chromosome::class.java)

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

        val ADAPTER = object : TypeAdapter<Chromosome>() {
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
        internal val LOG = LoggerFactory.getLogger(Strand::class.java)

        fun fromBool(isPlusStrand: Boolean) = if (isPlusStrand) Strand.PLUS else MINUS
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