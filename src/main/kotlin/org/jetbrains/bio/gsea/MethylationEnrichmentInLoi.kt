package org.jetbrains.bio.gsea

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.IntersectionNumberMetric
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream

// TODO: [make cmdline option] methylome coverage table is 1-based, convert to 0-based! loadInputRegionsAndMethylomeCovBackground:offsetIsOneBased

object MethylationEnrichmentInLoi {
    private val LOG = LoggerFactory.getLogger(MethylationEnrichmentInLoi::class.java)

    @JvmStatic
    fun main(args: Array<String>) {
        val metrics = listOf("overlap", "intersection")
        with(OptionParser()) {
            // Chrom.size
            acceptsAll(
                listOf("cs", "chrom.sizes"),
                "Chromosome sizes path, can be downloaded at" +
                        " http://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes"
            ).withRequiredArg()
                .required()
                .withValuesConvertedBy(PathConverter.exists())


            // Chromosome additional mapping
            acceptsAll(
                listOf("chrmap"),
                "Comma separated list of `chrInput:chrGenome` pairs, where `chrInput` is chromosome name in " +
                        "input file and `chrGenome` chromosome name in *.chrom.sizes file."
            ).withRequiredArg()
                .withValuesSeparatedBy(",")

            acceptsAll(
                listOf("b", "background"),
                "Covered methylome positions to use as random background. Strand is ignored. Must " +
                        "intersect with each of given input region. File should be TAB separated, 1st column is chromosome name, " +
                        "second column is 1-based offset. E.g. CpG offsets from DMC or proportions table."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            acceptsAll(
                listOf("genome-masked"),
                "Path to masked genome area file: Input regions intersecting masked area are skipped, masked area" +
                        " is also excluded from background. File is TAB-separated BED with at least 3 fields. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))

            acceptsAll(
                listOf("genome-allowed"),
                "Path to allowed genome area: Input regions/Background not-intersecting allowed area is skipped." +
                        " File is TAB-separated BED with at least 3 fields. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))

            acceptsAll(
                listOf("r", "regions"),
                "Regions to test file. File is TAB-separated BED with at least 3 fields. Strand is ignored." +
                        " E.g. regions could be DMRs."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))
                .required()

            acceptsAll(
                listOf("l", "loi"),
                "Loci of interest (LOI) to check enrichment. Could be single *.bed file or folder with *.bed" +
                        " files. Strand is ignored. E.g. loi could be CGIs file."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            acceptsAll(
                listOf("loi-filter"),
                "Filter LOI files ending with requested suffix string (not regexp)"
            )
                .withRequiredArg()

            acceptsAll(
                listOf("loi-intersect"),
                "Intersect LOI with given range 'chromosome:start-end'" +
                        " before over-representation test. 'start' offset 0-based, 'end' offset exclusive, i.e" +
                        " [start, end)"
            )
                .withRequiredArg()

            acceptsAll(
                listOf("metric"),
                "Metric is fun(a,b):Int. Supported values are: ${metrics.joinToString()}."
            )
                .withRequiredArg()
                .defaultsTo(metrics[0])

            // Output
            acceptsAll(
                listOf("o", "basename"),
                "Output files paths will start with given path prefix. Prefix could be relative" +
                        " or absolute path or it's part including folder name ending with '/'."
            )
                .withRequiredArg()
                .defaultsTo("regions_in_loi_regions_enrichment_")

            acceptsAll(
                listOf("h1"),
                "Alternative hypothesis: Input regions abundance in LOI is greater/less/different(two-sided) compared " +
                        "to simulated regions with similar lengths."
            )
                .withRequiredArg().ofType(PermutationAltHypothesis::class.java)
                .withValuesConvertedBy(PermutationAltHypothesis.converter())
                .defaultsTo(PermutationAltHypothesis.GREATER)

            acceptsAll(listOf("s", "simulations"), "Sampled regions sets number")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(100_000)

            acceptsAll(
                listOf("chunk-size"), "To reduce mem usage, simulation is done in several steps with" +
                        " provided simulations number per chunk. Use 0 to sample all sets at once."
            )
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(50_000)

            acceptsAll(listOf("retries"), "Regions shuffling max retries. Used because is not always possible to" +
                    " shuffle non-interested regions of given size.")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(1_000)


            accepts(
                "a-loi",
                "Metric is applied to (a,b) where by default 'a' is input/simulated regions regions set," +
                        " 'b' is LOI to check set. This option swaps a and b."
            )

            accepts(
                "replace",
                "Sample with replacement, i.e. sampled regions could intersect each other. By default w/o replacement"
            )

            acceptsAll(listOf("a-flanked"), "Flank 'a' ranges at both sides (non-negative dist in bp)")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(0)

            acceptsAll(listOf("parallelism"), "parallelism level")
                .withRequiredArg()
                .ofType(Int::class.java)

            accepts(
                "detailed",
                "Generated detailed stats report with metric value for each shuffle"
            )

            acceptsAll(
                listOf("m", "merge"),
                "Merge overlapped locations (e.g regions, loi) while loading files if applicable. Will speedup computations."
            )

            // Logging level:
            acceptsAll(listOf("d", "debug"), "Print all the debug info")

            val tool = BioinfToolsCLA.Tools.METHYLATION_ENRICHMENT_IN_LOI
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                LOG.info("CHROM.SIZES: $chromSizesPath")

                val chrMapStr = options.valuesOf("chrmap").filterIsInstance<String>()
                LOG.info("CHROM MAPPING: $chrMapStr")

                // val twoBitPath = options.valueOf("reference") as Path
                // LOG.info("GENOME REFERENCE: $twoBitPath")

                val genome = Genome.get(
                    chromSizesPath,
                    chrAltName2CanonicalMapping = parseChrNamesMapping(chrMapStr)
                )
                LOG.info("GENOME BUILD: ${genome.build}")

                val inputRegions = options.valueOf("regions") as Path // e.g. DMRs
                LOG.info("INPUT REGIONS: $inputRegions")

                val genomeMaskedAreaPath = options.valueOf("genome-masked") as Path?
                LOG.info("GENOME MASKED AREA: $genomeMaskedAreaPath")

                val genomeAllowedAreaPath = options.valueOf("genome-allowed") as Path?
                LOG.info("GENOME ALLOWED AREA: $genomeAllowedAreaPath")

                val methylomeBackground = options.valueOf("background") as Path
                LOG.info("METHYLOME BACKGROUND: $methylomeBackground")

                val loiFolderPath = options.valueOf("loi") as Path
                LOG.info("LOI TO TEST: $loiFolderPath")

                val loiNameSuffix = options.valueOf("loi-filter") as String?
                LOG.info("LOI FNAME SUFFIX: ${loiNameSuffix ?: "N/A"}")

                val loiLocationBasedFilter = (options.valueOf("loi-intersect") as String?)?.let {
                    EnrichmentInLoi.parseLocation(
                        it, genome, chromSizesPath
                    )
                }
                LOG.info("INTERSECT LOI WITH RANGE: ${loiLocationBasedFilter ?: "N/A"}")

                val baseName = options.valueOf("basename") as String
                val outputBaseName = FileSystems.getDefault().getPath(baseName).normalize().toAbsolutePath()
                LOG.info("OUTPUT_BASENAME: $outputBaseName")

                val hypAlt = options.valueOf("h1") as PermutationAltHypothesis
                LOG.info("Alt Hypothesis: $hypAlt")

                val simulationsNumber = options.valueOf("simulations") as Int
                LOG.info("SIMULATIONS: $simulationsNumber")

                val chunkSize = options.valueOf("chunk-size") as Int
                LOG.info("SIMULATIONS CHUNK SIZE: $chunkSize")

                val retries = options.valueOf("retries") as Int
                LOG.info("MAX RETRIES: $retries")

                val parallelism = options.valueOf("parallelism") as Int?
                configureParallelism(parallelism)
                LOG.info("THREADS: $parallelism (use #threads=${parallelismLevel()})")

                val detailedReport = options.has("detailed")
                LOG.info("DETAILED REPORT FLAG: $detailedReport")

                val samplingWithReplacement = options.has("replace")
                LOG.info("SAMPLE WITH REPLACEMENT: $samplingWithReplacement")

                val aSetIsRegions = !options.has("a-loi")
                if (aSetIsRegions) {
                    LOG.info("METRIC(a,b): a=INPUT/SIMULATED REGIONs, b=LOI")
                } else {
                    LOG.info("METRIC(a,b): a=LOI, b=INPUT/SIMULATED REGIONs")
                }

                val aSetFlankedBothSides = options.valueOf("a-flanked") as Int
                LOG.info("A FLANKED: $aSetFlankedBothSides")

                val mergeOverlapped = options.has("merge")
                LOG.info("MERGE OVERLAPPED: $mergeOverlapped")

                val metricStr = options.valueOf("metric") as String
                val metric = when (metricStr) {
                    "overlap" -> OverlapNumberMetric(aSetFlankedBothSides)
                    "intersection" -> IntersectionNumberMetric(aSetFlankedBothSides)
                    else -> throw IllegalArgumentException(
                        "Unsupported metric '${metricStr}', use one of: ${metrics.joinToString()}"
                    )
                }
                LOG.info("METRIC: ${metric.column}")

                doCalculations(
                    inputRegions, methylomeBackground, loiFolderPath, genome,
                    simulationsNumber,
                    if (chunkSize == 0) simulationsNumber else chunkSize,
                    outputBaseName,
                    metric,
                    detailedReport, retries, hypAlt,
                    aSetIsRegions,
                    mergeOverlapped,
                    loiLocationBasedFilter,
                    loiNameSuffix,
                    genomeMaskedAreaPath, genomeAllowedAreaPath,
                    samplingWithReplacement = samplingWithReplacement
                )
            }
        }
    }

    fun doCalculations(
        inputRegions: Path,
        methylomeBackgroundPath: Path,
        loiFolderPath: Path,
        genome: Genome,
        simulationsNumber: Int,
        chunkSize: Int,
        outputBasename: Path,
        metric: RegionsMetric,
        detailed: Boolean,
        maxRetries: Int,
        hypAlt: PermutationAltHypothesis,
        aSetIsRegions: Boolean,
        mergeOverlapped: Boolean,
        loiLocationBasedFilter: Location?,
        loiNameSuffix: String?,
        genomeMaskedAreaPath: Path?,
        genomeAllowedAreaPath: Path?,
        samplingWithReplacement: Boolean
    ) {
        val reportPath = "${outputBasename}${metric.column}.tsv".toPath()
        val detailedReportFolder = if (detailed) "${outputBasename}${metric.column}_stats".toPath() else null

        val loiFilter = loiLocationBasedFilter?.let {
            LocationsSortedList.create(genome.toQuery(), listOf(loiLocationBasedFilter))
        }

        val filesStream = if (loiFolderPath.isDirectory) Files.list(loiFolderPath) else Stream.of(loiFolderPath)
        val loiLabel2RangesList: List<Pair<String, LocationsList<out RangesList>>> = EnrichmentInLoi.collectLoiFrom(
            filesStream, genome, mergeOverlapped, loiFilter, loiNameSuffix
        )

        require(loiLabel2RangesList.isNotEmpty()) {
            "No LOI files passed file suffix filter."
        }

        RegionShuffleStatsFromMethylomeCoverage(
            genome,
            simulationsNumber, chunkSize,
            maxRetries
        ).calcStatistics(
            inputRegions, methylomeBackgroundPath,
            loiLabel2RangesList,
            detailedReportFolder,
            metric = metric,
            hypAlt = hypAlt,
            aSetIsRegions = aSetIsRegions,
            mergeOverlapped = mergeOverlapped,
            intersectionFilter = loiFilter,
            genomeMaskedAreaPath = genomeMaskedAreaPath,
            genomeAllowedAreaPath = genomeAllowedAreaPath,
            mergeRegionsToBg = false, // N/A
            samplingWithReplacement = samplingWithReplacement
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }
}
