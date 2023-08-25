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

// TODO: loci / regions rename
object MethylationEnrichmentInRegions {
    private val LOG = LoggerFactory.getLogger(MethylationEnrichmentInRegions::class.java)

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
                listOf("l", "loci"),
                "Methylome related loci file path in TAB separated BED, BED3 or BED4 format. Strand is ignored. E.g. DMRs."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))
                .required()

            acceptsAll(
                listOf("b", "background"),
                "Covered methylome positions to use as random background. Strand is ignored. Must " +
                        "intersect with each of given loci. File should be TAB separated, 1st column is chromosome name, " +
                        "second column is 1-based offset. E.g. CpG offsets from DMC or proportions table."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            acceptsAll(
                listOf("genome-masked"),
                "Mask genome regions: loci intersecting masked regions are skipped, masked regions removed from background." +
                        " File path in TAB separated BED with at least 3 fields. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))

            acceptsAll(
                listOf("genome-allowed"),
                "Allow only certain genome regions: input loci/background not-intersecting masked regions are skipped." +
                        " File path in TAB separated BED with at least 3 fields. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))


            acceptsAll(
                listOf("r", "regions"),
                "Regions *.bed file or folder with *.bed regions files. Strand is ignored. E.g. regions could be CGIs."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()


            acceptsAll(
                listOf("regions-filter"),
                "Filter regions files ending with requested suffix string (not regexp)"
            )
                .withRequiredArg()


            // TODO: optional save & load sampledRegions - to make results % reproducible w/o random effect

            acceptsAll(
                listOf("regions-intersect"),
                "Intersect regions with given range 'chromosome:start-end'" +
                        " before overrepresented regions test. 'start' offset 0-based, 'end' offset exclusive, i.e" +
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
                .defaultsTo("loi_regions_enrichment_")


            acceptsAll(
                listOf("h1"),
                "Alternative hypothesis: Loci abundance is greater/less/different(two-sided) compared " +
                        "to simulated regions with similar lengths in provided regions"
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

            acceptsAll(listOf("retries"), "Region shuffling max retries")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(1_000)


            accepts(
                "a-regions",
                "Metric is applied to (a,b) where by default 'a' is loci/simulated loci set, 'b' is region to check set. This" +
                        " option switches a and b definitions."
            )

            accepts(
                "replace",
                "Sample with replacement, i.e. sampled loci could intersect. By default w/o replacement"
            )

            acceptsAll(listOf("a-flanked"), "Flank 'a' regions at both sides (non-negative dist in bp)")
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
                "Merge overlapped locations while loading files if applicable. Will speedup computations."
            )

            // Logging level:
            acceptsAll(listOf("d", "debug"), "Print all the debug info")

            val tool = BioinfToolsCLA.Tools.METHYLATION_ENRICHMENT_IN_REGIONS
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

                val srcLoci = options.valueOf("loci") as Path // e.g. DMRs
                LOG.info("SOURCE_LOCI: $srcLoci")

                val genomeMaskedLociPath = options.valueOf("genome-masked") as Path?
                LOG.info("GENOME MASKED LOCI: $genomeMaskedLociPath")

                val genomeAllowedLociPath = options.valueOf("genome-allowed") as Path?
                LOG.info("GENOME ALLOWED LOCI: $genomeAllowedLociPath")

                val methylomeBackground = options.valueOf("background") as Path
                LOG.info("METHYLOME BACKGROUND: $methylomeBackground")

                val regionsFolderPath = options.valueOf("regions") as Path
                LOG.info("REGIONS_TO_TEST: $regionsFolderPath")
                val regionsNameSuffix = options.valueOf("regions-filter") as String?
                LOG.info("REGIONS_FNAME_SUFFIX: ${regionsNameSuffix ?: "N/A"}")

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

                val aSetIsLoi = !options.has("a-regions")
                if (aSetIsLoi) {
                    LOG.info("METRIC(a,b): a=REGION, b=LOCI/SIMULATED LOCI")
                } else {
                    LOG.info("METRIC(a,b): a=LOCI/SIMULATED LOCI, b=REGION")
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

                val regionsToTestFilterLocation = (options.valueOf("regions-intersect") as String?)?.let {
                    EnrichmentInRegions.parseLocation(
                        it, genome, chromSizesPath
                    )
                }
                LOG.info("INTERSECT REGIONS WITH RANGE: ${regionsToTestFilterLocation ?: "N/A"}")

                doCalculations(
                    srcLoci, methylomeBackground, regionsFolderPath, genome,
                    simulationsNumber,
                    if (chunkSize == 0) simulationsNumber else chunkSize,
                    outputBaseName,
                    metric,
                    detailedReport, retries, hypAlt,
                    aSetIsLoi,
                    mergeOverlapped,
                    regionsToTestFilterLocation,
                    regionsNameSuffix,
                    genomeMaskedLociPath, genomeAllowedLociPath,
                    samplingWithReplacement = samplingWithReplacement
                )
            }
        }
    }

    fun doCalculations(
        lociPath: Path,
        methylomeBackgroundPath: Path,
        regionsPath: Path,
        genome: Genome,
        simulationsNumber: Int,
        chunkSize: Int,
        outputBasename: Path,
        metric: RegionsMetric,
        detailed: Boolean,
        maxRetries: Int,
        hypAlt: PermutationAltHypothesis,
        aSetIsLoi: Boolean,
        mergeOverlapped: Boolean,
        regionsToTestFilterLocation: Location?,
        regionsNameSuffix: String?,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        samplingWithReplacement: Boolean
    ) {
        val reportPath = "${outputBasename}${metric.column}.tsv".toPath()
        val detailedReportFolder = if (detailed) "${outputBasename}${metric.column}_stats".toPath() else null

        val regionsToTestFilter = regionsToTestFilterLocation?.let {
            LocationsSortedList.create(
                genome.toQuery(),
                listOf(regionsToTestFilterLocation)
            )
        }

        val regionsAndRanges: List<Pair<String, LocationsList<out RangesList>>> =
            EnrichmentInRegions.collectRegionsFrom(
                if (regionsPath.isDirectory) Files.list(regionsPath) else Stream.of(regionsPath),
                genome, mergeOverlapped, regionsToTestFilter, regionsNameSuffix
            )

        require(regionsAndRanges.isNotEmpty()) {
            "No regions files passed file suffix filter."
        }

        RegionShuffleStatsFromMethylomeCoverage(
            genome,
            simulationsNumber, chunkSize,
            maxRetries
        ).calcStatistics(
            lociPath, methylomeBackgroundPath,
            regionsAndRanges,
            detailedReportFolder,
            metric = metric,
            hypAlt = hypAlt,
            aSetIsLoi = aSetIsLoi,
            mergeOverlapped = mergeOverlapped,
            intersectionFilter = regionsToTestFilter,
            genomeMaskedLociPath = genomeMaskedLociPath,
            genomeAllowedLociPath = genomeAllowedLociPath,
            addLoiToBg = false, // N/A
            samplingWithReplacement = samplingWithReplacement
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }
}
