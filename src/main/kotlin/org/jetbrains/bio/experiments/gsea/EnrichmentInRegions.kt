package org.jetbrains.bio.experiments.gsea

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.intersection.IntersectionMetric
import org.jetbrains.bio.genome.containers.intersection.IntersectionMetric.Companion.acceptMetricArg
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Collectors

object EnrichmentInRegions {
    private val LOG = LoggerFactory.getLogger(EnrichmentInRegions::class.java)

    fun doCalculations(
        loiPath: Path,
        backGroundRegions: Path?,
        regionsFolderPath: Path,
        genome: Genome,
        simulationsNumber: Int,
        chunkSize: Int,
        outputBasename: Path,
        metric: IntersectionMetric,
        detailed: Boolean,
        maxRetries: Int,
        hypAlt: PermutationAltHypothesis,
        aSetIsLoi: Boolean,
        aSetFlankedBothSides: Int,
        mergeOverlapped: Boolean
    ) {
        val reportPath = "${outputBasename}${metric.column}.tsv".toPath()
        val detailedReportFolder = if (detailed) "${outputBasename}${metric.column}_stats".toPath() else null

        RegionShuffleStats(
            genome,
            loiPath, backGroundRegions,
            simulationsNumber, chunkSize,
            maxRetries
        ).calcStatistics(
            collectFilesFrom(regionsFolderPath, genome, mergeOverlapped),
            detailedReportFolder,
            metric = metric, hypAlt = hypAlt, aSetFlankedBothSides = aSetFlankedBothSides, aSetIsLoi = aSetIsLoi,
            mergeOverlapped = mergeOverlapped
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }

    private fun collectFilesFrom(
        basePath: Path,
        genome: Genome,
        mergeOverlapped: Boolean
    ): List<Pair<String, LocationsMergingList>> {

        return if (basePath.isDirectory) {
            Files.list(basePath).map { path ->
                val name = path.fileName.toString()
                val locations = RegionShuffleStats.readLocations(path, genome, mergeOverlapped)
                name to locations
            }.collect(Collectors.toList())
        } else {
            val name = basePath.fileName.toString()
            val locations = RegionShuffleStats.readLocations(basePath, genome, mergeOverlapped)
            listOf(name to locations)
        }
    }

    @JvmStatic
    fun main(args: Array<String>) {
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
                listOf("l", "loi"),
                "Loci of interest (loi) file path in TAB separated BED, BED3 or BED4 format"
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))
                .required()

            acceptsAll(
                listOf("b", "background"),
                "Covered methylome regions to use as random background"
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())

            acceptsAll(
                listOf("r", "regions"),
                "Regions *.bed file or folder with *.bed regions files"
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            // TODO: optional save sampledRegions
            // TODO: optional load sampledRegions from ...

            // Output
            acceptsAll(
                listOf("o", "basename"),
                "Output files paths will start with given path prefix. Prefix could be relative" +
                        " or absolute path or it's part."
            )
                .withRequiredArg()
                .defaultsTo("loi_regions_enrichment_")

            acceptsAll(
                listOf("h1"),
                "Alternative hypothesis: Loi abundance is greater/less/different(two-sided) compared " +
                        "to simulated regions with similar lengths in provided regions"
            )
                .withRequiredArg().ofType(PermutationAltHypothesis::class.java)
                .withValuesConvertedBy(PermutationAltHypothesis.converter())
                .defaultsTo(PermutationAltHypothesis.TWO_SIDED)

            acceptsAll(listOf("s", "simulations"), "Sampled regions sets number")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(1000000)

            acceptsAll(listOf("chunk-size"), "To reduce mem usage, simulation is done in several steps with" +
                    " provided simulations number per chunk. Use 0 to sample all sets at once.")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(50000)

            acceptsAll(listOf("retries"), "Region shuffling max retries")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(1000)

            acceptMetricArg()

            accepts(
                "a-regions",
                "Metric is applied to (a,b) where by default 'a' is loi/simulated loi set, 'b' is region to check set. This" +
                        " option switches a and b definitions."
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

            val tool = BioinfToolsCLA.Tools.LOCI_ENRICHMENT_IN_REGIONS
            parse(args, description = "Given LOI enrichment in regions files") { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                LOG.info("CHROM.SIZES: $chromSizesPath")

                val chrMapStr = options.valuesOf("chrmap") as List<String>
                LOG.info("CHROM MAPPING: $chrMapStr")

                // val twoBitPath = options.valueOf("reference") as Path
                // LOG.info("GENOME REFERENCE: $twoBitPath")

                val genome = Genome.get(
                    chromSizesPath,
                    chrAltName2CanonicalMapping = parseChrNamesMapping(chrMapStr)
                )
                LOG.info("GENOME BUILD: ${genome.build}")

                val srcLoci = options.valueOf("loi") as Path
                LOG.info("SOURCE_LOCI: $srcLoci")

                val backGroundRegions = options.valueOf("background") as Path?
                LOG.info("BACKGROUND: $backGroundRegions")

                val regionsFolderPath = options.valueOf("regions") as Path
                LOG.info("REGIONS_TO_TEST: $regionsFolderPath")

                val baseName = options.valueOf("basename") as String
                val outputBaseName = FileSystems.getDefault().getPath(baseName).normalize().toAbsolutePath()
                LOG.info("OUTPUT_BASENAME: $outputBaseName")

                val metricName = options.valueOf("metric") as String
                LOG.info("METRIC: $metricName")
                val metric = IntersectionMetric.parse(metricName)

                val hypAlt = options.valueOf("h1") as PermutationAltHypothesis
                LOG.info("Alt Hypothesis: $hypAlt")

                val simulationsNumber = options.valueOf("simulations") as Int
                LOG.info("SIMULATIONS: $simulationsNumber")

                val chunkSize = options.valueOf("chunk-size") as Int
                LOG.info("SIMULATIONS CHUNK SIZE: $chunkSize")

                val retries = options.valueOf("retries") as Int
                LOG.info("MAX_RETRIES: $retries")

                val parallelism = options.valueOf("parallelism") as Int?
                LOG.info("THREADS: $parallelism")
                configureParallelism(parallelism)

                val detailedReport = options.has("detailed")
                LOG.info("DETAILED_REPORT_FLAG: $detailedReport")

                val aSetIsLoi = !options.has("a-regions")
                LOG.info("METRIC(REGION, LOI/SIMULATED LOI): $aSetIsLoi")

                val aSetFlankedBothSides = options.valueOf("a-flanked") as Int
                LOG.info("A FLANKED: $aSetFlankedBothSides")

                val mergeOverlapped = options.has("merge")
                LOG.info("MERGE OVERLAPPED: $mergeOverlapped")

                doCalculations(
                    srcLoci, backGroundRegions, regionsFolderPath, genome,
                    simulationsNumber,
                    if (chunkSize == 0) simulationsNumber else chunkSize,
                    outputBaseName, metric, detailedReport, retries, hypAlt,
                    aSetIsLoi = aSetIsLoi,
                    aSetFlankedBothSides = aSetFlankedBothSides,
                    mergeOverlapped = mergeOverlapped
                )
            }
        }
    }
}

fun parseChrNamesMapping(chrMapStr: List<String>) = chrMapStr.map {
    val pair = it.split(":", limit = 2)
    require(pair.size == 2) {
        "Chrom mapping should be in: '\$chrInput:\$chrGenome' format, but was: $pair"
    }
    pair[0].trim() to pair[1].trim()
}.toMap()