package org.jetbrains.bio.gsea

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.util.PathConverter
import org.jetbrains.bio.util.contains
import org.jetbrains.bio.util.parse
import org.jetbrains.bio.util.toPath
import org.slf4j.LoggerFactory

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

            // Force background to have 0-based offsets
            accepts(
                "bg-zero",
                "Background offsets in 0-based format instead of 1-based."
            )

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
                listOf("limit-intersect"),
                "Intersect LOI & simulated results with given range 'chromosome:start-end'" +
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
                "Metric will be applied to (a,b) where 'a' is LOI to check set, 'b' is input/simulated regions regions set. " +
                        "Without this option LOI will be 'b'."
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

                val samplingWithReplacement = options.has("replace")
                LOG.info("SAMPLE WITH REPLACEMENT: $samplingWithReplacement")

                val zeroBasedBg = "bg-zero" in options
                LOG.info("0-BASED OFFSETS IN BACKGROUND: $zeroBasedBg")

                val sharedOpts = EnrichmentInLoi.processSharedOptions(options, metrics, LOG)
                doCalculations(sharedOpts, samplingWithReplacement, zeroBasedBg)
            }
        }
    }

    fun doCalculations(
        sharedOpts: EnrichmentInLoi.SharedOptions,
        samplingWithReplacement: Boolean,
        zeroBasedBg: Boolean
    ) {
        val loiFolderPath = sharedOpts.loiFolderPath
        val outputBasename = sharedOpts.outputBaseName
        val limitResultsToSpecificLocation = sharedOpts.limitResultsToSpecificLocation
        val gq = sharedOpts.genome.toQuery()

        val reportPath = "${outputBasename}${sharedOpts.metric.column}.tsv".toPath()
        val detailedReportFolder =
            if (sharedOpts.detailedReport) "${outputBasename}${sharedOpts.metric.column}_stats".toPath() else null

        val limitResultsToSpecificLocationFilter: LocationsMergingList? = limitResultsToSpecificLocation?.let {
            LocationsMergingList.create(gq, listOf(limitResultsToSpecificLocation))
        }
        val loiInfos: List<LoiInfo> = EnrichmentInLoi.loiLocationBaseFilterList(
            sharedOpts,
            gq,
            limitResultsToSpecificLocationFilter,
            loiFolderPath
        )

        RegionShuffleStatsFromMethylomeCoverage(
            sharedOpts.simulationsNumber,
            sharedOpts.getChunkSize(), sharedOpts.retries,
            zeroBasedBg = zeroBasedBg,
        ).calcStatistics(
            gq,
            sharedOpts.inputRegions, // DMRs
            sharedOpts.backgroundPath, // Methylome background
            loiInfos,
            detailedReportFolder,
            metric = sharedOpts.metric,
            hypAlt = sharedOpts.hypAlt,
            aSetIsRegions = sharedOpts.aSetIsRegions,
            mergeOverlapped = sharedOpts.mergeOverlapped,
            intersectionFilter = limitResultsToSpecificLocationFilter,
            genomeMaskedAreaPath = sharedOpts.genomeMaskedAreaPath,
            genomeAllowedAreaPath = sharedOpts.genomeAllowedAreaPath,
            mergeRegionsToBg = false, // N/A
            samplingWithReplacement = samplingWithReplacement
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }
}
