package org.jetbrains.bio.gsea

import joptsimple.OptionParser
import joptsimple.OptionSet
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.IntersectionNumberMetric
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Collectors
import java.util.stream.Stream

// TODO: optional save/load sampled regions to make 100% reproducible

/**
 * The EnrichmentInLoi class is responsible for performing calculations related to enrichment of regions of interest (LOI)
 * in a given set of input regions.
 */
object EnrichmentInLoi {
    private val LOG = LoggerFactory.getLogger(EnrichmentInLoi::class.java)

    fun doCalculations(
        opts: SharedSamplingOptions,
        enrichmentOpts: SharedEnrichmentOptions,
        mergeInputRegionsToBg: Boolean, // add input regions to background if they are missing there
        backgroundIsMethylome: Boolean,
        zeroBasedBackground: Boolean,
        bedBgFlnk: Int,
    ) {
        val loiFolderPath = enrichmentOpts.loiFolderPath
        val outputBasename = opts.outputBaseName
        val truncateRangesToSpecificLocation = opts.truncateRangesToSpecificLocation

        val reportPath = "${outputBasename}${enrichmentOpts.metric.column}.tsv".toPath()
        val detailedReportFolder =
            if (enrichmentOpts.detailedReport) "${outputBasename}${enrichmentOpts.metric.column}_stats".toPath() else null

        val gq = opts.genome.toQuery()

        val truncateRangesToSpecificLocationFilter: LocationsMergingList? = truncateRangesToSpecificLocation?.let {
            LocationsMergingList.create(gq, listOf(truncateRangesToSpecificLocation))
        }
        val genomeAllowedAreaFilter = opts.genomeAllowedAreaPath?.let {
            RegionShuffleStats.readGenomeAreaFilter(it, gq)
        }
        val genomeMaskedAreaFilter = opts.genomeMaskedAreaPath?.let {
            RegionShuffleStats.readGenomeAreaFilter(it, gq)
        }

        val loiInfos: List<LoiInfo> = loiLocationBaseFilterList(
            opts,
            enrichmentOpts.loiNameSuffix,
            gq,
            genomeAllowedAreaFilter,
            genomeMaskedAreaFilter,
            truncateRangesToSpecificLocationFilter,
            loiFolderPath
        )

        RegionShuffleStatsFromBedCoverage(
            opts.simulationsNumber,
            opts.getChunkSize(), opts.setRetries, opts.regionRetries
        ).calcStatistics(
            gq,
            opts.inputRegions,
            enrichmentOpts.backgroundPath,
            loiInfos,
            detailedReportFolder,
            metric = enrichmentOpts.metric,
            hypAlt = enrichmentOpts.hypAlt,
            aSetIsRegions = enrichmentOpts.aSetIsRegions,
            mergeOverlapped = opts.mergeOverlapped,
            truncateFilter = truncateRangesToSpecificLocationFilter,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter,
            mergeInputRegionsToBg = mergeInputRegionsToBg,
            samplingWithReplacement = opts.samplingWithReplacement,
            backgroundIsMethylome = backgroundIsMethylome,
            zeroBasedBackground = zeroBasedBackground,
            bedBgFlnk = bedBgFlnk,
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }

    fun loiLocationBaseFilterList(
        opt: SharedSamplingOptions,
        loiNameSuffix: String?,
        gq: GenomeQuery,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
        truncateFilter: LocationsMergingList?,
        loiFolderPath: Path
    ): List<LoiInfo> {
        val filesStream = when {
            loiFolderPath.isDirectory -> Files.list(loiFolderPath)
            else -> Stream.of(loiFolderPath)
        }
        val loiInfosFiltered: List<LoiInfo> = collectLoiFrom(
            filesStream, gq, opt.mergeOverlapped,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter,
            truncateFilter = truncateFilter,
            loiNameSuffix = loiNameSuffix
        )

        require(loiInfosFiltered.isNotEmpty()) {
            "No LOI files passed file suffix filter."
        }
        return loiInfosFiltered
    }

    fun processInputRegions(
        inputRegionsPath: Path,
        gq: GenomeQuery,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
    ): List<Location> {
        val (inputRegions, recordsNumber) = RegionShuffleStats.readLocationsIgnoringStrand(inputRegionsPath, gq)
        LOG.info("Loaded input regions: ${inputRegions.size.formatLongNumber()} regions of ${recordsNumber.formatLongNumber()} lines.")

        return filterRegions(
            "Input Regions", inputRegions,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter
        )
    }

    fun filterRegions(
        regionTypeLabel: String,
        regions: List<Location>,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?
    ): List<Location> {
        // Here no need in intersecting input with `SharedSamplingOptions.truncateRangesToSpecificLocation`

        LOG.info("$regionTypeLabel: ${regions.size.formatLongNumber()} regions")

        val allowed = if (genomeAllowedAreaFilter != null) {
            LOG.info("$regionTypeLabel: Applying allowed area filter...")

            val filtered = regions.filter { genomeAllowedAreaFilter.intersects(it) }
            LOG.info("$regionTypeLabel: [DONE] Allowed-filtered regions: ${filtered.size.formatLongNumber()} regions of ${regions.size.formatLongNumber()}")

            filtered
        } else {
            regions
        }

        val allowedAndMasked = if (genomeMaskedAreaFilter != null) {
            LOG.info("$regionTypeLabel: Applying masked area filter...")

            // In bedtools shuffle rejects sampled regions if they intersect masked area => let's exclude input
            // regions intersecting masked area for the consistency with sampling.
            // Let's keep the background unchanged and verify intersection with masked area later after the candidate
            // region is sampled
            val filtered = allowed.filterNot { genomeMaskedAreaFilter.intersects(it) }

            LOG.info("$regionTypeLabel: [DONE] w/o masked: ${filtered.size.formatLongNumber()} regions of ${allowed.size.formatLongNumber()}")

            filtered
        } else {
            allowed
        }

        return allowedAndMasked
    }

    fun collectLoiFrom(
        filesStream: Stream<Path>,
        gq: GenomeQuery,
        mergeOverlapped: Boolean,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
        truncateFilter: LocationsMergingList?,
        loiNameSuffix: String?
    ): List<LoiInfo> = filesStream.filter { path ->
        loiNameSuffix == null || path.name.endsWith(loiNameSuffix)
    }.map { path ->
        val name = path.fileName.toString()

        val (locList, nLoadedLoi, nRecords) = RegionShuffleStats.readLocationsIgnoringStrand(
            path, gq, mergeOverlapped
        )

        val filteredLoi = filterRegions(
            "LOI [${path.fileName}]", locList.asLocationSequence().toList(), 
            genomeAllowedAreaFilter, genomeMaskedAreaFilter
        )

        var filteredLocList: LocationsList<out RangesList> = when {
            mergeOverlapped -> LocationsMergingList.create(gq, filteredLoi)
            else -> LocationsSortedList.create(gq, filteredLoi)
        }
        
        if (truncateFilter != null) {
            filteredLocList = filteredLocList.intersectRanges(truncateFilter)
        }
        // LOG.warn("DEBUG REGIONS=${locList.asLocationSequence().joinToString { it.toString() }}")

        LOG.info("LOI [$name] (all filters applied): ${filteredLocList.size.formatLongNumber()} regions of ${nLoadedLoi.formatLongNumber()} loaded region")

        LoiInfo(name, filteredLocList, processedLoiNumber = nLoadedLoi, recordsNumber = nRecords)
    }.collect(Collectors.toList())

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
                "Covered methylome BED file to use as random background. Strand is ignored. Must include input regions."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())

            accepts(
                "bg-methylome",
                "Treat background path as methylome. In this case background will " +
                        " be a set of merged regions, where each cytosine offset will be converted to a region" +
                        " [offset - flnk, offset + flnk]. `flnk` is flanking region radius and configured by" +
                        " `methylome-flnk` option."
            )
            // Force methylome to have 0-based offsets
            accepts(
                "methylome-zero",
                "If background is methylome and methylome offsets are in 0-based format instead of 1-based."
            )
            acceptsAll(listOf("methylome-flnk"), "Flanking region radius for methylome background.")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(50)

            accepts(
                "add-regions-to-bg",
                "Merge input regions into background before sampling if background doesn't cover all input regions."
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
                "Truncate LOI & simulated results by range 'chromosome:start-end'" +
                        " before over-representation test. 'start' offset 0-based, 'end' offset exclusive, i.e" +
                        " [start, end). Regions will be simulated from background and truncated before metric " +
                        " calculation. E.g. could be used to estimate over-representation in context of specific " +
                        " region, e.g. HLA locus. It is not the same as limit background only to HLA locus and " +
                        " sample input regions that also overlap the locus."
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
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(100_000)

            acceptsAll(
                listOf("chunk-size"), "To reduce mem usage, simulation is done in several steps with" +
                        " provided simulations number per chunk. Use 0 to sample all sets at once."
            )
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(50_000)

            acceptsAll(
                listOf("set_retries"), "Regions set sampling max retries. Used because is not always possible to" +
                        " sample the whole set.")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(1_000)

            acceptsAll(
                listOf("region_retries"), "Individual region sampling max retries. Used because is not always" +
                        " possible to sample region with given constraints on length.")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(10_000)

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
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(0)

            acceptsAll(listOf("parallelism"), "parallelism level")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)

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

            val tool = BioinfToolsCLA.Tools.ENRICHMENT_IN_LOI
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val mergeRegionsToBg = options.has("add-regions-to-bg")
                LOG.info("MERGE REGIONS INTO BACKGROUND: $mergeRegionsToBg")

                val sharedOpts = processSharedSamplingOptions(options, LOG)
                val sharedEnrichmentOpts = processSharedEnrichmentOptions(options, metrics, LOG)

                val backgroundIsMethylome = "bg-methylome" in options
                LOG.info("BACKGROUND IS METHYLOME: $backgroundIsMethylome")

                val bedBgFlnk: Int
                val zeroBasedBackground: Boolean
                if (backgroundIsMethylome) {
                    zeroBasedBackground = "methylome-zero" in options
                    LOG.info("0-BASED OFFSETS IN METHYLOME: $zeroBasedBackground")

                    bedBgFlnk = options.valueOf("methylome-flnk") as Int
                    LOG.info("BED BACKGROUND FLANKING RADIUS: ${bedBgFlnk.formatLongNumber()}")
                } else {
                    bedBgFlnk = 0
                    zeroBasedBackground = true
                }
                doCalculations(
                    sharedOpts, sharedEnrichmentOpts, mergeRegionsToBg,
                    backgroundIsMethylome = backgroundIsMethylome,
                    zeroBasedBackground = zeroBasedBackground,
                    bedBgFlnk = bedBgFlnk,
                )
            }
        }
    }

    fun parseLocation(
        intersectionFilterStr: String,
        genome: Genome,
        chromSizesPath: Path
    ): Location {
        var intersectionFilterStr1 = intersectionFilterStr
        intersectionFilterStr1 = intersectionFilterStr1
            // whitespaces:
            .replace(" ", "")
            .replace("\t", "")
            .replace("\n", "")
            .replace("\r", "")
            // number format separators
            .replace(".", "")
            .replace(",", "")

        val chunks = intersectionFilterStr1.split(':')
        require(chunks.size == 2) {
            "Intersection range should be in format: 'chromosome:start-end' but was: $intersectionFilterStr1"
        }
        val chrName = chunks[0]
        val offsetChunks = chunks[1].split('-')
        require(offsetChunks.size == 2) {
            "Intersection range should be in format: 'chromosome:start-end' but was: $intersectionFilterStr1"
        }
        val startOffset = offsetChunks[0].toIntOrNull()
        val endOffset = offsetChunks[1].toIntOrNull()
        require(startOffset != null && endOffset != null) {
            "Intersection range should be in format: 'chromosome:start-end' but was: $intersectionFilterStr1"
        }
        val gq = genome.toQuery()
        val chromosome = gq[chrName]
        requireNotNull(chromosome) {
            "Cannot find chromosome '$chrName' in genome $genome ($chromSizesPath)"
        }
        require(startOffset in chromosome.range) {
            "Start offset $startOffset is out of chromosome '${chrName}' ${chromosome.range} range"
        }
        require(endOffset in chromosome.range) {
            "End offset $endOffset is out of chromosome '${chrName}' ${chromosome.range} range"
        }
        return Location(startOffset, endOffset, chromosome)
    }

    class SharedSamplingOptions(
        val genome: Genome,
        val outputBaseName: Path,
        private val chunkSize: Int,
        val simulationsNumber: Int,
        val setRetries: Int,
        val regionRetries: Int,
        val truncateRangesToSpecificLocation: Location?,
        val inputRegions: Path,
        val mergeOverlapped: Boolean,
        val genomeMaskedAreaPath: Path?,
        val genomeAllowedAreaPath: Path?,
        val parallelism: Int?,
        val samplingWithReplacement: Boolean
    ) {
        fun getChunkSize(): Int = if (chunkSize == 0) simulationsNumber else chunkSize
    }
    class SharedEnrichmentOptions(
        val aSetIsRegions: Boolean,
        val metric: RegionsMetric,
        val detailedReport: Boolean,
        val hypAlt: PermutationAltHypothesis,
        val backgroundPath: Path?,
        val loiFolderPath: Path,
        val loiNameSuffix: String?,
    )

    fun processSharedSamplingOptions(
        options: OptionSet,
        logger: Logger
    ): SharedSamplingOptions {
        val chromSizesPath = options.valueOf("chrom.sizes") as Path
        logger.info("CHROM.SIZES: $chromSizesPath")

        val chrMapStr = options.valuesOf("chrmap").filterIsInstance<String>()
        logger.info("CHROM MAPPING: $chrMapStr")

        // val twoBitPath = options.valueOf("reference") as Path
        // LOG.info("GENOME REFERENCE: $twoBitPath")

        val genome = Genome.get(
            chromSizesPath,
            chrAltName2CanonicalMapping = parseChrNamesMapping(chrMapStr)
        )
        logger.info("GENOME BUILD: ${genome.build}")

        val baseName = options.valueOf("basename") as String
        val outputBaseName = FileSystems.getDefault().getPath(baseName).normalize().toAbsolutePath()
        logger.info("OUTPUT_BASENAME: $outputBaseName")

        val simulationsNumber = options.valueOf("simulations") as Int
        logger.info("SIMULATIONS: ${simulationsNumber.formatLongNumber()}")

        val chunkSize = options.valueOf("chunk-size") as Int
        logger.info("SIMULATIONS CHUNK SIZE: ${chunkSize.formatLongNumber()}")

        val setRetries = options.valueOf("set_retries") as Int
        logger.info("SET MAX RETRIES: ${setRetries.formatLongNumber()}")

        val regionRetries = options.valueOf("region_retries") as Int
        logger.info("REGION MAX RETRIES: ${regionRetries.formatLongNumber()}")

        val limitResultsToSpecificLocation = (options.valueOf("limit-intersect") as String?)?.let {
            parseLocation(
                it, genome, chromSizesPath
            )
        }
        // TODO: LOI, SAMPLED & INPUT , e.g. lightweitgh allowed regions ?
        // TODO: or just LOI ?
        logger.info("INTERSECT LOI & SAMPLED REGIONS WITH RANGE: ${limitResultsToSpecificLocation ?: "N/A"}")

        val inputRegions = options.valueOf("regions") as Path
        logger.info("INPUT REGIONS: $inputRegions")

        val mergeOverlapped = options.has("merge")
        logger.info("MERGE OVERLAPPED: $mergeOverlapped")

        val genomeMaskedAreaPath = options.valueOf("genome-masked") as Path?
        logger.info("GENOME MASKED AREA: $genomeMaskedAreaPath")
        val genomeAllowedAreaPath = options.valueOf("genome-allowed") as Path?
        logger.info("GENOME ALLOWED AREA: $genomeAllowedAreaPath")

        val parallelism = options.valueOf("parallelism") as Int?
        configureParallelism(parallelism)
        logger.info("THREADS: $parallelism (use #threads=${parallelismLevel()})")

        val samplingWithReplacement = options.has("replace")
        LOG.info("SAMPLE WITH REPLACEMENT: $samplingWithReplacement")

        return SharedSamplingOptions(
            genome = genome,
            outputBaseName = outputBaseName,
            chunkSize = chunkSize,
            simulationsNumber = simulationsNumber,
            setRetries = setRetries,
            regionRetries = regionRetries,
            truncateRangesToSpecificLocation = limitResultsToSpecificLocation,
            inputRegions = inputRegions,
            mergeOverlapped = mergeOverlapped,
            genomeMaskedAreaPath = genomeMaskedAreaPath,
            genomeAllowedAreaPath = genomeAllowedAreaPath,
            parallelism = parallelism,
            samplingWithReplacement = samplingWithReplacement
        );
    }

    fun processSharedEnrichmentOptions(
        options: OptionSet,
        metrics: List<String>,
        logger: Logger
    ): SharedEnrichmentOptions {
        val backGroundPath = options.valueOf("background") as Path?
        logger.info("BACKGROUND: $backGroundPath")

        val loiFolderPath = options.valueOf("loi") as Path
        logger.info("LOI TO TEST: $loiFolderPath")

        val loiNameSuffix = options.valueOf("loi-filter") as String?
        logger.info("LOI FNAME SUFFIX: ${loiNameSuffix ?: "N/A"}")

        val hypAlt = options.valueOf("h1") as PermutationAltHypothesis
        logger.info("Alt Hypothesis: $hypAlt")

        val detailedReport = options.has("detailed")
        logger.info("DETAILED_REPORT_FLAG: $detailedReport")

        val aSetIsRegions = !options.has("a-loi")
        if (aSetIsRegions) {
            logger.info("METRIC(a,b): a=INPUT/SIMULATED REGIONs, b=LOI")
        } else {
            logger.info("METRIC(a,b): a=LOI, b=INPUT/SIMULATED REGIONs")
        }

        val aSetFlankedBothSides = options.valueOf("a-flanked") as Int
        logger.info("A FLANKED: $aSetFlankedBothSides")

        val metricStr = options.valueOf("metric") as String
        val metric = when (metricStr) {
            "overlap" -> OverlapNumberMetric(aSetFlankedBothSides)
            "intersection" -> IntersectionNumberMetric(aSetFlankedBothSides)
            else -> throw IllegalArgumentException(
                "Unsupported metric '${metricStr}', use one of: ${metrics.joinToString()}"
            )
        }
        logger.info("METRIC: ${metric.column}")

        return SharedEnrichmentOptions(
            aSetIsRegions = aSetIsRegions,
            metric = metric,
            detailedReport = detailedReport,
            hypAlt = hypAlt,
            backgroundPath = backGroundPath,
            loiFolderPath = loiFolderPath,
            loiNameSuffix = loiNameSuffix,
        )
    }
}

fun parseChrNamesMapping(chrMapStr: List<String>) = chrMapStr.map {
    val pair = it.split(":", limit = 2)
    require(pair.size == 2) {
        "Chrom mapping should be in: '\$chrInput:\$chrGenome' format, but was: $pair"
    }
    pair[0].trim() to pair[1].trim()
}.toMap()

class LoiInfo(
    val label: String,
    val lociFiltered: LocationsList<out RangesList>,
    val processedLoiNumber: Int,
    val recordsNumber: Int
)