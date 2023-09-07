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
        sharedOpts: SharedOptions,
        mergeRegionsToBg: Boolean, // add input regions to background if they are missing there
    ) {
        val loiFolderPath = sharedOpts.loiFolderPath
        val outputBasename = sharedOpts.outputBaseName
        val limitResultsToSpecificLocation = sharedOpts.limitResultsToSpecificLocation

        val reportPath = "${outputBasename}${sharedOpts.metric.column}.tsv".toPath()
        val detailedReportFolder =
            if (sharedOpts.detailedReport) "${outputBasename}${sharedOpts.metric.column}_stats".toPath() else null

        val gq = sharedOpts.genome.toQuery()

        val limitResultsToSpecificLocationFilter: LocationsMergingList? = limitResultsToSpecificLocation?.let {
            LocationsMergingList.create(gq, listOf(limitResultsToSpecificLocation))
        }
        val loiInfos: List<LoiInfo> = loiLocationBaseFilterList(sharedOpts, gq, limitResultsToSpecificLocationFilter, loiFolderPath)

        RegionShuffleStatsFromBedCoverage(
            sharedOpts.simulationsNumber,
            sharedOpts.getChunkSize(), sharedOpts.retries
        ).calcStatistics(
            gq,
            sharedOpts.inputRegions,
            sharedOpts.backgroundPath,
            loiInfos,
            detailedReportFolder,
            metric = sharedOpts.metric,
            hypAlt = sharedOpts.hypAlt,
            aSetIsRegions = sharedOpts.aSetIsRegions,
            mergeOverlapped = sharedOpts.mergeOverlapped,
            intersectionFilter = limitResultsToSpecificLocationFilter,
            genomeMaskedAreaPath = sharedOpts.genomeMaskedAreaPath,
            genomeAllowedAreaPath = sharedOpts.genomeAllowedAreaPath,
            mergeRegionsToBg = mergeRegionsToBg,
            samplingWithReplacement = false
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }

    fun loiLocationBaseFilterList(
        sharedOpts: SharedOptions,
        gq: GenomeQuery,
        loiLocationBaseFilterList: LocationsMergingList?,
        loiFolderPath: Path
    ): List<LoiInfo> {
        val allowedGenomeFilter: LocationsMergingList? = RegionShuffleStats.makeAllowedRegionsFilter(
            sharedOpts.genomeMaskedAreaPath, sharedOpts.genomeAllowedAreaPath, gq
        )
        val loiFilter = when {
            loiLocationBaseFilterList == null -> allowedGenomeFilter
            allowedGenomeFilter == null -> loiLocationBaseFilterList
            else -> allowedGenomeFilter.intersectRanges(loiLocationBaseFilterList) as LocationsMergingList
        }

        val filesStream = if (loiFolderPath.isDirectory) Files.list(loiFolderPath) else Stream.of(loiFolderPath)
        val loiInfos: List<LoiInfo> = collectLoiFrom(
            filesStream, gq, sharedOpts.mergeOverlapped, loiFilter, sharedOpts.loiNameSuffix
        )

        require(loiInfos.isNotEmpty()) {
            "No LOI files passed file suffix filter."
        }
        return loiInfos
    }

    fun collectLoiFrom(
        filesStream: Stream<Path>,
        gq: GenomeQuery,
        mergeOverlapped: Boolean,
        loiFilter: LocationsMergingList?,
        loiNameSuffix: String?
    ): List<LoiInfo> = filesStream.filter { path ->
        loiNameSuffix == null || path.name.endsWith(loiNameSuffix)
    }.map { path ->
        val name = path.fileName.toString()

        val (locList, nLoadedLoi, nRecords) = RegionShuffleStats.readLocationsIgnoringStrand(
            path, gq, mergeOverlapped
        )
        val filteredLocList = if (loiFilter == null) {
            locList
        } else {
            val filteredIterator = locList.intersectRanges(loiFilter).locationIterator()
            when {
                mergeOverlapped -> LocationsMergingList.create(gq, filteredIterator)
                else -> LocationsSortedList.create(gq, filteredIterator)
            }
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
                "Metric is applied to (a,b) where by default 'a' is input/simulated regions regions set," +
                        " 'b' is LOI to check set. This option swaps a and b."
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

            val tool = BioinfToolsCLA.Tools.ENRICHMENT_IN_LOI
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val mergeRegionsToBg = options.has("add-regions-to-bg")
                LOG.info("MERGE REGIONS INTO BACKGROUND: $mergeRegionsToBg")

                val sharedOpts = processSharedOptions(options, metrics, LOG)

                doCalculations(sharedOpts, mergeRegionsToBg)
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

    class SharedOptions(
        val aSetIsRegions: Boolean,
        val mergeOverlapped: Boolean,
        val metric: RegionsMetric,
        val retries: Int,
        val detailedReport: Boolean,
        private val chunkSize: Int,
        val hypAlt: PermutationAltHypothesis,
        val parallelism: Int?,
        val simulationsNumber: Int,
        val outputBaseName: Path,
        val limitResultsToSpecificLocation: Location?,
        val backgroundPath: Path?,
        val genome: Genome,
        val inputRegions: Path,
        val loiFolderPath: Path,
        val loiNameSuffix: String?,
        val genomeMaskedAreaPath: Path?,
        val genomeAllowedAreaPath: Path?,
    ) {
        fun getChunkSize(): Int = if (chunkSize == 0) simulationsNumber else chunkSize
    }

    fun processSharedOptions(
        options: OptionSet,
        metrics: List<String>,
        logger: Logger
    ): SharedOptions {
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

        val inputRegions = options.valueOf("regions") as Path
        LOG.info("INPUT REGIONS: $inputRegions")

        val genomeMaskedAreaPath = options.valueOf("genome-masked") as Path?
        LOG.info("GENOME MASKED AREA: $genomeMaskedAreaPath")
        val genomeAllowedAreaPath = options.valueOf("genome-allowed") as Path?
        LOG.info("GENOME ALLOWED AREA: $genomeAllowedAreaPath")

        val backGroundPath = options.valueOf("background") as Path?
        LOG.info("BACKGROUND: $backGroundPath")

        val baseName = options.valueOf("basename") as String
        val outputBaseName = FileSystems.getDefault().getPath(baseName).normalize().toAbsolutePath()
        logger.info("OUTPUT_BASENAME: $outputBaseName")

        val loiFolderPath = options.valueOf("loi") as Path
        LOG.info("LOI TO TEST: $loiFolderPath")

        val loiNameSuffix = options.valueOf("loi-filter") as String?
        LOG.info("LOI FNAME SUFFIX: ${loiNameSuffix ?: "N/A"}")

        val limitResultsToSpecificLocation = (options.valueOf("limit-intersect") as String?)?.let {
            parseLocation(
                it, genome, chromSizesPath
            )
        }
        LOG.info("INTERSECT LOI WITH RANGE: ${limitResultsToSpecificLocation ?: "N/A"}")

        val parallelism = options.valueOf("parallelism") as Int?
        configureParallelism(parallelism)
        logger.info("THREADS: $parallelism (use #threads=${parallelismLevel()})")

        val hypAlt = options.valueOf("h1") as PermutationAltHypothesis
        LOG.info("Alt Hypothesis: $hypAlt")

        val simulationsNumber = options.valueOf("simulations") as Int
        LOG.info("SIMULATIONS: $simulationsNumber")

        val chunkSize = options.valueOf("chunk-size") as Int
        LOG.info("SIMULATIONS CHUNK SIZE: $chunkSize")

        val detailedReport = options.has("detailed")
        LOG.info("DETAILED_REPORT_FLAG: $detailedReport")

        val retries = options.valueOf("retries") as Int
        logger.info("MAX RETRIES: $retries")

        val aSetIsRegions = options.has("a-loi")
        if (aSetIsRegions) {
            logger.info("METRIC(a,b): a=INPUT/SIMULATED REGIONs, b=LOI")
        } else {
            logger.info("METRIC(a,b): a=LOI, b=INPUT/SIMULATED REGIONs")
        }

        val aSetFlankedBothSides = options.valueOf("a-flanked") as Int
        logger.info("A FLANKED: $aSetFlankedBothSides")

        val mergeOverlapped = options.has("merge")
        logger.info("MERGE OVERLAPPED: $mergeOverlapped")

        val metricStr = options.valueOf("metric") as String
        val metric = when (metricStr) {
            "overlap" -> OverlapNumberMetric(aSetFlankedBothSides)
            "intersection" -> IntersectionNumberMetric(aSetFlankedBothSides)
            else -> throw IllegalArgumentException(
                "Unsupported metric '${metricStr}', use one of: ${metrics.joinToString()}"
            )
        }
        logger.info("METRIC: ${metric.column}")

        return SharedOptions(
            aSetIsRegions = aSetIsRegions,
            mergeOverlapped = mergeOverlapped,
            metric = metric,
            retries = retries,
            detailedReport = detailedReport,
            chunkSize = chunkSize,
            hypAlt = hypAlt,
            parallelism = parallelism,
            simulationsNumber = simulationsNumber,
            outputBaseName = outputBaseName,
            limitResultsToSpecificLocation = limitResultsToSpecificLocation,
            backgroundPath = backGroundPath,
            genome = genome,
            inputRegions = inputRegions,
            loiFolderPath = loiFolderPath,
            loiNameSuffix = loiNameSuffix,
            genomeMaskedAreaPath = genomeMaskedAreaPath,
            genomeAllowedAreaPath = genomeAllowedAreaPath
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
    val loci: LocationsList<out RangesList>,
    val processedLoiNumber: Int,
    val recordsNumber: Int
)