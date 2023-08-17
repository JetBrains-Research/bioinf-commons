package org.jetbrains.bio.gse

import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import joptsimple.OptionParser
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.dataframe.*
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.*
import org.jetbrains.bio.genome.containers.intersection.IntersectionNumberMetric
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.coverage.binarySearchLeft
import org.jetbrains.bio.genome.methylome.Methylome
import org.jetbrains.bio.genome.methylome.MethylomeBuilder
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream
import kotlin.UnsupportedOperationException

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
            /*

                        accepts(
                            "add-loi-to-bg",
                            "Add loci of interest (loi) to background before sampling."
                        )

                        acceptsAll(
                            listOf("genome-masked"),
                            "Mask genome regions: loi intersecting masked regions are skipped, masked regions removed from background." +
                                    " File path in TAB separated BED with at least 3 fields. Strand is ignored."
                        )
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))

                        acceptsAll(
                            listOf("genome-allowed"),
                            "Allow only certain genome regions: loi/background not-intersecting masked regions are skipped." +
                                    " File path in TAB separated BED with at least 3 fields. Strand is ignored."
                        )
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))

            */
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


                                    // TODO: optional save sampledRegions
                                    // TODO: optional load sampledRegions from ...

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

                                    acceptsAll(
                                        listOf("chunk-size"), "To reduce mem usage, simulation is done in several steps with" +
                                                " provided simulations number per chunk. Use 0 to sample all sets at once."
                                    )
                                        .withRequiredArg()
                                        .ofType(Int::class.java)
                                        .defaultsTo(50000)

                                    acceptsAll(listOf("retries"), "Region shuffling max retries")
                                        .withRequiredArg()
                                        .ofType(Int::class.java)
                                        .defaultsTo(1000)

            /*
                                    accepts(
                                        "a-regions",
                                        "Metric is applied to (a,b) where by default 'a' is loi/simulated loi set, 'b' is region to check set. This" +
                                                " option switches a and b definitions."
                                    )
                                    */
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

                //TODO val genomeMaskedLociPath = options.valueOf("genome-masked") as Path?
                val genomeMaskedLociPath: Path? = null
                LOG.info("GENOME MASKED LOCI: $genomeMaskedLociPath")

                //TODO val genomeAllowedLociPath = options.valueOf("genome-allowed") as Path?
                val genomeAllowedLociPath: Path? = null
                LOG.info("GENOME ALLOWED LOCI: $genomeAllowedLociPath")

                val methylomeBackground = options.valueOf("background") as Path?
                LOG.info("METHYLOME BACKGROUND: $methylomeBackground")

                //TODO val addLoiToBg = options.has("add-loi-to-bg")
                val addLoiToBg = false
                LOG.info("ADD LOI TO BG: $addLoiToBg")

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
                LOG.info("MAX_RETRIES: $retries")

                val parallelism = options.valueOf("parallelism") as Int?
                LOG.info("THREADS: $parallelism")
                configureParallelism(parallelism)

                val detailedReport = options.has("detailed")
                LOG.info("DETAILED_REPORT_FLAG: $detailedReport")

                //TODO: val aSetIsLoi = !options.has("a-regions")
                val aSetIsLoi = true
                if (aSetIsLoi) {
                    LOG.info("METRIC(a,b): a=REGION, b=LOI/SIMULATED LOI")
                } else {
                    LOG.info("METRIC(a,b): a=LOI/SIMULATED LOI, b=REGION")
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
                    srcLoci, methylomeBackground, addLoiToBg, regionsFolderPath, genome,
                    simulationsNumber,
                    if (chunkSize == 0) simulationsNumber else chunkSize,
                    outputBaseName,
                    metric,
                    detailedReport, retries, hypAlt,
                    aSetIsLoi,
                    mergeOverlapped,
                    regionsToTestFilterLocation,
                    regionsNameSuffix,
                    genomeMaskedLociPath, genomeAllowedLociPath
                )
            }
        }
    }

    fun doCalculations(
        loiPath: Path,
        methylomeBackground: Path?,
        addLoiToBg: Boolean,
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
        genomeAllowedLociPath: Path?
    ) {
        val reportPath = "${outputBasename}${metric.column}.tsv".toPath()
        val detailedReportFolder = if (detailed) "${outputBasename}${metric.column}_stats".toPath() else null

        val regionsToTestFilter = if (regionsToTestFilterLocation != null) {
            LocationsSortedList.create(
                genome.toQuery(),
                listOf(regionsToTestFilterLocation)
            )
        } else {
            null
        }

        val regionsAndRanges: List<Pair<String, LocationsList<out RangesList>>> =
            EnrichmentInRegions.collectRegionsFrom(
                if (regionsPath.isDirectory) Files.list(regionsPath) else Stream.of(regionsPath),
                genome, mergeOverlapped, regionsToTestFilter, regionsNameSuffix
            )

        require(regionsAndRanges.isNotEmpty()) {
            "No regions files passed file suffix filter."
        }

        RegionShuffleStats1(
            genome,
            simulationsNumber, chunkSize,
            maxRetries
        ).calcStatistics(
            loiPath, methylomeBackground,
            regionsAndRanges,
            detailedReportFolder,
            metric = metric,
            hypAlt = hypAlt,
            aSetIsLoi = aSetIsLoi,
            mergeOverlapped = mergeOverlapped,
            intersectionFilter = regionsToTestFilter,
            genomeMaskedLociPath = genomeMaskedLociPath,
            genomeAllowedLociPath = genomeAllowedLociPath,
            addLoiToBg = addLoiToBg
        ).save(reportPath)

        LOG.info("Report saved to: $reportPath")
        if (detailedReportFolder != null) {
            LOG.info("Report details saved to: $detailedReportFolder")
        }
    }
}

class RegionShuffleStats1(
    private val genome: Genome,
    private val simulationsNumber: Int,
    private val chunkSize: Int,
    private val maxRetries: Int
) {

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStats1::class.java)

    }

    fun calcStatistics(
        srcRegionsPath: Path,
        methylomeBackground: Path?,
        regionLabelAndLociToTest: List<Pair<String, LocationsList<out RangesList>>>,
        outputFolderPath: Path? = null,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER,
        aSetIsLoi: Boolean = true,
        mergeOverlapped: Boolean = true,
        intersectionFilter: LocationsList<out RangesList>? = null,
        genomeMaskedLociPath: Path? = null,
        genomeAllowedLociPath: Path? = null,
        addLoiToBg: Boolean = false
    ): DataFrame {
        outputFolderPath?.createDirectories()
        val dumpDetails = outputFolderPath != null

        val gq = genome.toQuery()

        val (sourceLoci, bgLociList) = loadSrcLociAndBg(
            srcRegionsPath, methylomeBackground, addLoiToBg, genomeMaskedLociPath, genomeAllowedLociPath, gq
        )
        require(sourceLoci.isNotEmpty()) {
            "Loci file is empty or all loci were masked."
        }
        return DataFrame()
    }

    private fun loadSrcLociAndBg(
        srcRegionsPath: Path, backgroundRegionsPath: Path?,
        addLoiToBg: Boolean,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        gq: GenomeQuery
    ): Pair<List<Location>, LocationsMergingList> {
        val sourceLoci = RegionShuffleStats.readLocationsIgnoringStrand(srcRegionsPath, gq)
        LOG.info("Source loci: ${sourceLoci.size} regions")

        val bgLoci = if (backgroundRegionsPath == null) {
            // Whole Genome
            // TODO: make BG to be obligatory
            throw UnsupportedOperationException("Whole genome BG not supported!")
        } else {
            //TODO: (chr, BitList()) ?
            //TODO: or (chr, Dataframe()) ?

            // TODO: 1: load all . from BG
            // TODO: 2: [how bed shuffle works, what for prefix sum, see shuffleChromosomeRanges(..)
            // TODO: 3:
            //   for each k_cpg in cpg_lengths:
            //          find random position from all: chr x all positions
            //              here could be concept of prefix sum: chr, #items, so select random . => know chromosome & offset there
            //          make region using convert k_cpg to true length if it is in chromosome and not 'maskedGenomeMap'.
            //              if all ok => add to maskedGenomeMap
            val methCovData: GenomeStrandMap<TIntList> = genomeStrandMap(gq) { _, _ ->
                TIntArrayList()
            }
//            gq.genome.chromosomeNamesMap
//            Chromosome.getOrCreateChromosome()

            val header = false
            val format = CSVFormat.TDF.withQuoteMode(QuoteMode.MINIMAL).withCommentMarker('#')
                .withRecordSeparator("\n")!!

            LOG.info("Loading methylation background $backgroundRegionsPath ${backgroundRegionsPath.size}")
            val linesNumber = backgroundRegionsPath.bufferedReader().use {
                it.lines().mapToInt { line -> if (line[0] != format.commentMarker) 1 else 0 }.sum()
            }
            val rowsNumber = if (linesNumber == 0) 0 else linesNumber - (if (header) 1 else 0)
            LOG.info("Rows: ${rowsNumber.formatLongNumber()}")
            val progress = Progress { title = "Reading methylation coverage data $backgroundRegionsPath" }
                .bounded(rowsNumber.toLong())

            backgroundRegionsPath.bufferedReader().use {
                for ((i, row) in format.withSkipHeaderRecord(header).parse(it).withIndex()) {
                    val chrName = row[0]
                    val offset = row[1].toInt() - 1 // methylome coverage table is 1-based, convert to 0-based!

                    // XXX: Not the fastest impl, but work fast enough for WGBS all CpG table
                    val chr = Chromosome(gq.genome, chrName)
                    methCovData[chr, Strand.PLUS].add(offset.toInt())

                    progress.report()
                }

                progress.done()
            }

            if (addLoiToBg) {
                // merge source loci into background
                // TODO bgLocations.addAll(sourceLoci)
                throw UnsupportedOperationException("add loi to bg not supported!")
            }
            val bgLoci = LocationsMergingList.create(
                gq,
                // TODO bgLocations
                emptyList()
            )
            LOG.info("Background regions: ${bgLoci.size} regions")

            sourceLoci.forEach {
                val chrStrandData = methCovData[it.chromosome, it.strand]
                val idx1 = chrStrandData.binarySearchLeft(it.startOffset)
                val idx2 = chrStrandData.binarySearchLeft(it.endOffset)
                print("${it}, idx1=$idx1, idx2=$idx2")
//                require(bgLoci.includes(it)) {
//                    "Background $backgroundRegionsPath regions are required to include all loci of interest, but the " +
//                            "loci is missing in bg: ${it.toChromosomeRange()}"
//                }
            }
            bgLoci
        }

        // complementary to masked list
        val genomeMaskedComplementary = if (genomeMaskedLociPath != null) {
            val maskedGenome = RegionShuffleStats.readLocationsIgnoringStrand(genomeMaskedLociPath, gq)
            LOG.info("Genome masked loci: ${maskedGenome.size} regions")
            val maskedGenomeLocations = LocationsMergingList.create(gq, maskedGenome)
            LOG.info("Genome masked loci (merged): ${maskedGenomeLocations.size} regions")

            // complementary regions
            maskedGenomeLocations.apply { rl, chr, _ -> rl.complementaryRanges(chr.length) }
        } else {
            null
        }

        // allowed list
        val genomeAllowed = if (genomeAllowedLociPath != null) {
            val allowedGenome = genomeAllowedLociPath.let {
                RegionShuffleStats.readLocationsIgnoringStrand(it, gq)
            }
            LOG.info("Genome allowed loci: ${allowedGenome.size} regions")
            val allowedGenomeLocations = LocationsMergingList.create(gq, allowedGenome)
            LOG.info("Genome allowed loci (merged): ${allowedGenomeLocations.size} regions")
            allowedGenomeLocations
        } else {
            null
        }

        val allowedFilter = when {
            genomeAllowed == null -> genomeMaskedComplementary
            genomeMaskedComplementary == null -> genomeAllowed
            else -> genomeAllowed.intersectRanges(genomeMaskedComplementary) as LocationsMergingList
        }

        val (allowedBgList, allowedSourceLoci) = if (allowedFilter == null) {
            bgLoci to sourceLoci
        } else {
            val allowedBg = bgLoci.intersectRanges(allowedFilter) as LocationsMergingList
            val allowedSourceLoci = sourceLoci.filter { allowedFilter.includes(it) }

            LOG.info("Background regions (all restrictions applied): ${allowedBg.size} regions")
            LOG.info("Source loci (all restrictions applied): ${allowedSourceLoci.size} regions")
            allowedBg to allowedSourceLoci
        }

        return allowedSourceLoci to allowedBgList
    }
}