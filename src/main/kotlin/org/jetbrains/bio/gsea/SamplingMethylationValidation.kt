package org.jetbrains.bio.gsea

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.genome.sampling.shuffleChromosomeRanges
import org.jetbrains.bio.gsea.EnrichmentInLoi.processInputRegions
import org.jetbrains.bio.gsea.RegionShuffleStats.Companion.sampleRegions
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.ThreadLocalRandom
import kotlin.math.ceil
import kotlin.math.max
import kotlin.math.min
import kotlin.random.Random
import kotlin.random.asKotlinRandom

object SamplingMethylationValidation {
    private val LOG = LoggerFactory.getLogger(SamplingMethylationValidation::class.java)
    private val SUPPORTED_ALGS = listOf("METH", "BED")
    private const val DIST_CORRECTION_MAX_REGION_SIZE = 5000 // # TODO: customize?
    private const val DIST_CORRECTION_APPROACH_NAME = "dist"

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
                listOf("methylome"),
                "Covered methylome positions to use as random background. Strand is ignored. Must " +
                        "intersect with each of given input region. File should be TAB separated, 1st column is chromosome name, " +
                        "second column is 1-based offset. E.g. CpG offsets from DMC or proportions table."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            // Force methylome to have 0-based offsets
            accepts(
                "methylome-zero",
                "Methylome offsets in 0-based format instead of 1-based."
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
                listOf("limit-intersect"),
                "Truncate LOI & simulated results by range 'chromosome:start-end'" +
                        " before over-representation test. 'start' offset 0-based, 'end' offset exclusive, i.e" +
                        " [start, end). Regions will be simulated from background and truncated before metric " +
                        " calculation. E.g. could be used to estimate over-representation in context of specific " +
                        " region, e.g. HLA locus. It is not the same as limit background only to HLA locus and " +
                        " sample input regions that also overlap the locus."
            )
                .withRequiredArg()

            // Output
            acceptsAll(
                listOf("o", "basename"),
                "Output files paths will start with given path prefix. Prefix could be relative" +
                        " or absolute path or it's part including folder name ending with '/'."
            )
                .withRequiredArg()
                .defaultsTo("regions_in_loi_regions_enrichment_")

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
                "replace",
                "Sample with replacement, i.e. sampled regions could intersect each other. By default w/o replacement"
            )

            acceptsAll(listOf("parallelism"), "parallelism level")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)

            acceptsAll(
                listOf("m", "merge"),
                "Merge overlapped input regions while loading files if applicable. Will speedup computations."
            )

            // simulation
            acceptsAll(
                listOf("alg"),
                "Simulation algorithm. Supported options: ${SUPPORTED_ALGS.joinToString()}." +
                        " `METH` samples regions from methylome positions table provided." +
                        " The sampling preserves the distribution of cytosines in sampled regions as in input regions. " +
                        "`BED` samples regions from methylome converted to BED background. In this case background will" +
                        " be a set of merged regions, where each cytosine offset will be converted to a region" +
                        " [offset - flnk, offset + flnk]. `flnk` is flanking region radius. Such sampling preserves" +
                        " the distribution of regions genomic length taken from from input regions."
            )
                .withRequiredArg()

            acceptsAll(listOf("flnk"), "BED background ")
                .withRequiredArg()
                .withValuesConvertedBy(IntConverter)
                .defaultsTo(50)
            accepts(
                "add-regions-to-bg",
                "Merge input regions into background before sampling if background doesn't cover all input regions."
            )
            accepts(
                "ignore-regions-out-of-bg",
                "Do not use in analysis the input region that don't intersect background."
            )


            // Sampled distribution correction
            acceptsAll(
                listOf("length-correction"),
                "Optional setting that allow to do a correction of sampled length distribution. Use " +
                        "'dist' to make sampling result closer to input length distribution. Use int value as a " +
                        "threshold of maximum allowed length. Is supported only by `METH` simulation alg."
            )
                .withRequiredArg()

            // Logging level:
            acceptsAll(listOf("d", "debug"), "Print all the debug info")

            val tool = BioinfToolsCLA.Tools.METHYLATION_SAMPLING_VALIDATION_IN_LOI
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val  methylomePath = options.valueOf("methylome") as Path
                LOG.info("METHYLOME: $methylomePath")

                val zeroBasedMethylome = "methylome-zero" in options
                LOG.info("0-BASED OFFSETS IN METHYLOME: $zeroBasedMethylome")

                val algStr = options.valueOf("alg") as String
                LOG.info("ALG: $algStr")
                if (algStr !in SUPPORTED_ALGS) {
                    LOG.error("Unsupported algorithm name '$algStr', should be one of: ${SUPPORTED_ALGS.joinToString()}")
                }
                val sampleFromBEDBackground = (algStr == "BED")

                val bedBgFlnk: Int
                val mergeRegionsToBedBg: Boolean
                val ignoreRegionsOutOfBg: Boolean
                if (sampleFromBEDBackground) {
                    bedBgFlnk = options.valueOf("flnk") as Int
                    LOG.info("BED BACKGROUND FLANKING RADIUS: ${bedBgFlnk.formatLongNumber()}")

                    mergeRegionsToBedBg = options.has("add-regions-to-bg")
                    LOG.info("MERGE REGIONS INTO BED BACKGROUND: $mergeRegionsToBedBg")

                    ignoreRegionsOutOfBg = options.has("ignore-regions-out-of-bg")
                    LOG.info("IGNORE REGIONS OUT OF BACKGROUND: $ignoreRegionsOutOfBg")
                } else {
                    bedBgFlnk = 0
                    // if not only DMRs are used as input, it could be useful to support it:
                    mergeRegionsToBedBg = false
                    ignoreRegionsOutOfBg = false
                }

                val lengthCorrectionMethod = options.valueOf("length-correction") as String?
                LOG.info("SAMPLED REGIONS LENGTH CORRECTION: $lengthCorrectionMethod")
                require(lengthCorrectionMethod == null || !sampleFromBEDBackground) {
                    "Length correction is supported only for `METH` algorithm but was: $algStr"
                }

                val sharedOpts = EnrichmentInLoi.processSharedSamplingOptions(options, LOG)

                doCalculations(
                    sharedOpts,
                    methylomePath = methylomePath,
                    zeroBasedMethylome = zeroBasedMethylome,
                    sampleFromBEDBackground = sampleFromBEDBackground,
                    bedBgFlnk = bedBgFlnk,
                    mergeRegionsToBedBg = mergeRegionsToBedBg,
                    ignoreRegionsOutOfBg = ignoreRegionsOutOfBg,
                    lengthCorrectionMethod = lengthCorrectionMethod
                )
            }
        }
    }

    fun doCalculations(
        opts: EnrichmentInLoi.SharedSamplingOptions,
        methylomePath: Path,
        zeroBasedMethylome: Boolean,
        sampleFromBEDBackground: Boolean,
        bedBgFlnk: Int,
        mergeRegionsToBedBg: Boolean,
        ignoreRegionsOutOfBg: Boolean,
        lengthCorrectionMethod: String?
    ) {

        val truncateRangesToSpecificLocation = opts.truncateRangesToSpecificLocation
        val gq = opts.genome.toQuery()

        val detailedReportFolder = "${opts.outputBaseName}sampling${opts.simulationsNumber}_stats".toPath()

        //-------------------
        val truncateRangesToSpecificLocationFilter: LocationsMergingList? = truncateRangesToSpecificLocation?.let {
            LocationsMergingList.create(gq, listOf(truncateRangesToSpecificLocation))
        }
        val genomeAllowedAreaFilter = opts.genomeAllowedAreaPath?.let {
            RegionShuffleStats.readGenomeAreaFilter(it, gq)
        }
        val genomeMaskedAreaFilter = opts.genomeMaskedAreaPath?.let {
            RegionShuffleStats.readGenomeAreaFilter(it, gq)
        }

        //-------------------

        // Load Input regions, methylome + filter allowed regions / methylome
        var inputRegionsFiltered = processInputRegions(
            opts.inputRegions, gq,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter
        )

        val backgroundMethSampling = RegionShuffleStatsFromMethylomeCoverage.backgroundProviderFun(
            methylomePath, zeroBasedMethylome, gq,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter,
        )
        RegionShuffleStatsFromMethylomeCoverage.ensureInputRegionsMatchesBackgound(
            inputRegionsFiltered,
            backgroundMethSampling,
            methylomePath
        )

        // Make a filtered BED backgroundMethSampling from methylome if required
        val filteredBedBackgroundFromMethylome = if (sampleFromBEDBackground) {
            makeFilteredBEDBackgroundFromMethylome(
                gq,
                backgroundMethSampling.anchorMethylome,
                bedBgFlnk,
                inputRegionsFiltered,
                mergeRegionsToBedBg
            )
        } else {
            null
        }

        if (mergeRegionsToBedBg) {
            require(!ignoreRegionsOutOfBg) { "Cannot merge & ignore region at the same time." }
        }

        if (filteredBedBackgroundFromMethylome != null) {
            if (ignoreRegionsOutOfBg) {
                val nRegionsBeforeBgFilter = inputRegionsFiltered.size
                inputRegionsFiltered = inputRegionsFiltered.filter {
                    filteredBedBackgroundFromMethylome.intersects(it)
                }
                LOG.info("Filtering input regions that doesn't intersect background. Removed " +
                        "${nRegionsBeforeBgFilter-inputRegionsFiltered.size} of $nRegionsBeforeBgFilter regions")
            }

            // validate rest regions:
            inputRegionsFiltered.forEach {
                require(filteredBedBackgroundFromMethylome.intersects(it)) {
                    "BED Background regions are supposed to intersect all loci of interest, but the " +
                            "region is missing in bg: ${it.toChromosomeRange()}"
                }
            }
        }
        require(inputRegionsFiltered.isNotEmpty()) {
            "Regions file is empty or all regions were removed by filters."
        }

        val inputRegionsListFiltered = when {
            opts.mergeOverlapped -> LocationsMergingList.create(gq, inputRegionsFiltered)
            else -> LocationsSortedList.create(gq, inputRegionsFiltered)
        }
        val inputRegionsPreprocessed = inputRegionsListFiltered.asLocationSequence().toList()

        // Use full methylome, not anchor. Because finally sampled regions will be sampled from full methylome
        // with a point from anchor methylome
        val methylomeForRegionsCoverageHistograms = backgroundMethSampling.fullMethylome

        // Sample regions
        if (!sampleFromBEDBackground) {
            // Sample from methylome BG
            sampleAndCollectMetrics(
                inputRegionsPreprocessed,
                simulationsNumber = opts.simulationsNumber,
                chunkSize = opts.getChunkSize(),
                regionSetMaxRetries = opts.setRetries,
                singleRegionMaxRetries = opts.regionRetries,
                gq = gq,
                outputFolderPath = detailedReportFolder,
                truncateFilter = truncateRangesToSpecificLocationFilter,
                samplingWithReplacement = opts.samplingWithReplacement,
                methylomeForRegionsCoverageHistograms = methylomeForRegionsCoverageHistograms,
                samplingBackground = backgroundMethSampling,
                samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
                    val threadLocalRandom = ThreadLocalRandom.current().asKotlinRandom()
                    val res = RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
                        genomeQuery,
                        regions,
                        background,
                        singleRegionMaxRetries = singleRegionMaxRetries,
                        regionSetMaxRetries = regionSetMaxRetries,
                        withReplacement = withReplacement,
                        candidateFilterPredicate = getCandidateFilterPredicate(
                            lengthCorrectionMethod,
                            regions,
                            randomGen = threadLocalRandom
                        ),
                        genomeMaskedAreaFilter = genomeMaskedAreaFilter
                    )
                    res.first.map { it.on(Strand.PLUS) } to res.second
                },
            )
        } else {
            // Sample from BED bg
            sampleAndCollectMetrics(
                inputRegionsPreprocessed,
                simulationsNumber = opts.simulationsNumber,
                chunkSize = opts.getChunkSize(),
                regionSetMaxRetries = opts.setRetries,
                singleRegionMaxRetries = opts.regionRetries,
                gq = gq,
                outputFolderPath = detailedReportFolder,
                truncateFilter = truncateRangesToSpecificLocationFilter,
                samplingWithReplacement = opts.samplingWithReplacement,
                methylomeForRegionsCoverageHistograms = methylomeForRegionsCoverageHistograms,
                samplingBackground = filteredBedBackgroundFromMethylome!!,
                samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
                    val threadLocalRandom = ThreadLocalRandom.current().asKotlinRandom()
                    val res = shuffleChromosomeRanges(
                        genomeQuery,
                        regions,
                        background.asLocationSequence().map { it.toChromosomeRange() }.toList(),
                        regionSetMaxRetries = regionSetMaxRetries,
                        singleRegionMaxRetries = singleRegionMaxRetries,
                        withReplacement = withReplacement,
                        genomeMaskedAreaFilter = genomeMaskedAreaFilter,
                        randomGen = threadLocalRandom
                    )
                    res.first.map { it.on(Strand.PLUS) } to res.second
                }
            )
        }
    }

    fun getCandidateFilterPredicate(
        lengthCorrectionMethod: String?,
        inputRegionsPreprocessed: List<ChromosomeRange>,
        randomGen: Random
    ): ((ChromosomeRange) -> Double)? {
        return when (lengthCorrectionMethod) {
            null -> null // XXX: default
            DIST_CORRECTION_APPROACH_NAME -> makeFilterByLengthProbability(inputRegionsPreprocessed, randomGen)
            else -> {
                // By fixed threshold:
                val lenThreshold = lengthCorrectionMethod.toInt().toDouble() //e.g. UP: 2218, DOWN: 1261/1472

                val f: (ChromosomeRange) -> Double = { chrRange ->
                    val candLen = chrRange.length()
                    if (candLen <= lenThreshold) 0.0 else min(1.0, candLen / lenThreshold)
                }
                f
            }
        }
    }

    private fun makeFilterByLengthProbability(
        inputRegionsPreprocessed: List<ChromosomeRange>,
        randomGen: Random
    ): ((ChromosomeRange) -> Double) {
        // Get features of input regions distribution for fast filter calculation
        val resampleProb = determineResampleProbabilities(inputRegionsPreprocessed)
        val maxProcessedLen = resampleProb.size - 1 // last index of array

        // Filter
        return { chromosomeRange ->
            // XXX based on sampling dist
            val candLen = chromosomeRange.length()
            val prob = if (candLen > maxProcessedLen) {
                1.0
            } else {
                resampleProb[candLen]
            }
            val accepted = randomGen.nextDouble() >= prob
            if (accepted) 0.0 else prob
        }
    }

    private fun determineResampleProbabilities(inputRegionsPreprocessed: List<ChromosomeRange>): DoubleArray {
        val inputLen = inputRegionsPreprocessed.map { min(DIST_CORRECTION_MAX_REGION_SIZE, it.length()) }.toIntArray()
        val regionsWithMaxSupportedLen = inputLen.count { it == DIST_CORRECTION_MAX_REGION_SIZE }
        require((regionsWithMaxSupportedLen / inputLen.size) < 0.01) { // if >=1% of unsupported lengths
            // TODO: add an option to ignore it?
            "Too many regions (more than 5%) with length >= $DIST_CORRECTION_MAX_REGION_SIZE bp." +
                    " #Regions exceed threshold = $regionsWithMaxSupportedLen of ${inputLen.size}." +
                    " Consider changing threshold."
        }
        inputLen.sort() // sort from min to max
        val minLen = inputLen.first()
        val maxLen = min(DIST_CORRECTION_MAX_REGION_SIZE, inputLen.last())
        val resampleProb = DoubleArray(maxLen + 2)
        resampleProb[maxLen + 1] = 1.0
        for (idx in 0 until minLen) {
            resampleProb[idx] = 1.0
        }
        val medianIdx = inputLen.size / 2.0
        for (idx in minLen..maxLen) {
            val cntAbove = inputLen.count { it <= idx }
            val cntBelow = inputLen.count { idx <= it }
            val pval = min(cntAbove, cntBelow) / medianIdx
            resampleProb[idx] = 1.0 - pval
        }
        resampleProb[maxLen + 1] = 1.0
        return resampleProb
    }

    private fun makeFilteredBEDBackgroundFromMethylome(
        gq: GenomeQuery,
        methylomeCov: BasePairCoverage,
        bedBgFlnk: Int,
        inputRegions: List<Location>,
        mergeRegionsToBg: Boolean,
    ): LocationsMergingList {
        // Here no need in intersecting input or background with `SharedSamplingOptions.truncateRangesToSpecificLocation`

        val bgLoci = makeBEDBackgroundFromMethylome(
            gq, methylomeCov, bedBgFlnk,
            if (mergeRegionsToBg) inputRegions else emptyList()
        )

        return bgLoci
    }

    fun makeBEDBackgroundFromMethylome(
        gq: GenomeQuery,
        methylomeCov: BasePairCoverage,
        bedBgFlnk: Int,
        additionalLoci: List<Location>
    ): LocationsMergingList {
        LOG.info("Make BED Background from methylome")
        require(gq == methylomeCov.data.genomeQuery)

        val builder = LocationsMergingList.builder(gq)
        additionalLoci.forEach {
            builder.add(it)
        }
        gq.get().forEach { chr ->
            // We have too many cytosines in methylome, let's merge on the fly:
            val offsets = methylomeCov.data[chr]
            var prevStartOffset = 0
            var prevEndOffset = 0
            for (nextOffset in offsets) {
                // close prev location & start new
                val nextStart = max(0, nextOffset - bedBgFlnk)
                val nextEnd = min(nextOffset + bedBgFlnk, chr.length)
                if (prevEndOffset < nextStart) {
                    if (prevStartOffset != prevEndOffset) {
                        builder.add(Location(prevStartOffset, prevEndOffset, chr))
                    }
                    prevStartOffset = nextStart
                    prevEndOffset = nextEnd
                } else {
                    // elongate
                    prevEndOffset = nextEnd
                }
            }
            if (prevStartOffset != prevEndOffset) {
                builder.add(Location(prevStartOffset, prevEndOffset, chr))
            }
        }
        return builder.build()
    }

    private fun <T> sampleAndCollectMetrics(
        inputRegions: List<Location>,
        simulationsNumber: Int,
        chunkSize: Int,
        regionSetMaxRetries: Int,
        singleRegionMaxRetries: Int,
        gq: GenomeQuery,
        outputFolderPath: Path,
        truncateFilter: LocationsList<out RangesList>? = null,
        samplingWithReplacement: Boolean = false,
        methylomeForRegionsCoverageHistograms: BasePairCoverage,
        samplingBackground: T,
        samplingFun: (GenomeQuery, List<ChromosomeRange>, T, Int, Int, Boolean) -> Pair<List<Location>, IntHistogram>
    ) {
        outputFolderPath.createDirectories()

        var path: Path

        LOG.info("Calc input regions size distribution...")
        path = outputFolderPath / "input_regions.length.hist.tsv"
        IntHistogram.create(
            inputRegions.map { it.length() }.toIntArray()
        ).save(path, metric_col = "length")
        LOG.info("[DONE]: $path")

        LOG.info("Calc input regions overlap with background distribution...")
        path = outputFolderPath / "input_regions.bg_overlap.hist.tsv"
        IntHistogram.create(
            inputRegions.map { methylomeForRegionsCoverageHistograms.getCoverage(it) }.toIntArray()
        ).save(path, metric_col = "n_cpg")
        LOG.info("[DONE]: $path")

        val nChunks = ceil(simulationsNumber.toDouble() / chunkSize).toInt()

        LOG.info("Do simulation validation...")
        val progress = Progress { title = "Simulation & calc results stats progress (all chunks)" }.bounded(
            nChunks.toLong()
        )
        val sampledCummulativeLengthDist = IntHistogram()
        val sampledCummulativeBgOverlapDist = IntHistogram()
        val finalRegionAttemptsHist = IntHistogram()
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info(
                "Simulations: Chunk [${chunkId + 1} of $nChunks], simulations " +
                        "${start.formatLongNumber()}..${end.formatLongNumber()} of ${simulationsNumber.formatLongNumber()}"
            )

            val (sampledRegions: List<List<LocationsList<out RangesList>>>, regionAttemptsHist) =
                sampleRegions(
                    simulationsInChunk,
                    regionSetMaxRetries = regionSetMaxRetries,
                    singleRegionMaxRetries = singleRegionMaxRetries,
                    truncateFilter = truncateFilter,
                    srcLoci = inputRegions,
                    background = samplingBackground,
                    gq = gq,
                    parallelismLevel(),
                    samplingFun = samplingFun,
                    withReplacement = samplingWithReplacement
                )
            finalRegionAttemptsHist += regionAttemptsHist

            val chunkProgress = Progress { title = "Calculate stats for sampled chunk data" }.bounded(
                sampledRegions.sumOf { chunk -> chunk.sumOf { ll -> ll.size.toLong() } }
            )
            sampledRegions.forEach { chunk ->
                chunk.forEach { sampledSet ->
                    sampledSet.asLocationSequence().forEach { loc ->
                        sampledCummulativeLengthDist.increment(loc.length())
                        sampledCummulativeBgOverlapDist.increment(methylomeForRegionsCoverageHistograms.getCoverage(loc))

                        chunkProgress.report()
                    }
                }
            }
            chunkProgress.done()
            progress.report()
        }
        progress.done()

        // ----------------
        path = outputFolderPath / "sampled_${simulationsNumber}.region_attempts.hist.tsv"
        finalRegionAttemptsHist.save(path, metric_col = "attempt")
        LOG.info("[DONE]: $path")

        path = outputFolderPath / "sampled_${simulationsNumber}.length.hist.tsv"
        sampledCummulativeLengthDist.save(path, metric_col = "length")
        LOG.info("[DONE]: $path")

        path = outputFolderPath / "sampled_${simulationsNumber}.bg_overlap.hist.tsv"
        sampledCummulativeBgOverlapDist.save(path, metric_col = "n_cpg")
        LOG.info("[DONE]: $path")
    }
}