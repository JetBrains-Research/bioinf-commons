package org.jetbrains.bio.gsea

import gnu.trove.map.hash.TIntIntHashMap
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.*
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.genome.sampling.createMaskedArea
import org.jetbrains.bio.util.formatLongNumber
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.ThreadLocalRandom
import kotlin.math.max
import kotlin.math.min
import kotlin.random.asKotlinRandom

// TODO: optional save/load sampled regions to make 100% reproducible

class RegionShuffleStatsFromMethylomeCoverage(
    simulationsNumber: Int,
    chunkSize: Int,
    regionSetMaxRetries: Int,
    singleRegionMaxRetries: Int,
    val zeroBasedBg: Boolean
) : RegionShuffleStats(simulationsNumber, chunkSize, regionSetMaxRetries, singleRegionMaxRetries) {

    fun calcStatistics(
        gq: GenomeQuery,
        inputRegionsPath: Path,
        backgroundPath: Path?,
        loiInfos: List<LoiInfo>,
        outputFolderPath: Path?,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis,
        aSetIsRegions: Boolean,
        mergeOverlapped: Boolean,
        truncateFilter: LocationsList<out RangesList>?,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
        samplingWithReplacement: Boolean,
        lengthCorrectionMethod: String?
    ): DataFrame {
        requireNotNull(backgroundPath) { "Background should be provided" }

        return doCalcStatistics(
            inputRegionsPath,
            gq,
            genomeAllowedAreaFilter,
            genomeMaskedAreaFilter,
            loiInfos,
            outputFolderPath,
            metric,
            hypAlt,
            aSetIsRegions,
            mergeOverlapped,
            truncateFilter,
            samplingWithReplacement = samplingWithReplacement,
            backgroundProvider = { _, inputRegionsFiltered ->
                backgroundProviderFun(
                    inputRegionsFiltered, backgroundPath, zeroBasedBg, gq,
                    genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                    genomeMaskedAreaFilter = genomeMaskedAreaFilter
                )
            },
            loiOverlapWithBgFun = {loiFiltered, background ->
                loiFiltered.asLocationSequence().count {
                    background.fullMethylome.getCoverage(it) > 0
                }
            },
            samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
                val threadLocalRandom = ThreadLocalRandom.current().asKotlinRandom()

                val res = shuffleChromosomeRanges(
                    genomeQuery,
                    regions,
                    background,
                    regionSetMaxRetries = regionSetMaxRetries,
                    singleRegionMaxRetries = singleRegionMaxRetries,
                    withReplacement = withReplacement,
                    genomeMaskedAreaFilter = genomeMaskedAreaFilter,
                    candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                        lengthCorrectionMethod,
                        regions,
                        randomGen = threadLocalRandom
                    )
                )
                res.first.map { it.on(Strand.PLUS) } to res.second
            }
        )
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStatsFromBedCoverage::class.java)

        fun shuffleChromosomeRanges(
            genomeQuery: GenomeQuery,
            regions: List<ChromosomeRange>,
            background: MethylomeSamplingBackground,
            regionSetMaxRetries: Int = 100,
            singleRegionMaxRetries: Int = 100,
            endPositionShift: Int = 2, // For if 'k' is index of latest CpG in candidate, use 'k+2' to get 'exclusive' bound after the latest CpG bounds
            withReplacement: Boolean,
            genomeMaskedAreaFilter: LocationsMergingList?,
            candidateFilterPredicate: ((ChromosomeRange) -> Double)?
        ): Pair<List<ChromosomeRange>, IntHistogram> {

            // Use coverage of input regions for sampling instead of region lengths in basepairs
            // make it one, not in loop
            val expectedCoverages = background.calcExpCoverage(regions)

            // Chromosomes have different amount of data => sampling should take into account that chromosomes probablities
            // should depend on # of observed data on each chromosome
            require(genomeQuery == background.fullMethylome.data.genomeQuery) {
                "Background was made for different genome query. This query=[$genomeQuery]," +
                        " background genome query=[${background.fullMethylome.data.genomeQuery}]"
            }

            val threadLocalRandom = ThreadLocalRandom.current().asKotlinRandom()
            for (i in 1..regionSetMaxRetries) {
                val (sampledRegions, attemptsHistogram) = tryShuffle(
                    genomeQuery,
                    background,
                    expectedCoverages,
                    singleRegionMaxRetries,
                    endPositionShift,
                    withReplacement,
                    genomeMaskedAreaFilter = genomeMaskedAreaFilter,
                    candidateFilterPredicate = candidateFilterPredicate,
                    randomGen = threadLocalRandom
                )

                if (sampledRegions != null) {
                    // TODO: collect stats how many retries is done
                    return sampledRegions to attemptsHistogram
                }
                // TODO: collect stats how many retries is done
                // TODO: collect stats how many retries is done
                // TODO: collect stats how many retries is done
            }
            throw RuntimeException("Too many shuffle attempts, max limit is: $regionSetMaxRetries")
        }



        fun backgroundProviderFun(
            inputRegionsFiltered: List<Location>,
            backgroundPath: Path?,
            zeroBasedBg: Boolean,
            gq: GenomeQuery,
            genomeAllowedAreaFilter: LocationsMergingList?,
            genomeMaskedAreaFilter: LocationsMergingList?,
        ): MethylomeSamplingBackground {

            val methylomeCov = loadMethylomeCovBackground(backgroundPath!!, zeroBasedBg, gq)
            val anchorMethylomeCov = filterMethylomeCovBackground(
                methylomeCov,
                genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                // also exclude, such points cannot be an anchor and will be dismissed after sampling!
                genomeMaskedAreaFilter = genomeMaskedAreaFilter
            )
            ensureInputRegionsMatchesBackgound(inputRegionsFiltered, anchorMethylomeCov, backgroundPath)
            return MethylomeSamplingBackground(methylomeCov, anchorMethylomeCov)
        }

        fun tryShuffle(
            gq: GenomeQuery,
            background: MethylomeSamplingBackground,
            expectedCoverages: IntArray, // number of covered CpG seen in input regions
            singleRegionMaxRetries: Int = 100, // >=1
            endPositionShift: Int, // E.g. if 'k' is index of latest CpG in candidate, use 'k+2' (endPositionShift = 2) to get 'exclusive' bound after latest CpG bounds
            withReplacement: Boolean,
            genomeMaskedAreaFilter: LocationsMergingList?,
            candidateFilterPredicate: ((ChromosomeRange) -> Double)?,
            randomGen: kotlin.random.Random  // Need Random.nextLong(bound) method that is available in Java 17 or in Kotlin
        ): Pair<List<ChromosomeRange>?, IntHistogram> {
            // Trick: used one storage for masked & seen intervals that is forbidden to intersect
            val maskedArea: GenomeMap<MutableList<Range>>? = createMaskedArea(genomeMaskedAreaFilter, withReplacement, gq)

            val result = ArrayList<ChromosomeRange>(expectedCoverages.size)
            val attemptsHistogram = IntHistogram()

            val numOfCandidatesForAnchor: Long = background.achorPrefixSum.last()

            for (i in expectedCoverages.indices) {
                val requiredCoveredPositionsNum = expectedCoverages[i]
                require(requiredCoveredPositionsNum > 0) { "Uncovered regions not supported, i=$i" }

                var nextRange: ChromosomeRange? = null

                var smalestPenalty: Double = Double.POSITIVE_INFINITY
                var bestCandidate: ChromosomeRange? = null

                var lastRangeTry = 0
                for (rangeTry in 1..singleRegionMaxRetries) {
                    lastRangeTry = rangeTry
                    val rVal = randomGen.nextLong(numOfCandidatesForAnchor)
                    val rndStartShiftInCpGs = randomGen.nextInt(requiredCoveredPositionsNum)

                    var candidate: ChromosomeRange? = generateShuffledRegion(
                        rVal,
                        background,
                        requiredCoveredPositionsNum,
                        endPositionShift,
                        rndStartShiftInCpGs = rndStartShiftInCpGs,
                        maskedArea = maskedArea,
                    )
                    if (candidate != null) {
                        if (candidateFilterPredicate != null) {
                            val penalty = candidateFilterPredicate(candidate)
                            if (penalty != 0.0) {
                                if (penalty < smalestPenalty) {
                                    bestCandidate = candidate
                                    smalestPenalty = penalty
                                }
                                // do not accept a current candidate, looking for better choise
                                candidate = null
                            }
                        }
                    }

                    if (candidate != null) {
                        // best candidate found!
                        nextRange = candidate
                        break
                    }
                    // else retry
                }
                attemptsHistogram.increment(lastRangeTry)

                // All available attempts done or result found:
                if (nextRange == null) {
                    // XXX: Cannot find IDEAL candidate => take the most probable
                    nextRange = bestCandidate
                }

                // success
                if (nextRange == null) {
                    // Cannot sample next range in given attempts number
                    return null to attemptsHistogram
                }

                if (!withReplacement) {
                    maskedArea!![nextRange.chromosome].add(nextRange.toRange())
                }
                result.add(nextRange)
            }

            return result to attemptsHistogram
        }

        fun generateShuffledRegion(
            rVal: Long,
            background: MethylomeSamplingBackground,
            requiredCoveredPositionsNum: Int,
            endPositionShift: Int,
            rndStartShiftInCpGs: Int,
            maskedArea: GenomeMap<MutableList<Range>>?,
        ): ChromosomeRange? {
            require(endPositionShift > 0)

            val anchorPrefixSum = background.achorPrefixSum
            val anchorChrs = background.achorPrefixChrs
            // `prefixSum` Used because coverage is grouped by chromosome and chrs covered differently,
            // so the probability of chromosome is weighted.
            // So it allows us to select chromosome with weighted probability


            //NB: assume no duplicated values in prefix sum!!!!
            val index = anchorPrefixSum.binarySearch(rVal)
            val j = if (index >= 0) { // chromosome index
                index + 1
            } else {
                -index - 1
            }

            // Find anchor in 'allowed' methylome subset
            val sampledChr = anchorChrs[j] // same for anchor and full meth by design
            val achorMethChrData = background.anchorMethylome.data[sampledChr]
            val anchorStartOffsetIdx = achorMethChrData.size() + (rVal - anchorPrefixSum[j]).toInt()
            val anchorStartOffset = achorMethChrData[anchorStartOffsetIdx]

            // Map anchor id to id in full methylome
            val mappedAnchorStartOffsetIdx = background.mapAnchorOffsetToFull(sampledChr, anchorStartOffsetIdx)
            val fullMethChrData = background.fullMethylome.data[sampledChr]
            require(anchorStartOffset == fullMethChrData[mappedAnchorStartOffsetIdx]) {
                val mappedOffset = fullMethChrData[mappedAnchorStartOffsetIdx]
                "Mapping failed: $anchorStartOffset (idx:$anchorStartOffsetIdx) != $mappedOffset (idx: $mappedAnchorStartOffsetIdx)"
            }

            // [Not BED Tools]: just improvement so not to start a sampled region always directly in
            // a background position, but just intersect background position, called 'anchor'

            // XXX: Option 1
            // val startOffsetIdx = mappedAnchorStartOffsetIdx // XXX: like default BED Tools behavior

            // XXX: Option 2
            // shift at the beginning of chr, let's do not overflow + make that `mappedAnchorStartOffsetIdx` remained in the interval
            val startOffsetIdx = max(0, mappedAnchorStartOffsetIdx - rndStartShiftInCpGs)

            // Start offset final
            val startOffset = fullMethChrData[startOffsetIdx]

            // End offset:
            val endOffsetInclusiveIdx = startOffsetIdx + requiredCoveredPositionsNum - 1
            val endOffsetNextCoveredIdx = startOffsetIdx + requiredCoveredPositionsNum

            if (endOffsetInclusiveIdx < fullMethChrData.size()) {
                val endOffsetInclusiveWithShift = fullMethChrData[endOffsetInclusiveIdx] + endPositionShift
                val endOffset = if (endOffsetNextCoveredIdx < fullMethChrData.size()) {
                    // Due to shift could overlap with next position from background:
                    min(endOffsetInclusiveWithShift, fullMethChrData[endOffsetNextCoveredIdx])
                } else {
                    endOffsetInclusiveWithShift
                }

                val candidate = ChromosomeRange(startOffset, endOffset, sampledChr)

                if (candidate.endOffset <= candidate.chromosome.length) {
                    val isCandidateValid = if (maskedArea == null) {
                        // All candidates are valid, i.e., sampling with replacement and not masked genome
                        true
                    } else {
                        // Doesn't intersect already seen regions or initially masked genome
                        !candidate.intersectMap(maskedArea)
                    }

                    if (isCandidateValid) {
                        return candidate
                    }
                }
            }

            return null
        }

        fun calculatePrefixSum(background: BasePairCoverage): Pair<ArrayList<Chromosome>, LongArray> {
            // Here prefix sum is a sequence:
            //  * number of covered dots on chr1,
            //  * number of covered dots on chr1 + chr2,
            //  * number of covered dots on chr1 + chr2 + chr3 ....
            //  * ...
            // Used because coverage is grouped by chromosome and chrs covered differently,
            // so the probability of chromosome is weighted

            // XXX Make prefix sum list w/o duplicated values (i.e not covered chromosomes), otherwise shuffle doesn't work correctly
            val prefixSumList = arrayListOf<Long>()
            val backgroundChrs = arrayListOf<Chromosome>()
            var s = 0L
            for (chr in background.data.genomeQuery.get()) {
                val chrItems = background.data[chr].size()
                if (chrItems > 0) {
                    s += chrItems
                    prefixSumList.add(s)
                    backgroundChrs.add(chr)
                }
            }
            val prefixSum = prefixSumList.toLongArray()
            return Pair(backgroundChrs, prefixSum)
        }

        fun loadMethylomeCovBackground(
            backgroundRegionsPath: Path,
            zeroBasedBg: Boolean,
            gq: GenomeQuery
        ): BasePairCoverage {
            val methCovData = BasePairCoverage.loadFromTSV(
                gq, backgroundRegionsPath,
                offsetIsOneBased = !zeroBasedBg,
                progress = true
            )
            LOG.info("Background coverage: ${methCovData.depth.formatLongNumber()} offsets")

            return methCovData
        }

        fun filterMethylomeCovBackground(
            methCovData: BasePairCoverage,
            genomeAllowedAreaFilter: LocationsMergingList?,
            genomeMaskedAreaFilter: LocationsMergingList?,
        ): BasePairCoverage {
            // Here no need in intersecting background with `SharedSamplingOptions.truncateRangesToSpecificLocation`

            var anchorMethCovData = methCovData
            if (genomeAllowedAreaFilter != null) {
                LOG.info("Background coverage: Applying allowed area filter...")

                // Use meth covered in allowed areas to place a random anchor and then sample from full methylome
                // E.g. if DMR intersects a narrow open chromatin region,
                // we will not be able to place there a long candidate DMR, altough if there is enough adjacent CpG
                // to place intersecting DMRs.
                // Or we will be able to place DMR but it size could be very long because it goes only through CpG in
                // open chromantin
                anchorMethCovData = methCovData.filter(
                    genomeAllowedAreaFilter, progress = true, includeRegions = true, ignoreRegionsOnMinusStrand = true
                )

                LOG.info("Background coverage (allowed filter applied): ${anchorMethCovData.depth.formatLongNumber()} offsets of ${methCovData.depth.formatLongNumber()}")
            }

            if (genomeMaskedAreaFilter != null) {
                LOG.info("Background coverage: Applying masked area filter...")

                // These points are anchor points, if they are in `genomeMaskedAreaFilter`, then resulting regions
                // will be dismissed after sampling => we could remove them now.
                // But later checks still required because region started in valid location could intersect some masked.

                val beforeFilteringDepth = anchorMethCovData.depth
                anchorMethCovData = anchorMethCovData.filter(
                    genomeMaskedAreaFilter, progress = true, includeRegions = false, ignoreRegionsOnMinusStrand = true
                )

                LOG.info("Background coverage (w/o masked): ${anchorMethCovData.depth.formatLongNumber()} offsets of ${beforeFilteringDepth.formatLongNumber()}")
            }

            return anchorMethCovData
        }


        fun ensureInputRegionsMatchesBackgound(
            inputRegions: List<Location>,
            methCovData: BasePairCoverage,
            backgroundRegionsPath: Path
        ) {
            inputRegions.forEach { loc ->
                require(methCovData.getCoverage(loc) > 0) {
                    "Background $backgroundRegionsPath coverage should cover all input regions, but the " +
                            "region is missing in background: ${loc.toChromosomeRange()}"
                }
            }
            LOG.info("[OK] Regions matches methylome coverage")
        }

    }
}
class MethylomeSamplingBackground(
    val fullMethylome: BasePairCoverage,
    val anchorMethylome: BasePairCoverage,
) {

    val achorPrefixSum: LongArray
    val achorPrefixChrs: List<Chromosome>
    private val mapping: Map<Chromosome, TIntIntHashMap>?

    init {
        val res = RegionShuffleStatsFromMethylomeCoverage.calculatePrefixSum(anchorMethylome)
        achorPrefixChrs = res.first
        achorPrefixSum = res.second

        mapping = if (fullMethylome == anchorMethylome) {
            null
        } else {
            anchorMethylome.buildIndexMappingTo(fullMethylome)
        }
    }

    fun mapAnchorOffsetToFull(chr: Chromosome, anchorOffsetIdx: Int) = when (mapping) {
        null -> anchorOffsetIdx // same anchor methylome and full methylome
        else -> mapping[chr]!![anchorOffsetIdx]
    }

    fun calcExpCoverage(regions: List<ChromosomeRange>) = regions.map { fullMethylome.getCoverage(it.on(Strand.PLUS)) }.toIntArray()
}