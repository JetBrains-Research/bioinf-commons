package org.jetbrains.bio.gsea

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.genome.sampling.hasIntersection
import org.jetbrains.bio.util.formatLongNumber
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.ThreadLocalRandom
import kotlin.math.min

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
        intersectionFilter: LocationsList<out RangesList>?,
        genomeMaskedAreaPath: Path?,
        genomeAllowedAreaPath: Path?,
        samplingWithReplacement: Boolean,
        lengthCorrectionMethod: String?
    ): DataFrame {
        requireNotNull(backgroundPath) { "Background should be provided" }

        return doCalcStatistics(
            gq,
            loiInfos,
            outputFolderPath,
            metric,
            hypAlt,
            aSetIsRegions,
            mergeOverlapped,
            intersectionFilter,
            samplingWithReplacement = samplingWithReplacement,
            inputRegionsAndBackgroundProvider = { _ ->
                inputRegionsAndBackgroundProviderFun(inputRegionsPath, backgroundPath, zeroBasedBg, gq,
                    genomeMaskedAreaPath, genomeAllowedAreaPath)
            },
            loiOverlapWithBgFun = {loiFiltered, background ->
                loiFiltered.asLocationSequence().count {
                    background.getCoverage(it) > 0
                }
            },
            samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
                val res = shuffleChromosomeRanges(
                    genomeQuery,
                    regions,
                    background,
                    regionSetMaxRetries = regionSetMaxRetries,
                    singleRegionMaxRetries = singleRegionMaxRetries,
                    withReplacement = withReplacement,
                    candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                        lengthCorrectionMethod,
                        regions
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
            background: BasePairCoverage,
            regionSetMaxRetries: Int = 100,
            singleRegionMaxRetries: Int = 100,
            endPositionShift: Int = 2, // For if 'k' is index of latest CpG in candidate, use 'k+2' to get 'exclusive' bound after the latest CpG bounds
            withReplacement: Boolean,
            candidateFilterPredicate: ((ChromosomeRange) -> Double)?
        ): Pair<List<ChromosomeRange>, IntHistogram> {

            // Use coverage of input regions for sampling instead of region lengths in basepairs
            val expectedCoverages = regions.map { background.getCoverage(it.on(Strand.PLUS)) }.toIntArray()

            // Chromosomes have different amount of data => sampling should take into account that chromosomes probablities
            // should depend on # of observed data on each chromosome
            require(genomeQuery == background.data.genomeQuery) {
                "Background was made for different genome query. This query=[$genomeQuery]," +
                        " background genome query=[${background.data.genomeQuery}]"
            }

            val (backgroundChrs, prefixSum) = calculatePrefixSum(background)

            for (i in 1..regionSetMaxRetries) {
                val result = tryShuffle(
                    genomeQuery,
                    background,
                    backgroundChrs,
                    expectedCoverages,
                    prefixSum,
                    singleRegionMaxRetries,
                    endPositionShift,
                    withReplacement,
                    candidateFilterPredicate = candidateFilterPredicate
                )

                if (result != null) {
                    // TODO: collect stats how many retries is done
                    return result
                }
            }
            throw RuntimeException("Too many shuffle attempts, max limit is: $regionSetMaxRetries")
        }

        fun inputRegionsAndBackgroundProviderFun(
            inputRegionsPath: Path,
            backgroundPath: Path?,
            zeroBasedBg: Boolean,
            gq: GenomeQuery,
            genomeMaskedAreaPath: Path?,
            genomeAllowedAreaPath: Path?
        ): Pair<List<Location>, BasePairCoverage> {
            val inputRegions = loadInputRegions(inputRegionsPath, gq)
            val methylomeCov = loadMethylomeCovBackground(backgroundPath!!, zeroBasedBg, gq)
            val (inputRegionsFiltered, methylomeCovFiltered) = filterInputRegionsAndMethylomeCovBackground(
                inputRegions, methylomeCov, genomeMaskedAreaPath,
                genomeAllowedAreaPath, gq
            )
            ensureInputRegionsMatchesBackgound(inputRegionsFiltered, methylomeCovFiltered, backgroundPath)
            return inputRegionsFiltered to methylomeCovFiltered
        }

        private fun tryShuffle(
            genomeQuery: GenomeQuery,
            background: BasePairCoverage,
            backgroundChromosomes: List<Chromosome>,
            expectedCoverages: IntArray,
            prefixSum: LongArray, // should be w/o duplicated values
            singleRegionMaxRetries: Int = 100,
            endPositionShift: Int, // E.g. if 'k' is index of latest CpG in candidate, use 'k+2' (endPositionShift = 2) to get 'exclusive' bound after latest CpG bounds
            withReplacement: Boolean,
            candidateFilterPredicate: ((ChromosomeRange) -> Double)?
        ): Pair<List<ChromosomeRange>, IntHistogram>? {
            val maskedGenomeMap = when {
                withReplacement -> null
                // mask positions already used for sampled loci:
                else -> genomeMap<MutableList<Range>>(genomeQuery) {
                    ArrayList()
                }
            }
            val result = ArrayList<ChromosomeRange>(expectedCoverages.size)
            val attemptsHistogram = IntHistogram()

            val r = ThreadLocalRandom.current()

            val sum = prefixSum.last()

            for (i in expectedCoverages.indices) {
                val requiredCoveredPositionsNum = expectedCoverages[i]
                var nextRange: ChromosomeRange? = null

                var smalestPenalty: Double = Double.POSITIVE_INFINITY
                var bestCandidate: ChromosomeRange? = null

                var lastRangeTry = 0
                for (rangeTry in 1..singleRegionMaxRetries) {
                    lastRangeTry = rangeTry
                    val rVal = r.nextLong(sum)


                    var candidate: ChromosomeRange? = generateShuffledRegion(
                        prefixSum,
                        rVal,
                        backgroundChromosomes,
                        background,
                        requiredCoveredPositionsNum,
                        endPositionShift,
                        withReplacement =withReplacement,
                        maskedGenomeMap,
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

                // All available attempts done or result found:
                if (nextRange == null) {
                    // XXX: Cannot find IDEAL candidate => take the most probable
                    nextRange = bestCandidate
                }

                if (nextRange == null) {
                    // Cannot sample next range in given attempts number
                    return null
                }

                // success
                if (!withReplacement) {
                    maskedGenomeMap!![nextRange.chromosome].add(nextRange.toRange())
                }
                result.add(nextRange)
                attemptsHistogram.increment(lastRangeTry)
            }

            return result to attemptsHistogram
        }

        fun generateShuffledRegion(
            prefixSum: LongArray,
            rVal: Long,
            backgroundChromosomes: List<Chromosome>,
            background: BasePairCoverage,
            requiredCoveredPositionsNum: Int,
            endPositionShift: Int,
            withReplacement: Boolean,
            maskedGenomeMap: GenomeMap<MutableList<Range>>?,
        ): ChromosomeRange? {
            require(endPositionShift > 0)

            //NB: assume no duplicated values in prefix sum!!!!
            val index = prefixSum.binarySearch(rVal)
            val j = if (index >= 0) {
                index + 1
            } else {
                -index - 1
            }


            val jChr = backgroundChromosomes[j]
            val jChrBackground = background.data[jChr]
            val startOffsetIdx = jChrBackground.size() + (rVal - prefixSum[j]).toInt()

            val startOffset = jChrBackground[startOffsetIdx]

            val endOffsetInclusiveIdx = startOffsetIdx + requiredCoveredPositionsNum - 1
            val endOffsetNextCoveredIdx = startOffsetIdx + requiredCoveredPositionsNum

            if (endOffsetInclusiveIdx < jChrBackground.size()) {
                val endOffsetInclusiveWithShift = jChrBackground[endOffsetInclusiveIdx] + endPositionShift
                val endOffset = if (endOffsetNextCoveredIdx < jChrBackground.size()) {
                    // Due to shift could overlap with next position from background:
                    min(endOffsetInclusiveWithShift, jChrBackground[endOffsetNextCoveredIdx])
                } else {
                    endOffsetInclusiveWithShift
                }

                val candidate = ChromosomeRange(startOffset, endOffset, jChr)

                if (candidate.endOffset <= candidate.chromosome.length) {
                    val isCandidateValid = if (withReplacement) {
                        // N/A
                        true
                    } else {
                        // Doesn't intersect already seen regions
                        !hasIntersection(maskedGenomeMap!!, candidate)
                    }

                    if (isCandidateValid) {
                        return candidate
                    }
                }
            }

            return null
        }

        fun calculatePrefixSum(background: BasePairCoverage): Pair<ArrayList<Chromosome>, LongArray> {
            // XXX Make prefix sum list w/o duplicated values (i.e not covered chromosomes), otherwise shuffule doesn't work correctly
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

        fun loadInputRegions(
            inputRegionsPath: Path,
            gq: GenomeQuery
        ): List<Location> {
            val inputRegions = readLocationsIgnoringStrand(inputRegionsPath, gq).first
            LOG.info("Input regions: ${inputRegions.size.formatLongNumber()} regions")

            return inputRegions
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

        fun filterInputRegionsAndMethylomeCovBackground(
            inputRegions: List<Location>,
            methCovData: BasePairCoverage,
            genomeMaskedAreaPath: Path?,
            genomeAllowedAreaPath: Path?,
            gq: GenomeQuery
        ): Pair<List<Location>, BasePairCoverage> {
            val allowedGenomeFilter = makeAllowedRegionsFilter(
                genomeMaskedAreaPath, genomeAllowedAreaPath, gq
            )

            @Suppress("FoldInitializerAndIfToElvis")
            if (allowedGenomeFilter == null) {
                return inputRegions to methCovData
            }

            LOG.info("Applying allowed regions filters...")

            val allowedMethCovData = methCovData.filter(allowedGenomeFilter, progress = true, includeRegions = true, ignoreRegionsOnMinusStrand = true)
            LOG.info("Background coverage (all filters applied): ${allowedMethCovData.depth.formatLongNumber()} offsets of ${methCovData.depth.formatLongNumber()}")

            val allowedInputRegions = inputRegions.filter { allowedGenomeFilter.includes(it) }
            LOG.info("Input regions (all filters applied): ${allowedInputRegions.size.formatLongNumber()} regions of ${inputRegions.size.formatLongNumber()}")

            return allowedInputRegions to allowedMethCovData
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