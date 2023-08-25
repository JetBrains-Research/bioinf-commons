package org.jetbrains.bio.gsea

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.*
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.genome.sampling.hasIntersection
import org.jetbrains.bio.util.formatLongNumber
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.ThreadLocalRandom
import kotlin.math.min

class RegionShuffleStatsFromMethylomeCoverage(
    genome: Genome,
    simulationsNumber: Int,
    chunkSize: Int,
    maxRetries: Int
) : RegionShuffleStats(genome, simulationsNumber, chunkSize, maxRetries) {

    override fun calcStatistics(
        srcRegionsPath: Path,
        backgroundPath: Path?,
        regionLabelAndLociToTest: List<Pair<String, LocationsList<out RangesList>>>,
        outputFolderPath: Path?,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis,
        aSetIsLoi: Boolean,
        mergeOverlapped: Boolean,
        intersectionFilter: LocationsList<out RangesList>?,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        addLoiToBg: Boolean,
        samplingWithReplacement: Boolean
    ): DataFrame {
        requireNotNull(backgroundPath) { "Background should be provided" }

        return doCalcStatistics(
            regionLabelAndLociToTest,
            outputFolderPath,
            metric,
            hypAlt,
            aSetIsLoi,
            mergeOverlapped,
            intersectionFilter,
            srcLociAndBackgroundProvider = { gq ->
                loadSrcLociAndMethylomeCovBackground(
                    srcRegionsPath, backgroundPath, genomeMaskedLociPath, genomeAllowedLociPath, gq
                )
            },
            samplingFun = { gq, regions, background, maxRetries, withReplacement ->
                shuffleChromosomeRanges(
                    gq,
                    regions,
                    background,
                    maxRetries = maxRetries,
                    withReplacement = withReplacement
                ).map { it.on(Strand.PLUS) }
            }
        )
    }

    fun loadSrcLociAndMethylomeCovBackground(
        srcRegionsPath: Path,
        backgroundRegionsPath: Path,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        gq: GenomeQuery
    ): Pair<List<Location>, BasePairCoverage> {
        val sourceLoci = readLocationsIgnoringStrand(srcRegionsPath, gq)
        LOG.info("Source loci: ${sourceLoci.size} regions")

        val methCovData = BasePairCoverage.loadFromTSV(
            gq, backgroundRegionsPath,
            offsetIsOneBased = true,  // TODO: [make cmdline option] methylome coverage table is 1-based, convert to 0-based!
            progress = true
        )
        LOG.info("Background coverage: ${methCovData.depth.formatLongNumber()} offsets")

        sourceLoci.forEach {
            require(methCovData.getCoverage(it) > 0) {
                "Background $backgroundRegionsPath coverage should cover all input loci, but the " +
                        "loci is missing in background: ${it.toChromosomeRange()}"
            }
        }
        LOG.info("[OK] Loci matches methylome coverage")

        val allowedGenomeFilter = makeAllowedRegionsFilter(
            genomeMaskedLociPath, genomeAllowedLociPath, gq
        )

        val (allowedMethCovData, allowedSourceLoci) = if (allowedGenomeFilter == null) {
            methCovData to sourceLoci
        } else {
            LOG.info("Applying allowed regions filters...")

            val allowedMethCovData = methCovData.filterExclude(allowedGenomeFilter, progress = true)
            LOG.info("Background coverage (all restrictions applied): ${allowedMethCovData.depth.formatLongNumber()} offsets")

            val allowedSourceLoci = sourceLoci.filter { allowedGenomeFilter.includes(it) }
            LOG.info("Source loci (all restrictions applied): ${allowedSourceLoci.size.formatLongNumber()} regions")

            allowedMethCovData to allowedSourceLoci
        }
        return allowedSourceLoci to allowedMethCovData
    }


    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStatsFromBedCoverage::class.java)

        fun shuffleChromosomeRanges(
            genomeQuery: GenomeQuery,
            regions: List<ChromosomeRange>,
            background: BasePairCoverage,
            maxRetries: Int = 100,
            endPositionShift: Int = 2, // For if 'k' is index of latest CpG in candidate, use 'k+2' to get 'exclusive' bound after latest CpG bounds
            withReplacement: Boolean
        ): List<ChromosomeRange> {

            // Use coverage of input regions for sampling instead of region lengths in basepairs
            val expectedCoverages = regions.map { background.getCoverage(it.on(Strand.PLUS)) }.toIntArray()

            // Chromosome have different amount of data => sampling should take into account that chromosomes probablities
            // should depend on # of observed data on each chromosome
            require(genomeQuery == background.data.genomeQuery) {
                "Background was made for different genome query. This query=[$genomeQuery]," +
                        " background genome query=[${background.data.genomeQuery}]"
            }

            val (backgroundChrs, prefixSum) = calculatePrefixSum(background)

            for (i in 1..maxRetries) {
                val result = tryShuffle(
                    genomeQuery,
                    background,
                    backgroundChrs,
                    expectedCoverages,
                    prefixSum,
                    maxRetries,
                    endPositionShift,
                    withReplacement
                )
                if (result != null) {
                    // TODO: collect stats how many retries is done
                    return result
                }
            }
            throw RuntimeException("Too many shuffle attempts")
        }

        private fun tryShuffle(
            genomeQuery: GenomeQuery,
            background: BasePairCoverage,
            backgroundChromosomes: List<Chromosome>,
            expectedCoverages: IntArray,
            prefixSum: LongArray, // should be w/o duplicated values
            maximalReTries: Int = 100,
            endPositionShift: Int, // E.g. if 'k' is index of latest CpG in candidate, use 'k+2' (endPositionShift = 2) to get 'exclusive' bound after latest CpG bounds
            withReplacement: Boolean
        ): List<ChromosomeRange>? {
            val maskedGenomeMap = when {
                withReplacement -> null
                // mask positions already used for sampled loci:
                else -> genomeMap<MutableList<Range>>(genomeQuery) {
                    ArrayList()
                }
            }
            val result = ArrayList<ChromosomeRange>(expectedCoverages.size)

            val r = ThreadLocalRandom.current()

            val sum = prefixSum.last()

            for (i in expectedCoverages.indices) {
                val requiredCoveredPositionsNum = expectedCoverages[i]
                var nextRange: ChromosomeRange? = null

                for (rangeTry in 1..maximalReTries) {
                    val rVal = r.nextLong(sum)

                    val candidate: ChromosomeRange? = generateShuffledRegion(
                        prefixSum,
                        rVal,
                        backgroundChromosomes,
                        background,
                        requiredCoveredPositionsNum,
                        endPositionShift,
                        withReplacement,
                        maskedGenomeMap
                    )

                    if (candidate != null) {
                        nextRange = candidate
                        break
                    }
                    // else retry
                }

                if (nextRange == null) {
                    // Cannot sample next range in given attempts number
                    return null
                }

                result.add(nextRange)
            }

            return result
        }

        fun generateShuffledRegion(
            prefixSum: LongArray,
            rVal: Long,
            backgroundChromosomes: List<Chromosome>,
            background: BasePairCoverage,
            requiredCoveredPositionsNum: Int,
            endPositionShift: Int,
            withReplacement: Boolean,
            maskedGenomeMap: GenomeMap<MutableList<Range>>?
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
                        // success
                        if (!withReplacement) {
                            maskedGenomeMap!![candidate.chromosome].add(candidate.toRange())
                        }
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

    }
}