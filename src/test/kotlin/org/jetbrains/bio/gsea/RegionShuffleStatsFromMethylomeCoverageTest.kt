package org.jetbrains.bio.gsea

import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.util.withTempBedFile
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.fail

class RegionShuffleStatsFromMethylomeCoverageTest {
    val COVERAGE_EX = listOf(
        chr1 to 1,
        chr1 to 5,
        chr1 to 9,
        chr1 to 10,
        chr1 to 20,
        chr1 to 22,
        chr1 to 24,
        chr1 to 25,
        chr1 to 26,
        chr1 to 28,
        chr2 to 3,
    )
    val COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR = listOf(
        chr3 to 1,
        chr3 to 5,
        chr3 to 9,
        chr3 to 10,
        chr3 to 20,
        chr3 to 22,
        chr3 to 24,
        chr3 to 25,
        chr3 to 26,
        chr3 to 28,
    )

    @Test
    fun shuffleChromosomeRangesWithReplacement() {
        val background = createBasePairCoverage(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, gq)
        val regions = listOf(
            ChromosomeRange(2, 7, chr3),
            ChromosomeRange(2, 27, chr3),
            ChromosomeRange(10, 20, chr3),
            ChromosomeRange(22, 27, chr3),
        )
        val newRegions = RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
            gq, regions, background, singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = true,
            candidateFilterPredicate = null
        ).first

        assertEquals(regions.size, newRegions.size)
        assertEquals(
            regions.map { background.getCoverage(it.on(Strand.PLUS)) }.sorted(),
            newRegions.map { background.getCoverage(it.on(Strand.PLUS)) }.sorted(),
        )
    }

    @Test
    fun shuffleChromosomeRangesWithoutReplacement() {
        val background = createBasePairCoverage(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, gq)
        val regions = listOf(
            ChromosomeRange(2, 7, chr3),
            ChromosomeRange(10, 20, chr3),
            ChromosomeRange(22, 27, chr3),
        )
        val shuffled = RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
            gq, regions, background, singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
            candidateFilterPredicate = null
        ).first

        assertEquals(regions.size, shuffled.size)
        assertEquals(
            regions.map { background.getCoverage(it.on(Strand.PLUS)) }.sorted(),
            shuffled.map { background.getCoverage(it.on(Strand.PLUS)) }.sorted(),
        )

        for (i1 in shuffled.indices) {
            for (i2 in shuffled.indices) {
                if (i1 == i2) {
                    continue
                }
                if (shuffled[i1].chromosome == shuffled[i2].chromosome &&
                    shuffled[i1].toRange() intersects shuffled[i2].toRange()
                ) {
                    fail("Regions must not intersect")
                }
            }
        }
    }

    @Test
    fun shuffleChromosomeRangesWithoutReplacementNotPossible() {
        val background = createBasePairCoverage(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, gq)
        val regions = listOf(
            ChromosomeRange(2, 27, chr3),
            ChromosomeRange(22, 27, chr3),
        )
        assertFailsWith(RuntimeException::class, message = "Too many shuffle attempts") {
            RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
                gq, regions, background, singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
                candidateFilterPredicate = null
            ).first
        }
    }

    @Test
    fun shuffleChromosomeRangesWrongBgGenome() {
        val background = createBasePairCoverage(
            COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, GenomeQuery(gq.genome, "chr3")
        )
        val regions = listOf(
            ChromosomeRange(2, 27, chr3),
            ChromosomeRange(22, 27, chr3),
        )
        assertFailsWith(
            IllegalArgumentException::class,
            message = "Background was made for different genome query. This query=[to1], background genome query=[to1[chr3]]"
        ) {
            RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
                gq, regions, background, singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
                candidateFilterPredicate = null
            ).first
        }
    }

    @Test
    fun generateShuffledRegionWithReplace() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 3, ChromosomeRange(1, 10, chr1), -1),
                Element(rVal = 2, nPos = 1, ChromosomeRange(9, 10, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 8, nPos = 2, ChromosomeRange(26, 30, chr1), -1),
                Element(rVal = 8, nPos = 3, null, -1),
            ),
            endPositionShift = 2,
            withReplacement = true
        )
    }

    @Test
    fun generateShuffledRegionWithReplaceWithEmptyChrs() {
        checkShuffledRegion(
            COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr3), -1),
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr3), -1),
                Element(rVal = 0, nPos = 3, ChromosomeRange(1, 10, chr3), -1),
                Element(rVal = 2, nPos = 1, ChromosomeRange(9, 10, chr3), -1),
                Element(rVal = 9, nPos = 1, ChromosomeRange(28, 30, chr3), -1),
                Element(rVal = 9, nPos = 2, null, -1),
            ),
            endPositionShift = 2,
            withReplacement = true
        )
    }


    @Test
    fun generateShuffledRegionWithoutReplace() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                // 1st success:
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), 1),
                // same interval
                Element(rVal = 0, nPos = 1, null, 1),
                // container interval
                Element(rVal = 0, nPos = 3, null, 1),

                // 2nd success:
                Element(rVal = 2, nPos = 4, ChromosomeRange(9, 24, chr1), 2),
                // intersect left:
                Element(rVal = 1, nPos = 2, null, 2),
                // intersect right:
                Element(rVal = 4, nPos = 2, null, 2),
                // out of bounds:
                Element(rVal = 8, nPos = 3, null, 2),

                // 3rd success:
                Element(rVal = 8, nPos = 2, ChromosomeRange(26, 30, chr1), 3),

                // 4th success:
                Element(rVal = 10, nPos = 1, ChromosomeRange(3, 5, chr2), 4),
            ),
            endPositionShift = 2,
            withReplacement = false
        )
    }

    @Test
    fun generateShuffledRegionWithCustomShift() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 2, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 9, nPos = 1, ChromosomeRange(28, 29, chr1), -1),
            ),
            endPositionShift = 1,
            withReplacement = true
        )

        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 4, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 9, nPos = 1, ChromosomeRange(28, 31, chr1), -1),
            ),
            endPositionShift = 3,
            withReplacement = true
        )

        assertFailsWith(IllegalArgumentException::class) {
            checkShuffledRegion(
                COVERAGE_EX,
                listOf(
                    Element(rVal = 0, nPos = 1, ChromosomeRange(1, 1, chr1), -1),
                ),
                endPositionShift = 0,
                withReplacement = true
            )
        }
    }

    @Test
    fun generateShuffledRegionLenExceedBounds() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 8, nPos = 3, null, -1),
                Element(rVal = 10, nPos = 2, null, -1),
            ),
            endPositionShift = 2,
            withReplacement = true
        )
    }

    private fun checkShuffledRegion(
        coverageContent: List<Pair<Chromosome, Int>>,
        elementsInfo: List<Element>,
        withReplacement: Boolean,
        endPositionShift: Int
    ) {
        val maskedGenomeMap = when {
            withReplacement -> null
            // mask positions already used for sampled loci:
            else -> genomeMap<MutableList<Range>>(gq) {
                ArrayList()
            }
        }

        val background = createBasePairCoverage(coverageContent, gq)

        val (backgroundChrs, prefixSum) = RegionShuffleStatsFromMethylomeCoverage.calculatePrefixSum(
            background
        )

        val results = elementsInfo.map { info ->
            require((info.rVal >= 0) && (info.rVal < prefixSum.last()))


            val actualRange = RegionShuffleStatsFromMethylomeCoverage.generateShuffledRegion(
                prefixSum,
                info.rVal.toLong(),
                backgroundChrs,
                background,
                info.nPos,
                endPositionShift,
                withReplacement = withReplacement,
                maskedGenomeMap
            )
            val actualMaskedPositions = maskedGenomeMap?.let { m ->
                gq.get().sumOf { m.get(it).size }
            } ?: -1
            info.copy(actualRange = actualRange, actualMaskedPositions = actualMaskedPositions)
        }

        results.forEach { r ->
            assertEquals(r.expectedRange, r.actualRange, message = "Mismatch in : rVal=${r.rVal}, nPos=${r.nPos}")
            assertEquals(
                r.expectedMaskedPositions,
                r.actualMaskedPositions,
                message = "Mismatch in : rVal=${r.rVal}, nPos=${r.nPos}"
            )
            if (r.actualRange != null) {
                val actualNPos = background.getCoverage(r.actualRange.on(Strand.PLUS))
                assertEquals(
                    r.nPos,
                    actualNPos,
                    message = "Mismatch in : rVal=${r.rVal}, nPos=${r.nPos}"
                )
            }
        }
    }

    private fun createBasePairCoverage(
        coverageContent: List<Pair<Chromosome, Int>>,
        genomeQuery: GenomeQuery
    ): BasePairCoverage {
        val builder = BasePairCoverage.builder(genomeQuery, false)
        coverageContent.forEach { (chr, offset) ->
            builder.process(chr, offset)
        }
        return builder.build(true)
    }

    @Test
    fun calculatePrefixSum() {
        checkCalculatePrefixSum(COVERAGE_EX, 0, chr1, 0)
        checkCalculatePrefixSum(COVERAGE_EX, 3, chr1, 0)
        checkCalculatePrefixSum(COVERAGE_EX, 9, chr1, 0)
        assertFailsWith(IndexOutOfBoundsException::class, message = "Index 1 out of bounds for length 1") {
            checkCalculatePrefixSum(COVERAGE_EX, 11, chr2, 1)
        }
    }

    @Test
    fun calculatePrefixSumWithNonCovChrs() {
        checkCalculatePrefixSum(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, 0, chr3, 0)
        checkCalculatePrefixSum(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, 9, chr3, 0)
        assertFailsWith(IndexOutOfBoundsException::class, message = "Index 1 out of bounds for length 1") {
            checkCalculatePrefixSum(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, 10, chr3, 0)
        }
    }

    @Test
    fun loadInputRegionsAndMethylomeCovBackground() {
        val lociGenomeAllowed = listOf(
            Location(10, 100, RegionShuffleStatsTest.chr1),
        )
        val lociGenomeMasked = listOf(
            Location(28, 80, RegionShuffleStatsTest.chr1),
        )

        val lociInputRegions = listOf(
            Location(20, 22, RegionShuffleStatsTest.chr1),
            Location(24, 27, RegionShuffleStatsTest.chr1),
            // intersects ignored area
            Location(9, 14, RegionShuffleStatsTest.chr1),
            Location(27, 60, RegionShuffleStatsTest.chr1),
            // masked:
            Location(62, 70, RegionShuffleStatsTest.chr1),
            // outside allowed:
            Location(2, 5, RegionShuffleStatsTest.chr1),
            Location(110, 200, RegionShuffleStatsTest.chr1),
        )

        val (inputRegionsFiltered, cov) = withTempBedFile(lociGenomeAllowed) { genomeAllowedLociPath ->
            withTempBedFile(lociGenomeMasked) { genomeMaskedLociPath ->
                withTempBedFile(lociInputRegions) { inputRegionsPath ->
                    withTempFile("bgcov", ".tsv") { bgRegionsPath ->
                        createBasePairCoverage(COVERAGE_EX, gq).saveToTSV(bgRegionsPath)

                        RegionShuffleStatsFromMethylomeCoverage.inputRegionsAndBackgroundProviderFun(
                            inputRegionsPath,
                            bgRegionsPath,
                            zeroBasedBg = false,
                            gq,
                            genomeMaskedLociPath,
                            genomeAllowedLociPath
                        )
                    }
                }
            }
        }

        assertEquals(
            listOf(
                Location(20,22, chr1),
                Location(24,27, chr1)
            ).sorted(),
            inputRegionsFiltered.toList().sorted()
        )

        assertEquals(6, cov.depth)
        assertEquals("{10, 20, 22, 24, 25, 26}", cov.data[chr1].toString())
    }

    private fun checkCalculatePrefixSum(
        coverageContent: List<Pair<Chromosome, Int>>,
        rVal: Long,
        chr: Chromosome,
        idx: Int
    ) {
        val background = createBasePairCoverage(coverageContent, gq)

        val (backgroundChrs, prefixSum) = RegionShuffleStatsFromMethylomeCoverage.calculatePrefixSum(
            background
        )
        val index = prefixSum.binarySearch(rVal)
        val j = if (index >= 0) {
            index + 1
        } else {
            -index - 1
        }
        assertEquals(chr, backgroundChrs[j])
        assertEquals(idx, j)
    }


    data class Element(
        val rVal: Int,
        val nPos: Int,
        val expectedRange: ChromosomeRange?,
        val expectedMaskedPositions: Int = -1,
        val actualRange: ChromosomeRange? = null,
        val actualMaskedPositions: Int = -1,
    )

    companion object {
        internal val gq: GenomeQuery = GenomeQuery(Genome["to1"])
        internal val chr1: Chromosome = gq.get()[0]
        internal val chr2: Chromosome = gq.get()[1]
        internal val chr3: Chromosome = gq.get()[2]
    }
}