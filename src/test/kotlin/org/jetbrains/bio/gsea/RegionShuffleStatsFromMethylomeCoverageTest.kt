package org.jetbrains.bio.gsea

import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.BasePairCoverage
import org.jetbrains.bio.genome.sampling.createMaskedArea
import org.jetbrains.bio.gsea.EnrichmentInLoi.processInputRegions
import org.jetbrains.bio.util.withTempBedFile
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.math.floor
import kotlin.random.Random
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
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
    val COVERAGE_EX2 = listOf(
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
        chr1 to 99,
    )
    val COVERAGE_CGI = listOf(
        chr1 to 1,
        chr1 to 9,
        chr1 to 20,
        chr1 to 32,
        chr1 to 42,
        chr1 to 52,
        chr1 to 100,
        chr1 to 102,
        chr1 to 104,
        chr1 to 106,
        chr1 to 108,
        chr1 to 110,
        chr1 to 112,
        chr1 to 150,
        chr1 to 170,
        chr1 to 210,
        chr1 to 412,
        chr1 to 414,
        chr1 to 416,
        chr1 to 420,
        chr1 to 422,
        chr1 to 424,
        chr1 to 430,
        chr1 to 500,
        chr1 to 900,
        chr1 to 980,
        chr1 to 1000,
        chr2 to 3,
        chr2 to 30,
        chr2 to 100,
        chr2 to 200,
        chr2 to 300,
        chr2 to 600,
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
            gq, regions,
            MethylomeSamplingBackground(background, background),
            singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = true,
            candidateFilterPredicate = null, genomeMaskedAreaFilter = null
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
            gq, regions, MethylomeSamplingBackground(background, background),
            singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
            candidateFilterPredicate = null, genomeMaskedAreaFilter = null
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
        Tests.assertThrowsWithMessage(
            expectedThrowable = RuntimeException::class.java,
            message = "Too many shuffle attempts, max limit is: 4"
        ) {
            RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
                gq, regions, MethylomeSamplingBackground(background, background),
                singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
                candidateFilterPredicate = null, genomeMaskedAreaFilter = null
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
        Tests.assertThrowsWithMessage(
            expectedThrowable = IllegalArgumentException::class.java,
            message = "Background was made for different genome query. This query=[to1], background genome query=[to1[chr3]]"
        ) {
            RegionShuffleStatsFromMethylomeCoverage.shuffleChromosomeRanges(
                gq, regions, MethylomeSamplingBackground(background, background),
                singleRegionMaxRetries = 4, regionSetMaxRetries = 4, endPositionShift = 2, withReplacement = false,
                candidateFilterPredicate = null, genomeMaskedAreaFilter = null
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
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.0
        )
    }

    @Test
    fun generateShuffledRegionWithShift() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 3, ChromosomeRange(1, 10, chr1), -1),
                Element(rVal = 2, nPos = 1, ChromosomeRange(9, 10, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 8, nPos = 2, ChromosomeRange(25, 28, chr1), -1),
                Element(rVal = 8, nPos = 3, ChromosomeRange(24, 28, chr1), -1),
            ),
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.8
        )

        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 3, ChromosomeRange(1, 10, chr1), -1),
                Element(rVal = 2, nPos = 1, ChromosomeRange(9, 10, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 8, nPos = 2, ChromosomeRange(25, 28, chr1), -1),
                Element(rVal = 8, nPos = 3, ChromosomeRange(25, 30, chr1), -1),
            ),
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.5
        )

        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), -1),
                Element(rVal = 0, nPos = 3, ChromosomeRange(1, 10, chr1), -1),
                Element(rVal = 2, nPos = 1, ChromosomeRange(5, 7, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(22, 24, chr1), -1),
                Element(rVal = 8, nPos = 2, ChromosomeRange(24, 26, chr1), -1),
                Element(rVal = 8, nPos = 3, ChromosomeRange(22, 26, chr1), -1),
            ),
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 1.0
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
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.0
        )
    }


    @Test
    fun generateShuffledRegionWithoutReplace() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                // 1st success:
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 3, chr1), 1), //TODO: 0
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
            withReplacement = false,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.0
        )
    }

    @Test
    fun generateShuffledRegionWithCustomEndShift() {
        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 2, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 9, nPos = 1, ChromosomeRange(28, 29, chr1), -1),
            ),
            withReplacement = true,
            endPositionShift = 1,
            rndStartShiftInCpGsFractionFraction = 0.0
        )

        checkShuffledRegion(
            COVERAGE_EX,
            listOf(
                Element(rVal = 0, nPos = 1, ChromosomeRange(1, 4, chr1), -1),
                Element(rVal = 6, nPos = 1, ChromosomeRange(24, 25, chr1), -1),
                Element(rVal = 9, nPos = 1, ChromosomeRange(28, 31, chr1), -1),
            ),
            withReplacement = true,
            endPositionShift = 3,
            rndStartShiftInCpGsFractionFraction = 0.0
        )

        Tests.assertThrowsWithMessage(IllegalArgumentException::class.java) {
            checkShuffledRegion(
                COVERAGE_EX,
                listOf(
                    Element(rVal = 0, nPos = 1, ChromosomeRange(1, 1, chr1), -1),
                ),
                withReplacement = true,
                endPositionShift = 0,
                rndStartShiftInCpGsFractionFraction = 0.0
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
            withReplacement = true,
            endPositionShift = 2,
            rndStartShiftInCpGsFractionFraction = 0.0
        )
    }

    @Test
    fun testTryShuffleRegionsWithMaskedRegionWithoutReplacement() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            maskedGenome = listOf(
                Location(23, 26, chr1),
                Location(1, 2, chr2),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[1, 22), chr1:[26, 28)], {14=1, 2=1, 1=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 20
        )
    }

    @Test
    fun testTryShuffleRegionsWithMaskedRegionWithReplacement() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            maskedGenome = listOf(
                Location(23, 26, chr1),
                Location(1, 2, chr2),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[1, 22), chr1:[9, 10)], {2=1, 1=2}",
            withReplacement = true,
            endPositionShift = 2,
            singleRegionMaxRetries = 20
        )
    }

    @Test
    fun testTryShuffleRegionsWithMaskedRegionNotPossible() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            maskedGenome = listOf(
                Location(23, 26, chr1),
                Location(1, 2, chr2),
            ),
            expectedResult = null,
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 13
        )
    }

    @Test
    fun testTryShuffleRegionsNotCovered() {
        Tests.assertThrowsWithMessage(
            java.lang.IllegalArgumentException::class.java,
            "Uncovered regions not supported, i=1"
        ) {
            checkTryShuffledRegions(
                COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR,
                listOf(
                    ChromosomeRange(1, 3, chr3),
                    ChromosomeRange(9, 24, chr1),
                ),
                expectedResult = "",
                withReplacement = false,
                endPositionShift = 2,
            )
        }
    }
    @Test
    fun testTryShuffleRegionsFromSingleChr() {
        checkTryShuffledRegions(
            COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR,
            listOf(
                ChromosomeRange(1, 3, chr3),
                ChromosomeRange(9, 24, chr3),
                ChromosomeRange(28, 29, chr3),
            ),
            expectedResult = "[chr3:[26, 28), chr3:[10, 25), chr3:[25, 26)], {2=2, 1=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 2
        )
    }

    @Test
    fun testTryShuffleRegionsFromChrDataZeroPos() {
        checkTryShuffledRegions(
            coverageContent = listOf(chr3 to 1,),
            inputRegions = listOf(ChromosomeRange(1, 3, chr3),),
            expectedResult = "[chr3:[1, 3)], {1=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 2
        )
    }
    @Test
    fun testTryShuffleRegionsWithEndPosShift1() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = "[chr2:[3, 4), chr1:[9, 25), chr1:[1, 2)], {3=1, 1=2}",
            withReplacement = false,
            endPositionShift = 1,
            singleRegionMaxRetries = 3
        )
    }
    @Test
    fun testTryShuffleRegionsWithEndPosShift2() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[9, 25), chr1:[1, 3)], {3=1, 1=2}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 3
        )
    }

    @Test
    fun testTryShuffleRegionsWithEndPosShift5() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = "[chr2:[3, 8), chr1:[9, 25), chr1:[1, 5)], {3=1, 1=2}",
            withReplacement = false,
            endPositionShift = 5,
            singleRegionMaxRetries = 3
        )
    }

    @Test
    fun testTryShuffleRegionsWithoutReplacement() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(9, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[22, 28), chr1:[9, 10)], {1=3}",
            withReplacement = false,
            endPositionShift = 2,
        )
    }

    @Test
    fun testTryShuffleRegionsWithoutReplacementNotPossible() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(9, 24, chr1),
                ChromosomeRange(24, 28, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = null,
            withReplacement = false,
            endPositionShift = 2,
        )
    }
    @Test
    fun testTryShuffleRegionsWithReplacement() {
        checkTryShuffledRegions(
            COVERAGE_EX,
            listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(9, 24, chr1),
                ChromosomeRange(24, 28, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[22, 28), chr1:[1, 10), chr1:[9, 10)], {1=4}",
            withReplacement = true,
            endPositionShift = 2,
        )
    }

    @Test
    fun testTryShuffleRegionsWithLenCorrectionNo() {
        val regions = listOf(
            ChromosomeRange(102, 112, chr1),
            ChromosomeRange(20, 22, chr1),
            ChromosomeRange(414, 426, chr1),
        )

        checkTryShuffledRegions(
            COVERAGE_CGI,
            regions,
            candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                null,
                regions,  Random(22)
            ),
            expectedResult = "[chr1:[424, 982), chr1:[414, 416), chr1:[1, 44)], {2=1, 1=2}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 3
        )
    }

    @Test
    fun testTryShuffleRegionsWithLenCorrectionFixed16() {
        val regions = listOf(
            ChromosomeRange(102, 112, chr1),
            ChromosomeRange(20, 22, chr1),
            ChromosomeRange(414, 426, chr1),
        )

        checkTryShuffledRegions(
            COVERAGE_CGI,
            regions,
            candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                "16",
                regions, Random(22)
            ),
            expectedResult = "[chr1:[104, 114), chr2:[100, 102), chr1:[414, 426)], {10=1, 5=1, 1=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 10
        )
    }

    @Test
    fun testTryShuffleRegionsWithLenCorrectionFixed100() {
        val regions = listOf(
            ChromosomeRange(102, 112, chr1),
            ChromosomeRange(20, 22, chr1),
            ChromosomeRange(414, 426, chr1),
        )

        checkTryShuffledRegions(
            COVERAGE_CGI,
            regions,
            candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                "100",
                regions, Random(22)
            ),
            expectedResult = "[chr1:[1, 44), chr1:[110, 112), chr1:[52, 108)], {4=1, 2=1, 1=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 10
        )
    }


    @Test
    fun testTryShuffleRegionsWithLenCorrectionDist() {
        val regions = listOf(
            ChromosomeRange(102, 112, chr1),
            ChromosomeRange(20, 22, chr1),
            ChromosomeRange(414, 426, chr1),
        )

        // TODO: why not stable?
        checkTryShuffledRegions(
            COVERAGE_CGI,
            regions,
            candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                "dist",
                regions,  Random(22)
            ),
            expectedResult = "[chr1:[104, 114), chr1:[100, 102), chr1:[414, 426)], {9=1, 5=1, 2=1}",
            withReplacement = false,
            endPositionShift = 2,
            singleRegionMaxRetries = 10
        )
    }

    @Test
    fun testTryShuffleRegionsWithLenCorrectionDistWithReplacement() {
        val regions = listOf(
            ChromosomeRange(102, 112, chr1),
            ChromosomeRange(20, 22, chr1),
            ChromosomeRange(414, 426, chr1),
        )

        checkTryShuffledRegions(
            COVERAGE_CGI,
            regions,
            candidateFilterPredicate = SamplingMethylationValidation.getCandidateFilterPredicate(
                "dist",
                regions, Random(22)
            ),
            expectedResult = "[chr1:[104, 114), chr1:[100, 102), chr1:[104, 114)], {10=1, 5=1, 2=1}",
            withReplacement = true,
            endPositionShift = 2,
            singleRegionMaxRetries = 10
        )
    }

    private fun checkShuffledRegion(
        coverageContent: List<Pair<Chromosome, Int>>,
        elementsInfo: List<Element>,
        withReplacement: Boolean,
        endPositionShift: Int,
        rndStartShiftInCpGsFractionFraction: Double,
        genomeMaskedAreaFilter: LocationsMergingList? = null
    ) {
        // Trick: used one storage for masked & seen intervals that is forbidden to intersect
        val maskedArea: GenomeMap<MutableList<Range>>? = createMaskedArea(genomeMaskedAreaFilter, withReplacement, gq)

        val fullMethylome = createBasePairCoverage(coverageContent, gq)
        val background = MethylomeSamplingBackground(fullMethylome, fullMethylome)

        val prefixSum = background.achorPrefixSum
        val results = elementsInfo.map { info ->
            require((info.rVal >= 0) && (info.rVal < prefixSum.last()))

            val rndStartShiftInCpGsFraction = floor(info.nPos * rndStartShiftInCpGsFractionFraction).toInt()
            val nextRange = RegionShuffleStatsFromMethylomeCoverage.generateShuffledRegion(
                info.rVal.toLong(),
                background,
                info.nPos,
                endPositionShift,
                rndStartShiftInCpGsFraction,
                maskedArea,
            )

            // Emulation of `tryShuffle` to check the sequence of accepted regions and rejected due to masking:
            if (!withReplacement && nextRange != null) {
                maskedArea!![nextRange.chromosome].add(nextRange.toRange())
            }

            val actualMaskedPositions = maskedArea?.let { m ->
                gq.get().sumOf { m[it].size }
            } ?: -1
            info.copy(actualRange = nextRange, actualMaskedPositions = actualMaskedPositions)
        }

        results.forEach { r ->
            assertEquals(r.expectedRange, r.actualRange, message = "Mismatch in : rVal=${r.rVal}, nPos=${r.nPos}")
            assertEquals(
                r.expectedMaskedPositions,
                r.actualMaskedPositions,
                message = "Mismatch in : rVal=${r.rVal}, nPos=${r.nPos}"
            )
            if (r.actualRange != null) {
                val actualNPos = background.fullMethylome.getCoverage(r.actualRange.on(Strand.PLUS))
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
        Tests.assertThrowsWithMessage(IndexOutOfBoundsException::class.java, "Index 2 out of bounds for length 2") {
            checkCalculatePrefixSum(COVERAGE_EX, 11, chr2, 1)
        }
    }

    @Test
    fun calculatePrefixSumWithNonCovChrs() {
        checkCalculatePrefixSum(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, 0, chr3, 0)
        checkCalculatePrefixSum(COVERAGE_WITH_NON_COVERED_FIRST_AND_LAST_CHR, 9, chr3, 0)
        Tests.assertThrowsWithMessage(IndexOutOfBoundsException::class.java, "Index 1 out of bounds for length 1") {
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
            Location(90, 105, RegionShuffleStatsTest.chr1),
            // intersects masked area
            Location(27, 60, RegionShuffleStatsTest.chr1),
            // masked:
            Location(62, 70, RegionShuffleStatsTest.chr1),
            // outside allowed:
            Location(2, 5, RegionShuffleStatsTest.chr1),
            Location(110, 200, RegionShuffleStatsTest.chr1),
        )

        val genomeAllowedAreaFilter = withTempBedFile(lociGenomeAllowed) { genomeAllowedLociPath ->
            RegionShuffleStats.readGenomeAreaFilter(genomeAllowedLociPath, gq)
        }
        
        val genomeMaskedAreaFilter = withTempBedFile(lociGenomeMasked) { genomeMaskedLociPath ->
            RegionShuffleStats.readGenomeAreaFilter(genomeMaskedLociPath, gq)
        }
        
        val inputRegionsFiltered = withTempBedFile(lociInputRegions) { inputRegionsPath ->
            // Load Input regions, methylome + filter allowed regions / methylome
            val inputRegionsFiltered = processInputRegions(
                inputRegionsPath, gq,
                genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                genomeMaskedAreaFilter = genomeMaskedAreaFilter
            )
            require(inputRegionsFiltered.isNotEmpty()) {
                "Regions file is empty or all regions were removed by filters."
            }
            inputRegionsFiltered
        }
        assertEquals(
            listOf(
                Location(9,14, chr1),
                Location(20,22, chr1),
                Location(24,27, chr1),
                Location(90,105, chr1)
            ).sorted(),
            inputRegionsFiltered.toList().sorted()
        )

        val cov = withTempFile("bgcov", ".tsv") { bgRegionsPath ->
            createBasePairCoverage(COVERAGE_EX2, gq).saveToTSV(bgRegionsPath)

            RegionShuffleStatsFromMethylomeCoverage.backgroundProviderFun(
                bgRegionsPath,
                zeroBasedBg = false,
                gq,
                genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                genomeMaskedAreaFilter = genomeMaskedAreaFilter
            )
        }

        assertEquals(12, cov.fullMethylome.depth)
        assertEquals("{1, 5, 9, 10, 20, 22, 24, 25, 26, 28, 99}", cov.fullMethylome.data[chr1].toString())
        
        assertEquals(7, cov.anchorMethylome.depth)
        assertEquals("{10, 20, 22, 24, 25, 26, 99}", cov.anchorMethylome.data[chr1].toString())
    }

    @Test
    fun loadInputRegionsAndMethylomeCovBackgroundNotAllRegionsIntersectBg() {
        val inputRegionsFiltered = listOf(
            Location(24,27, chr1),
            Location(90,105, chr1)
        ).sorted()

        val lociGenomeAllowed = listOf(
            Location(10, 100, RegionShuffleStatsTest.chr1),
        )
        val genomeAllowedAreaFilter = withTempBedFile(lociGenomeAllowed) { genomeAllowedLociPath ->
            RegionShuffleStats.readGenomeAreaFilter(genomeAllowedLociPath, gq)
        }

        withTempFile("bgcov", ".tsv") { bgRegionsPath ->
            createBasePairCoverage(COVERAGE_EX, gq).saveToTSV(bgRegionsPath)

            val bg = RegionShuffleStatsFromMethylomeCoverage.backgroundProviderFun(
                bgRegionsPath,
                zeroBasedBg = false,
                gq,
                genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                genomeMaskedAreaFilter = null
            )

            Tests.assertThrowsWithMessage(
                java.lang.IllegalArgumentException::class.java,
                "overage should cover all input regions, but the region is missing in background: chr1:[90, 105)",
                partialMessageMatch = true
            ) {
                RegionShuffleStatsFromMethylomeCoverage.ensureInputRegionsMatchesBackgound(
                    inputRegionsFiltered, bg, bgRegionsPath
                )
            }
        }
    }

    @Test
    fun testBuildIndexMappingTo() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(0, 2, chr1),
                Location(3, 21, chr1, Strand.MINUS),
                Location(24, 30, chr1),
            )
        )
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 9)
            .process(chr1, 0)
            .process(chr1, 1)
            .process(chr1, 2)
            .process(chr1, 8)
            .process(chr1, 5)
            .process(chr1, 24)
            .process(chr1, 25)
            .process(chr1, 29)
            .process(chr1, 30)
            .process(chr1, 31)
            .build(unique = false)

        val filteredCov = cov.filter(filter, includeRegions = false)

        val filteredDepth = filteredCov.depth.toInt()
        require(filteredDepth == 6) { "Actual: $filteredDepth"}

        val bg = MethylomeSamplingBackground(cov, filteredCov)

        for (i in 0 until filteredDepth) {
            val j = bg.mapAnchorOffsetToFull(chr1, i)
            assertEquals(bg.anchorMethylome.data[chr1][i], bg.fullMethylome.data[chr1][j], "Offset: $i -> $j")
        }

        assertEquals(2, bg.mapAnchorOffsetToFull(chr1, 0))
        assertEquals(3, bg.mapAnchorOffsetToFull(chr1, 1))
        assertEquals(4, bg.mapAnchorOffsetToFull(chr1, 2))
        assertEquals(5, bg.mapAnchorOffsetToFull(chr1, 3))
        assertEquals(9, bg.mapAnchorOffsetToFull(chr1, 4))
        assertEquals(10, bg.mapAnchorOffsetToFull(chr1, 5))

        assertEquals(1, bg.achorPrefixChrs.size)
        assertEquals(listOf(chr1), bg.achorPrefixChrs)
        assertEquals(6, bg.achorPrefixSum.last())
    }

    @Test
    fun testBuildIndexMappingToIdentical() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 9)
            .process(chr1, 0)
            .process(chr1, 5)
            .build(unique = false)

        val bg = MethylomeSamplingBackground(cov, cov)
        assertEquals(0, bg.mapAnchorOffsetToFull(chr1, 0))
        assertEquals(1, bg.mapAnchorOffsetToFull(chr1, 1))
        assertEquals(2, bg.mapAnchorOffsetToFull(chr1, 2))

        assertEquals(1, bg.achorPrefixChrs.size)
        assertEquals(listOf(chr1), bg.achorPrefixChrs)
        assertEquals(3, bg.achorPrefixSum.last())
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


    private fun checkTryShuffledRegions(
        coverageContent: List<Pair<Chromosome, Int>>,
        inputRegions: List<ChromosomeRange>,
        expectedResult: String?,
        withReplacement: Boolean,
        endPositionShift: Int,
        maskedGenome: List<Location>? = null,
        candidateFilterPredicate: ((ChromosomeRange) -> Double)? = null,
        singleRegionMaxRetries: Int = 1,
    ) {
        // Trick: used one storage for masked & seen intervals that is forbidden to intersect
        val maskedGenomeFilter = if (maskedGenome == null) null else LocationsMergingList.Companion.create(gq, maskedGenome)

        val fullMethylome = createBasePairCoverage(coverageContent, gq)
        val background = MethylomeSamplingBackground(fullMethylome, fullMethylome)

        val fakeRandom = Random(10)
        val (actualRegions, actualAttemptsHistogram) = RegionShuffleStatsFromMethylomeCoverage.tryShuffle(
            gq,
            background,
            background.calcExpCoverage(inputRegions),
            singleRegionMaxRetries = singleRegionMaxRetries,
            endPositionShift,
            withReplacement,
            genomeMaskedAreaFilter = maskedGenomeFilter,
            candidateFilterPredicate = candidateFilterPredicate,
            randomGen = fakeRandom
        )

        if (expectedResult != null) {
            assertNotNull(actualRegions, message = "Result: $actualRegions, ${actualAttemptsHistogram.data}")
            assertEquals(expectedResult, "$actualRegions, ${actualAttemptsHistogram.data}")
            assertEquals(inputRegions.size, actualRegions.size)
        } else {
            kotlin.test.assertNull(actualRegions, message = "Result: $actualRegions, ${actualAttemptsHistogram.data}")
        }
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