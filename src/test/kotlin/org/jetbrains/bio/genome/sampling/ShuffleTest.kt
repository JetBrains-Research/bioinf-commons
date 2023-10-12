package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.util.formatLongNumber
import org.junit.Test
import kotlin.random.Random
import kotlin.test.*


class ShuffleTest {
    val genomeQuery = GenomeQuery(Genome["to1"])
    val chromosomes = genomeQuery.get()
    val chr1 = chromosomes[0]
    val chr2 = chromosomes[1]
    val chr3 = chromosomes[2]

    private fun checkShuffled(regions: List<ChromosomeRange>, shuffled: List<ChromosomeRange>) {
        val l1 = regions.map { it.length() }.sorted().toList()
        val l2 = shuffled.map { it.length() }.sorted().toList()
        assertEquals(l1, l2, "Shuffled regions must have same length set")
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
    fun shuffleRegionsWithBackgroud() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val background = listOf(
            ChromosomeRange(100, 300, chromosomes[0]),
            ChromosomeRange(1000, 2000, chromosomes[1]),
            ChromosomeRange(2000, 3000, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(
            genomeQuery, regions, background, withReplacement = false, genomeMaskedAreaFilter = null,
            randomGen = Random(10)
        ).first

        checkShuffled(regions, shuffled)
        for (r in shuffled) {
            var intersectsBackground = false
            for (b in background) {
                if (ChromosomeRange.intersects(r, b)) {
                    intersectsBackground = true
                }
            }
            assertTrue(intersectsBackground, "All region start must intersect backgroud.")
        }
    }

    @Test
    fun shuffleRegionsWithBackgroudWithReplacement() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val background = listOf(
            ChromosomeRange(100, 300, chromosomes[0]),
            ChromosomeRange(1000, 2000, chromosomes[1]),
            ChromosomeRange(2000, 3000, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(
            genomeQuery, regions, background, withReplacement = true, genomeMaskedAreaFilter = null,
            randomGen = Random(10)
        ).first

        val l1 = regions.map { it.length() }.sorted().toList()
        val l2 = shuffled.map { it.length() }.sorted().toList()
        assertEquals(l1, l2, "Shuffled regions must have same length set")

        for (r in shuffled) {
            var inBackground = false
            for (b in background) {
                if (r.chromosome == b.chromosome && r.startOffset in b.toRange()) {
                    inBackground = true
                }
            }
            assertTrue(inBackground, "All region start must be in backgroud.")
        }
    }

    @Test
    fun shuffleRegions() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(
            genomeQuery, regions, withReplacement = false, genomeMaskedAreaFilter = null,
            randomGen = Random(10)
        ).first
        checkShuffled(regions, shuffled)
    }

    @Test
    fun shuffleRegionsWithReplacement() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(
            genomeQuery, regions, withReplacement = true, genomeMaskedAreaFilter = null,
            randomGen = Random(10)
        ).first
        checkShuffled(regions, shuffled)
    }

    @Test
    fun testTryShuffleRegionsWithMaskedRegionWithoutReplacement() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),           
            background = listOf(
                ChromosomeRange(5, 40, chr1),
                ChromosomeRange(0, 5, chr2),
            ),
            maskedGenome = listOf(
                Location(23, 26, chr1),
                Location(1, 2, chr2),
            ),
            expectedResult = "[chr1:[31, 33), chr1:[0, 19), chr2:[2, 3)], {2=2, 1=1}",
            withReplacement = false,
            singleRegionMaxRetries = 10
        )
    }
    @Test
    fun testTryShuffleRegionsWithMaskedRegionWithReplacement() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            background = listOf(
                ChromosomeRange(5, 40, chr1),
            ),
            maskedGenome = listOf(
                Location(23, 26, chr1),
                Location(1, 2, chr2),
            ),
            expectedResult = "[chr1:[11, 13), chr1:[3, 22), chr1:[27, 28)], {2=2, 1=1}",
            withReplacement = true,
            singleRegionMaxRetries = 10
        )
    }

    @Test
    fun testTryShuffleRegionsNotPossibleNotEnoughLength() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            background = listOf(
                ChromosomeRange(5, 10, chr1),
            ),
            maskedGenome = listOf(
                Location(0, 5, chr1),
                Location(10, 100, chr1),
            ),
            expectedResult = null,
            withReplacement = true,
            singleRegionMaxRetries = 30
        )
    }

    @Test
    fun testTryShuffleRegionsWithMaskedRegionNotPossibleByMasking() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 5, chr1), // out of background
                ChromosomeRange(28, 29, chr1), // masked
            ),
            background = listOf(
                ChromosomeRange(50, 80, chr1),
            ),
            maskedGenome = listOf(
                Location(10, 100, chr1),
            ),
            expectedResult = null,
            withReplacement = true,
            singleRegionMaxRetries = 30
        )
    }

    @Test
    fun testTryShuffleRegionsEmptyNotSupported() {
        Tests.assertThrowsWithMessage(
            IllegalArgumentException::class.java, "Empty region not supported"
        ) {
            checkTryShuffledRegions(
                inputRegions = listOf(
                    ChromosomeRange(1, 1, chr1),
                ),
                background = listOf(
                    ChromosomeRange(50, 80, chr1),
                ),
                expectedResult = "",
                withReplacement = true,
                singleRegionMaxRetries = 30
            )
        }
    }

    @Test
    fun testTryShuffleRegionsInputNotChecked() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 5, chr1), // out of background
                ChromosomeRange(28, 29, chr1), // masked
            ),
            background = listOf(
                ChromosomeRange(50, 80, chr1),
            ),
            maskedGenome = listOf(
                Location(10, 60, chr1),
                Location(70, 100, chr1),
            ),
            expectedResult = "[chr1:[63, 67), chr1:[67, 68)], {3=1, 2=1}",
            withReplacement = true,
            singleRegionMaxRetries = 30
        )
    }

    @Test
    fun testTryShuffleRegionsFromSingleChr() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr3),
                ChromosomeRange(9, 24, chr3),
                ChromosomeRange(28, 29, chr3),
            ),
            background = listOf(
                ChromosomeRange(100, 300, chr3),
                ChromosomeRange(1000, 2000, chr3),
                ChromosomeRange(2000, 3000, chr3)
            ),
            expectedResult = "[chr3:[1778, 1780), chr3:[1452, 1467), chr3:[1275, 1276)], {1=3}",
            withReplacement = false,
            singleRegionMaxRetries = 2
        )
    }

    @Test
    fun testTryShuffleRegionsFromChrDataZeroPos() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 30, chr3),
            ),
            background = listOf(
                ChromosomeRange(0, 1, chr3),
            ),
            expectedResult = "[chr3:[0, 29)], {1=1}",
            withReplacement = false,
            singleRegionMaxRetries = 2
        )
    }

    @Test
    fun testTryShuffleRegionsWithoutReplacement() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            background = listOf(
                ChromosomeRange(35, 40, chr1),
                ChromosomeRange(0, 5, chr2),
            ),
            expectedResult = "[chr2:[3, 5), chr1:[20, 39), chr2:[2, 3)], {3=1, 1=2}",
            withReplacement = false,
            singleRegionMaxRetries = 10
        )
    }
    
    @Test
    fun testTryShuffleRegionsWithReplacement() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 3, chr1),
                ChromosomeRange(5, 24, chr1),
                ChromosomeRange(28, 29, chr1),
            ),
            background = listOf(
                ChromosomeRange(35, 40, chr1),
                ChromosomeRange(0, 5, chr2),
            ),
            expectedResult = "[chr2:[3, 5), chr2:[0, 19), chr2:[0, 1)], {1=3}",
            withReplacement = true,
            singleRegionMaxRetries = 10
        )
    }

    @Test
    fun testTryShuffleRegionsWithoutReplacementWithoutBackground() {
        print(genomeQuery.get().map { "$it -> ${it.length.formatLongNumber()}" })
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 500_000, chr1),
                ChromosomeRange(1_000_000, 9_000_000, chr1),
                ChromosomeRange(5_000_000, 5_800_000, chr1),
                ChromosomeRange(6_000_000, 6_800_000, chr1),
            ),
            background = null,
            expectedResult = "[chr1:[676448, 1176447), chr1:[1469293, 9469293), chr2:[0, 800000), chrX:[0, 800000)], {33=1, 4=1, 1=2}",
            withReplacement = false,
            singleRegionMaxRetries = 50
        )
    }

    @Test
    fun testTryShuffleRegionsWithReplacementWithoutBackground() {
        checkTryShuffledRegions(
            inputRegions = listOf(
                ChromosomeRange(1, 500_000, chr1),
                ChromosomeRange(1_000_000, 9_000_000, chr1),
                ChromosomeRange(5_000_000, 5_800_000, chr1),
                ChromosomeRange(6_000_000, 6_800_000, chr1),
            ),
            background = null,
            expectedResult = "[chr1:[676448, 1176447), chr1:[0, 8000000), chr1:[3261304, 4061304), chr1:[7765070, 8565070)], {2=2, 1=2}",
            withReplacement = true,
            singleRegionMaxRetries = 10
        )
    }

    private fun checkTryShuffledRegions(
        inputRegions: List<ChromosomeRange>,
        background: List<ChromosomeRange>?,
        expectedResult: String?,
        withReplacement: Boolean,
        maskedGenome: List<Location>? = null,
        singleRegionMaxRetries: Int = 1,
    ) {
        // Trick: used one storage for masked & seen intervals that is forbidden to intersect
        val maskedGenomeFilter = if (maskedGenome == null) null else LocationsMergingList.Companion.create(
            genomeQuery, maskedGenome)

        val backgroundRegions = background ?: genomeQuery.get().map {
            ChromosomeRange(0, it.length, it)
        }
        
        val lengths = inputRegions.map { it.length() }.toIntArray()
        val prefixSum = createPrefixSumFor(backgroundRegions)

        val fakeRandom = Random(10)
        val (actualRegions, actualAttemptsHistogram) = tryShuffle(
            genomeQuery,
            backgroundRegions,
            lengths,
            prefixSum,
            singleRegionMaxRetries = singleRegionMaxRetries,
            withReplacement,
            genomeMaskedAreaFilter = maskedGenomeFilter,
            randomGen = fakeRandom
        )
        if (expectedResult != null) {
            assertNotNull(actualRegions, message = "Result: $actualRegions, ${actualAttemptsHistogram.data}")
            assertEquals(expectedResult, "$actualRegions, ${actualAttemptsHistogram.data}")
            assertEquals(inputRegions.size, actualRegions.size)
        } else {
            assertNull(actualRegions, message = "Result: $actualRegions, ${actualAttemptsHistogram.data}")
        }
    }
}