package org.jetbrains.bio.genome.coverage

import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.util.withTempFile
import org.jetbrains.bio.util.write
import org.junit.Test
import kotlin.test.assertEquals

class BasePairCoverageTest {
    private val gq: GenomeQuery = GenomeQuery(Genome["to1"])
    private val chr1: Chromosome = gq.get()[0]
    private val chr2: Chromosome = gq.get()[1]

    @Test
    fun testLoadFromFileZeroBased() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_TWO_COLS)

            BasePairCoverage.loadFromTSV(gq, trackPath, false)
        }

        assertEquals(4, cov.depth)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(2, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testLoadFromFileOneBased() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_TWO_COLS)

            BasePairCoverage.loadFromTSV(gq, trackPath, true)
        }

        assertEquals(4, cov.depth)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testLoadFromFileZeroBasedManyCols() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_MANY_COLS)

            BasePairCoverage.loadFromTSV(gq, trackPath, false)
        }

        assertEquals(4, cov.depth)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testLoadFromFileZeroGrchChrs() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_TWO_COLS_GRCH)

            BasePairCoverage.loadFromTSV(gq, trackPath, false)
        }

        assertEquals(4, cov.depth)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testLoadFromFileZeroBasedNotUnique() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_TWO_COLS_NOT_UNIQUE)

            BasePairCoverage.loadFromTSV(gq, trackPath, false)
        }

        assertEquals(4, cov.depth)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testLoadFromFileZeroWithUnsupportedChrs() {
        withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_WRONG_CHR)

            Tests.assertThrowsWithMessage(
                IllegalArgumentException::class.java,
                "Unknown chromosome 'chrA' for genome: 'Test Organism: to1 [test]'",
            ) {
                BasePairCoverage.loadFromTSV(gq, trackPath, false)
            }
        }
    }

    @Test
    fun testLoadFromFileZeroWithUnsupportedChrsIgnored() {
        val cov = withTempFile("cov", ".txt") { trackPath ->
            trackPath.write(COVERAGE_WRONG_CHR)

            BasePairCoverage.loadFromTSV(gq, trackPath, false, failOnMissingChromosomes = false)
        }
        assertEquals(2, cov.depth)
        assertEquals(1, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(1, cov.getCoverage(Location(4, 6, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))

    }

    @Test
    fun testEmptyCoverage() {
        val cov = BasePairCoverage.builder(gq, false).build(false)
        assertEquals(0, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
        assertEquals(0, cov.getCoverage(Location(0, 10, chr1, Strand.MINUS)))
        assertEquals(0, cov.depth)
    }

    @Test
    fun testCoverageWithDuplicates() {
        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 2)
            .process(chr1, 3)
            .process(chr1, 3)
            .process(chr1, 3)
            .process(chr1, 5)
            .process(chr2, 3)
            .process(chr2, 3)
            .process(chr2, 5)
            .build(false)
        assertEquals(9, cov.depth)
    }

    @Test
    fun testUniqueCoverageWithDuplicates() {
        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 2)
            .process(chr1, 3)
            .process(chr1, 3)
            .process(chr1, 3)
            .process(chr1, 5)
            .process(chr2, 3)
            .process(chr2, 3)
            .process(chr2, 5)
            .build(true)
        assertEquals(5, cov.depth)
    }

    @Test
    fun testBuilderRangeCheck() {
        Tests.assertThrowsWithMessage(
            IllegalArgumentException::class.java,
            "Zero-based offset should be >= 1, but was: -1",
        ) {
            BasePairCoverage.builder(gq, offsetIsOneBased = false).process(chr1, -1)
        }
        Tests.assertThrowsWithMessage(
            IllegalArgumentException::class.java,
            "One-based offset should be < 10000000 (chr1 size), but was: 10000000",
        ) {
            BasePairCoverage.builder(gq, offsetIsOneBased = false).process(chr1, chr1.length)
        }

        Tests.assertThrowsWithMessage(
            IllegalArgumentException::class.java,
            "One-based offset should be >= 1, but was: 0",
        ) {
            BasePairCoverage.builder(gq, offsetIsOneBased = true).process(chr1, 0)
        }
        Tests.assertThrowsWithMessage(
            IllegalArgumentException::class.java,
            "One-based offset should be <= 10000000 (chr1 size), but was: 10000001",
        ) {
            BasePairCoverage.builder(gq, offsetIsOneBased = true).process(chr1, chr1.length + 1)
        }

        assertEquals(
            1,
            BasePairCoverage.builder(gq, offsetIsOneBased = false).process(chr1, 0).build(false).depth
        )
        assertEquals(
            1,
            BasePairCoverage.builder(gq, offsetIsOneBased = false).process(chr1, chr1.length - 1)
                .build(false).depth
        )
        assertEquals(
            1,
            BasePairCoverage.builder(gq, offsetIsOneBased = true).process(chr1, chr1.length)
                .build(false).depth
        )
    }

    @Test
    fun testSimpleTagsCoverage() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testSimpleTagsCoverageWrongChr() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)
        assertEquals(0, cov.getCoverage(Location(0, 10, chr2, Strand.PLUS)))
    }

    @Test
    fun testSimpleTagsCoverageMultipleChr() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .process(chr2, 5)
            .build(unique = false)
        assertEquals(1, cov.getCoverage(Location(0, 10, chr2, Strand.PLUS)))
    }

    @Test
    fun testStartIndex() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(1, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testEndIndex() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(0, 9, chr1, Strand.PLUS)))
    }

    @Test
    fun testStartIndexOneBased() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 1)
            .process(chr1, 6)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(1, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testEndIndexOneBased() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 1)
            .process(chr1, 6)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(0, 9, chr1, Strand.PLUS)))
    }

    @Test
    fun testStartIndexNotFirstOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(5, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testStartIndexOneBasedNotFirstOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 1)
            .process(chr1, 6)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(2, cov.getCoverage(Location(5, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testEndIndexNotLastOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(0, 5, chr1, Strand.PLUS)))
    }

    @Test
    fun testEndIndexOneBasedNotLastOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 1)
            .process(chr1, 6)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(0, 5, chr1, Strand.PLUS)))
    }

    @Test
    fun testStartIndexOutLastOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(9, 10, chr1, Strand.PLUS)))
    }

    @Test
    fun testStartIndexOneBasedLastOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 1)
            .process(chr1, 6)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(9, 10, chr1, Strand.PLUS)))
    }


    @Test
    fun testEndIndexFirstOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 2)
            .process(chr1, 5)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testEndIndexOneBasedFirstOffset() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 3)
            .process(chr1, 6)
            .build(unique = false)

        assertEquals(1, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testEmptyRangeOutOfCoverageLeft() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 3)
            .build(unique = false)

        assertEquals(0, cov.getCoverage(Location(1, 3, chr1, Strand.PLUS)))
    }

    @Test
    fun testEmptyRangeOutOfCoverageMiddle() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 0)
            .process(chr1, 5)
            .build(unique = false)

        assertEquals(0, cov.getCoverage(Location(1, 5, chr1, Strand.PLUS)))
    }

    @Test
    fun testEmptyRangeOneBasedOutOfCoverageMiddle() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 1)
            .process(chr1, 6)
            .build(unique = false)

        assertEquals(0, cov.getCoverage(Location(1, 5, chr1, Strand.PLUS)))
    }

    @Test
    fun testEmptyRangeOutOfCoverageRight() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 9)
            .build(unique = false)

        assertEquals(0, cov.getCoverage(Location(10, 11, chr1, Strand.PLUS)))
    }

    @Test
    fun testEmptyRangeOneBasedOffsetOutOfCoverageRight() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = true)
            .process(chr1, 10)
            .build(unique = false)

        assertEquals(0, cov.getCoverage(Location(10, 11, chr1, Strand.PLUS)))
    }

    @Test
    fun testSortingTagsCoverage() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 9)
            .process(chr1, 0)
            .process(chr1, 5)
            .build(unique = false)
        assertEquals(3, cov.getCoverage(Location(0, 10, chr1, Strand.PLUS)))
    }

    companion object {
        val COVERAGE_TWO_COLS = """
            chr1	9
            chr1	1
            chr2	3
            chr1	5
        """.trimIndent().trim()

        val COVERAGE_TWO_COLS_NOT_UNIQUE = """
            chr1	9
            chr1	0
            chr1	0
            chr2	3
            chr2	3
            chr2	3
            chr1	5
        """.trimIndent().trim()

        val COVERAGE_TWO_COLS_GRCH = """
            1	9
            1	0
            2	3
            1	5
        """.trimIndent().trim()

        val COVERAGE_MANY_COLS = """
            chr1	9	-	0.10
            chr1	0	+	0.2
            chr2	3	+	0.77
            chr1	4	+	1
        """.trimIndent().trim()

        val COVERAGE_WRONG_CHR = """
            chrA	9	-	0.10
            chrB	0	+	0.2
            chr2	3	+	0.77
            chr1	4	+	1
        """.trimIndent().trim()

    }

    @Test
    fun testSerializationToNpz() {
        val builder = BasePairCoverage.builder(gq, offsetIsOneBased = false)
        gq.get().forEachIndexed { j, chromosome ->
            for (i in 0..99) {
                builder.process(chromosome, i + 50)
            }
        }

        val coverage = builder.build(false)
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.saveToNpz(coveragePath)
            val loaded = Coverage.load(coveragePath, gq)
            assertEquals(BasePairCoverage::class.java, loaded::class.java)
            assertEquals(coverage.depth, (loaded as BasePairCoverage).depth)
            assertEquals(coverage.data, loaded.data)
        }
    }

    @Test
    fun testSerializationToTsv() {
        val builder = BasePairCoverage.builder(gq, offsetIsOneBased = false)
        gq.get().forEachIndexed { j, chromosome ->
            for (i in 0..99) {
                builder.process(chromosome, i + 50)
            }
        }

        for (offsetIsOneBased in listOf(false, true)) {
            val coverage = builder.build(false)
            withTempFile("coverage", ".tsv") { coveragePath ->
                coverage.saveToTSV(coveragePath, offsetIsOneBased = offsetIsOneBased)
                val loaded = BasePairCoverage.loadFromTSV(gq, coveragePath, offsetIsOneBased = offsetIsOneBased)
                assertEquals(BasePairCoverage::class.java, loaded::class.java)
                assertEquals(coverage.depth, loaded.depth)
                assertEquals(coverage.data, loaded.data)
            }
        }
    }

    @Test
    fun testFilterExclude() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(0, 2, chr1),
                Location(3, 21, chr1, Strand.MINUS),
                Location(24, 30, chr1),

                Location(11, 15, chr2),
                Location(16, 50, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 5)
            .process(chr1, 10)
            .process(chr1, 20)
            .process(chr1, 23)
            .process(chr1, 30)
            .process(chr2, 3)
            .process(chr2, 10)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = false)

        assertEquals(9, cov.depth)
        assertEquals("{2, 5, 10, 20, 23, 30}", cov.data[chr1].toString())
        assertEquals("{3, 10, 15}", cov.data[chr2].toString())
    }

    @Test
    fun testFilterExcludeStranded() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(0, 2, chr1),
                Location(3, 21, chr1, Strand.MINUS),
                Location(24, 30, chr1),

                Location(11, 15, chr2),
                Location(16, 50, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 5)
            .process(chr1, 10)
            .process(chr1, 20)
            .process(chr1, 23)
            .process(chr1, 30)
            .process(chr2, 3)
            .process(chr2, 10)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = false, ignoreRegionsOnMinusStrand = false)

        assertEquals(6, cov.depth)
        assertEquals("{2, 23, 30}", cov.data[chr1].toString())
        assertEquals("{3, 10, 15}", cov.data[chr2].toString())
    }


    @Test
    fun testFilterExcludeAll() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(1, 20, chr1),
                Location(0, 20, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 5)
            .process(chr1, 10)
            .process(chr2, 3)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = false)

        assertEquals(0, cov.depth)
    }


    @Test
    fun testFilterInclude() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(0, 2, chr1),
                Location(3, 21, chr1, Strand.MINUS),
                Location(24, 30, chr1),

                Location(11, 15, chr2),
                Location(16, 50, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 1)
            .process(chr1, 2)
            .process(chr1, 10)
            .process(chr1, 20)
            .process(chr1, 23)
            .process(chr1, 30)
            .process(chr1, 32)
            .process(chr2, 3)
            .process(chr2, 10)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = true)

        assertEquals(1, cov.depth)
        assertEquals("{1}", cov.data[chr1].toString())
        assertEquals("{}", cov.data[chr2].toString())
    }

    @Test
    fun testFilterIncludeStranded() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(0, 2, chr1),
                Location(3, 21, chr1, Strand.MINUS),
                Location(24, 30, chr1),

                Location(11, 15, chr2),
                Location(16, 50, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 1)
            .process(chr1, 2)
            .process(chr1, 10)
            .process(chr1, 20)
            .process(chr1, 23)
            .process(chr1, 30)
            .process(chr1, 32)
            .process(chr2, 3)
            .process(chr2, 10)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = true, ignoreRegionsOnMinusStrand = false)

        assertEquals(3, cov.depth)
        assertEquals("{1, 10, 20}", cov.data[chr1].toString())
        assertEquals("{}", cov.data[chr2].toString())
    }


    @Test
    fun testFilterIncludeAll() {
        val filter = LocationsMergingList.create(
            gq,
            listOf(
                Location(1, 20, chr1),
                Location(3, 16, chr2),
            )
        )

        val cov = BasePairCoverage.builder(gq, false)
            .process(chr1, 2)
            .process(chr1, 5)
            .process(chr1, 10)
            .process(chr2, 3)
            .process(chr2, 15)
            .build(false)
            .filter(filter, includeRegions = true)

        assertEquals(5, cov.depth)
        assertEquals("{2, 5, 10}", cov.data[chr1].toString())
        assertEquals("{3, 15}", cov.data[chr2].toString())
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

        Tests.assertThrowsWithMessage(
            java.lang.IllegalArgumentException::class.java,
            message = "Subset methylome should be part of full, but subset offset[idx=0]:0 not found in original methylome and doesn't match to offset[idx=0]:2",
        ) {
            cov.buildIndexMappingTo(filteredCov)
        }

        val mapping = filteredCov.buildIndexMappingTo(cov)
        for (i in 0 until filteredDepth) {
            val j = mapping[chr1]!![i]
            assertEquals(filteredCov.data[chr1][i], cov.data[chr1][j], "Offset: $i -> $j")
        }

        assertEquals(gq.get().size, mapping.size)
        assertEquals(0, mapping[chr2]!!.size())

        assertEquals(6, mapping[chr1]!!.size())
        assertEquals(2, mapping[chr1]!![0])
        assertEquals(3, mapping[chr1]!![1])
        assertEquals(4, mapping[chr1]!![2])
        assertEquals(5, mapping[chr1]!![3])
        assertEquals(9, mapping[chr1]!![4])
        assertEquals(10, mapping[chr1]!![5])
    }

    @Test
    fun testBuildIndexMappingToIdentical() {
        val cov = BasePairCoverage.builder(gq, offsetIsOneBased = false)
            .process(chr1, 9)
            .process(chr1, 0)
            .process(chr1, 5)
            .build(unique = false)

        val mapping = cov.buildIndexMappingTo(cov)
        assertEquals(gq.get().size, mapping.size)
        assertEquals(3, mapping[chr1]!!.size())
        assertEquals(0, mapping[chr2]!!.size())
        assertEquals(0, mapping[chr1]!![0])
        assertEquals(1, mapping[chr1]!![1])
        assertEquals(2, mapping[chr1]!![2])
    }
}