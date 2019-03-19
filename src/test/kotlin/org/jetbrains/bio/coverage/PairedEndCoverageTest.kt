package org.jetbrains.bio.coverage

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertNotIn
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.fail

class PairedEndCoverageTest {

    @Test
    fun testEmptyTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery).build(false)
        assertEquals(
            0,
            coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS))
        )
        assertEquals(
            0,
            coverage.averageInsertSize
        )
    }

    @Test
    fun testSimpleTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 90, 0, 10
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 89, 9, 11
                ).build(unique = false)
        assertEquals(
            3,
            coverage.getCoverage(Location(45, 55, chromosome1, Strand.PLUS))
        )
        assertEquals(
            0,
            coverage.getCoverage(Location(45, 55, chromosome1, Strand.MINUS))
        )
        assertEquals(
            95,
            coverage.averageInsertSize
        )
    }

    @Test
    fun testStartIndex() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 90, 0, 10
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 89, 9, 11
                ).build(unique = false)
        assertEquals(
            2,
            coverage.getCoverage(Location(51, 55, chromosome1, Strand.PLUS))
        )
    }

    @Test
    fun testSortingTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 89, 9, 11
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 90, 0, 10
                ).build(unique = false)
        assertEquals(
            3,
            coverage.getCoverage(Location(45, 55, chromosome1, Strand.PLUS))
        )
        assertEquals(
            95,
            coverage.averageInsertSize
        )
    }

    @Test
    fun testMultiplePointsTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 90, 0, 10
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 89, 9, 11
                ).process(
                    chromosome1, 90, 0, 10
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 89, 9, 11
                ).process(
                    chromosome1, 90, 0, 10
                ).process(
                    chromosome1, 80, 5, 20
                ).process(
                    chromosome1, 89, 9, 11
                ).build(unique = false)
        assertEquals(
            9,
            coverage.getCoverage(Location(45, 55, chromosome1, Strand.PLUS))
        )
        assertEquals(
            95,
            coverage.averageInsertSize
        )
    }

    @Test
    fun testWrongOutOfBoundsLeftTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 90, 0, 10
                ).build(unique = false)
        assertEquals(
            0,
            coverage.getCoverage(Location(55, 65, chromosome1, Strand.PLUS))
        )
    }

    @Test
    fun testWrongOutOfBoundsRightTagsCoverage() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .process(
                    chromosome1, 90, 0, 10
                ).build(unique = false)
        assertEquals(
            0,
            coverage.getCoverage(Location(35, 45, chromosome1, Strand.PLUS))
        )
    }

    @Test
    fun testGetTagsSimple() {
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, 5, 13, 23, 1, 111, 7, 4, 5, 50)
                .build(unique = false)

        Assert.assertArrayEquals(
            intArrayOf(5, 5, 7, 13, 23),
            coverage.getTags(ChromosomeRange(5, 50, chromosome1))
        )
        Assert.assertArrayEquals(
            intArrayOf(13, 23),
            coverage.getTags(ChromosomeRange(10, 25, chromosome1))
        )
        Assert.assertArrayEquals(
            intArrayOf(5, 5, 7, 13, 23, 50),
            coverage.getTags(ChromosomeRange(5, 55, chromosome1))
        )
    }

    @Test
    fun testGetEqualTags() {
        val tags = intArrayOf(5, 5, 5, 5, 5, 5, 5, 5, 5)
        val coverage = PairedEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, *tags)
                .build(unique = false)

        Assert.assertArrayEquals(
            tags,
            coverage.getTags(ChromosomeRange(5, 50, chromosome1))
        )
    }

    @Test
    fun testSerialization() {
        val coverage = generateCoverage()
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.save(coveragePath)
            val loaded = Coverage.load(coveragePath, genomeQuery)
            assertEquals(PairedEndCoverage::class.java, loaded::class.java)
            assertEquals(coverage.data, (loaded as PairedEndCoverage).data)
            assertEquals(coverage.averageInsertSize, loaded.averageInsertSize)
        }
    }

    @Test
    fun testPartialLoading() {
        val coverage = generateCoverage()
        assertIn(chromosome2, coverage.genomeQuery.get())
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.save(coveragePath)
            val loaded = Coverage.load(
                coveragePath,
                GenomeQuery(Genome["to1"], chromosome1.name)
            )
            assertEquals(PairedEndCoverage::class.java, loaded::class.java)
            assertNotIn(chromosome2, loaded.genomeQuery.get())
            assertEquals(coverage.averageInsertSize, (loaded as PairedEndCoverage).averageInsertSize)
        }
    }

    @Test
    fun testReversePartialLoading() {
        val partialGenomeQuery = GenomeQuery(Genome["to1"], chromosome1.name)
        val coverage = generateCoverage(partialGenomeQuery)
        assertNotIn(chromosome2, coverage.genomeQuery.get())
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.save(coveragePath)
            try {
                Coverage.load(coveragePath, genomeQuery)
                fail("Loading partial coverage with full genome completed successfully.")
            } catch (e: IllegalStateException) {
                val message = e.message
                assertNotNull(message)
                assertIn(
                    "Cache file $coveragePath doesn't contain ${chromosome2.name}",
                    message
                )
            }
        }
    }

    // TODO[dievsky] test loading from BAM files?

    companion object {
        internal var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
        internal var chromosome1: Chromosome = genomeQuery.get()[0]
        internal var chromosome2: Chromosome = genomeQuery.get()[1]

        private fun generateCoverage(gq: GenomeQuery = genomeQuery): PairedEndCoverage {
            val builder = PairedEndCoverage.builder(gq)
            gq.get().forEachIndexed { j, chromosome ->
                for (i in 0..99) {
                    // generate some valid values
                    val leftmost = 8 * i + 1
                    val rightmost = 10 * i + 3
                    val len = i + 1
                    val pos = rightmost - len
                    val pnext = leftmost
                    builder.process(chromosome, pos, pnext, len)
                }
            }

            return builder.build(false)
        }
    }
}

internal fun PairedEndCoverage.Builder.putAll(
        chromosome: Chromosome,
        vararg offsets: Int
): PairedEndCoverage.Builder {
    data[chromosome].addAll(offsets)
    return this
}