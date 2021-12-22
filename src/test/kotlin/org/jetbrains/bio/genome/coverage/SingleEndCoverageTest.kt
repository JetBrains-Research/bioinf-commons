package org.jetbrains.bio.genome.coverage

import gnu.trove.list.array.TIntArrayList
import kotlinx.support.jdk7.use
import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertNotIn
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.coverage.FragmentSize.detectFragmentSize
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.processReads
import org.jetbrains.bio.genome.format.toBedEntry
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert.assertArrayEquals
import org.junit.Test
import kotlin.test.*

class SingleEndCoverageTest {

    @Test
    fun testEmptyTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery).build(false)
        assertEquals(0, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
        assertEquals(0, coverage.detectedFragment)
    }

    @Test
    fun testSimpleTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(3, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testStartIndex() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(2, coverage.getCoverage(Location(1, 10, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testSortingTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(3, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testMultiplePointsTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(5, 100, chromosome1, Strand.PLUS))
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .process(Location(9, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(9, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testWrongStrandTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(1, 100, chromosome1, Strand.MINUS))
            .build(unique = false).withFragment(0)
        assertEquals(0, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testWrongOutOfBoundsLeftTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(0, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(0, coverage.getCoverage(Location(10, 20, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testWrongOutOfBoundsRightTagsCoverage() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .process(Location(30, 100, chromosome1, Strand.PLUS))
            .build(unique = false).withFragment(0)
        assertEquals(0, coverage.getCoverage(Location(10, 20, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testGetTagsSimple() {
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .putAll(chromosome1, Strand.PLUS, 5, 13, 23, 1, 111, 7, 4, 5, 50)
            .build(unique = false).withFragment(0)

        assertArrayEquals(
            intArrayOf(5, 5, 7, 13, 23),
            coverage.getTags(Location(5, 50, chromosome1, Strand.PLUS))
        )
        assertArrayEquals(
            intArrayOf(13, 23),
            coverage.getTags(Location(10, 25, chromosome1, Strand.PLUS))
        )
        assertArrayEquals(
            intArrayOf(5, 5, 7, 13, 23, 50),
            coverage.getTags(Location(5, 55, chromosome1, Strand.PLUS))
        )
    }

    @Test
    fun testGetEqualTags() {
        val tags = intArrayOf(5, 5, 5, 5, 5, 5, 5, 5, 5)
        val coverage = SingleEndCoverage.builder(genomeQuery)
            .putAll(chromosome1, Strand.PLUS, *tags)
            .build(unique = false).withFragment(0)

        assertArrayEquals(tags, coverage.getTags(Location(5, 50, chromosome1, Strand.PLUS)))
    }

    @Test
    fun testSerialization() {
        val builder = SingleEndCoverage.builder(genomeQuery)
        genomeQuery.get().forEachIndexed { j, chromosome ->
            for (i in 0..99) {
                val strand = if ((i + j) % 2 == 0) Strand.PLUS else Strand.MINUS
                builder.process(Location(i, i + 50, chromosome, strand))
            }
        }

        val coverage = builder.build(false).withFragment(147)
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.save(coveragePath)
            val loaded = Coverage.load(coveragePath, genomeQuery)
            assertEquals(SingleEndCoverage::class.java, loaded::class.java)
            assertEquals(
                coverage.detectedFragment,
                (loaded as SingleEndCoverage).detectedFragment
            )
            assertEquals(coverage.data, loaded.data)
        }
    }

    @Test
    fun testPartialLoading() {
        val coverage = generateCoverage()
        withTempFile("coverage", ".cov") { coveragePath ->
            coverage.save(coveragePath)
            assertTrue(chromosome2 in coverage.genomeQuery.get())
            val loaded = Coverage.load(coveragePath, GenomeQuery(Genome["to1"], chromosome1.name))
            assertEquals(SingleEndCoverage::class.java, loaded::class.java)
            assertFalse(chromosome2 in loaded.genomeQuery.get())
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
                    "File $coveragePath doesn't contain data for ${PairedEndCoverageTest.chromosome2.name}.",
                    message
                )
            }
        }
    }

    @Test
    fun testSortInComputeTagsCoverage() {
        withTempFile("track", ".bed") { trackPath ->
            val bedFormat = BedFormat()
            bedFormat.print(trackPath).use { bedPrinter ->
                bedPrinter.print(Location(0, 3, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(2, 4, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(4, 7, chromosome1, Strand.PLUS).toBedEntry())

                // Not sorted!
                bedPrinter.print(Location(4, 7, chromosome2, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(2, 9, chromosome2, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 11, chromosome2, Strand.PLUS).toBedEntry())
            }

            val coverage = SingleEndCoverage.builder(genomeQuery).apply {
                processReads(genomeQuery, trackPath) {
                    this.process(it)
                }
            }.build(true).withFragment(0)
            assertEquals(3, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
            assertEquals(3, coverage.getCoverage(Location(0, 10, chromosome2, Strand.PLUS)))
        }
    }

    @Test
    fun testUniqueTags() {
        withTempFile("track", ".bed") { trackPath ->
            val bedFormat = BedFormat()
            bedFormat.print(trackPath).use { bedPrinter ->
                bedPrinter.print(Location(0, 1, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 1, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 1, chromosome1, Strand.PLUS).toBedEntry())
            }

            val coverage = SingleEndCoverage.builder(genomeQuery).apply {
                processReads(genomeQuery, trackPath) {
                    this.process(it)
                }
            }.build(true).withFragment(0)

            assertEquals(1, coverage.getCoverage(Location(0, 10, chromosome1, Strand.PLUS)))
        }
    }

    @Test
    fun testFragmentSize() {
        withTempFile("track", ".bed") { trackPath ->
            val bedFormat = BedFormat()
            bedFormat.print(trackPath).use { bedPrinter ->
                bedPrinter.print(Location(0, 100, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 50, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 25, chromosome1, Strand.PLUS).toBedEntry())
                bedPrinter.print(Location(0, 101, chromosome1, Strand.MINUS).toBedEntry())
                bedPrinter.print(Location(0, 151, chromosome1, Strand.MINUS).toBedEntry())
            }

            val coverage = SingleEndCoverage.builder(genomeQuery).apply {
                processReads(genomeQuery, trackPath) {
                    this.process(it)
                }
            }.build(true).withFragment(150)
            assertArrayEquals(
                intArrayOf(0),
                coverage.getTags(Location(0, 100, chromosome1, Strand.PLUS))
            )
            assertArrayEquals(
                intArrayOf(0),
                coverage.getTags(Location(0, 500, chromosome1, Strand.PLUS))
            )
            assertArrayEquals(
                intArrayOf(100, 150),
                coverage.getTags(Location(0, 100, chromosome1, Strand.MINUS))
            )
        }
    }


    companion object {
        internal var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
        internal var chromosome1: Chromosome = genomeQuery.get()[0]
        internal var chromosome2: Chromosome = genomeQuery.get()[1]

        private fun generateCoverage(gq: GenomeQuery = genomeQuery): SingleEndCoverage {
            val builder = SingleEndCoverage.builder(gq)
            gq.get().forEachIndexed { j, chromosome ->
                for (i in 0..99) {
                    val strand = if ((i + j) % 2 == 0) Strand.PLUS else Strand.MINUS
                    builder.process(Location(i, i + 50, chromosome, strand))
                }
            }

            return builder.build(false).withFragment(0)
        }
    }
}

/**
 * Simplifies reading and writing test code
 */
internal fun SingleEndCoverage.withFragment(fragment: Int): SingleEndCoverage = withFragment(FixedFragment(fragment))


internal fun SingleEndCoverage.Builder.putAll(
    chromosome: Chromosome,
    strand: Strand,
    vararg offsets: Int
): SingleEndCoverage.Builder {
    data[chromosome, strand].addAll(offsets)
    return this
}

class BinarySearchLeftTest {
    @Test
    fun empty() = assertEquals(0, TIntArrayList().binarySearchLeft(42))

    @Test
    fun singletonLower() {
        assertEquals(1, TIntArrayList(intArrayOf(0)).binarySearchLeft(42))
    }

    @Test
    fun singletonGreaterEqual() {
        assertEquals(0, TIntArrayList(intArrayOf(42)).binarySearchLeft(0))
        assertEquals(0, TIntArrayList(intArrayOf(42)).binarySearchLeft(42))
    }

    @Test
    fun duplicates() {
        assertEquals(0, TIntArrayList(intArrayOf(42, 42, 42)).binarySearchLeft(42))
        assertEquals(1, TIntArrayList(intArrayOf(0, 42, 42, 42)).binarySearchLeft(42))
    }
}

