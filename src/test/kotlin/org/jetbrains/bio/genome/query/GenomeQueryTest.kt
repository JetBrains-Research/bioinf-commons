package org.jetbrains.bio.genome.query

import org.apache.log4j.SimpleLayout
import org.apache.log4j.WriterAppender
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.*
import org.junit.Test
import java.io.ByteArrayOutputStream
import kotlin.test.*

class GenomeQueryTest {

    @Test
    fun testNames() {
        val genomeQuery = GenomeQuery(Genome["to1"], "chr1", "chr2")
        assertNotNull(genomeQuery.restriction)
        assertTrue("chr1" in genomeQuery.restriction!!)
        assertTrue("chr2" in genomeQuery.restriction!!)
        assertFalse("chr3" in genomeQuery.restriction!!)
    }

    @Test
    fun testOnly() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        val same = genomeQuery.only(genomeQuery.get().map { it.name })
        assertSame(genomeQuery, same)
        assertEquals("to1", same.id)
        assertEquals("to1[chr1]", genomeQuery.only(listOf("chr1")).id)
    }

    @Test
    fun testChromSizes() {
        withTempDirectory("foo") { dir ->
            val chromSizesPath = dir / "to1.chrom.sizes"
            chromSizesPath.bufferedWriter().use {
                it.write("chr1\t100000\n")
                it.write("chr10\t100000\n")
                it.write("chr100\t100000\n")
            }
            val g0 = Genome["to1"]
            assertEquals("[chr1, chr2, chr3, chrX, chrM]", g0.chromosomes.map { it.name }.toString())
            val g1 = Genome["to1", chromSizesPath]
            assertEquals("[chr1, chr10, chr100]", g1.chromosomes.map { it.name }.toString())
            assertNotEquals(g0, g1)
        }
    }


    @Test
    fun testOnlySameGenome() {
        withTempDirectory("foo") { dir ->
            val chromSizesPath = dir / "foo.chrom.sizes"
            chromSizesPath.bufferedWriter().use {
                it.write("fooBarBaz\t10\n")
                it.write("chr1\t100000\n")
                it.write("chr10\t100000\n")
                it.write("chr100\t100000\n")
            }
            val genome = Genome["to1", chromSizesPath]
            val genomeQuery = GenomeQuery(genome)
            val restricted = genomeQuery.only(listOf("chr1"))
            assertSame(genome, genomeQuery.genome)
            assertSame(genome, restricted.genome)
            assertNotEquals(genome, GenomeQuery(Genome["to1"]).genome)
        }
    }


    @Test
    fun testNamesAllChromosomes() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        assertNull(genomeQuery.restriction)
    }

    @Test
    fun testGetShortName() {
        val genomeQuery = GenomeQuery(Genome["to1"], "chr1", "chr2")
        assertEquals("to1[chr1,chr2]", genomeQuery.id)
    }

    @Test
    fun testChromosomeNamesMap() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        for (chromosome in genomeQuery.get()) {
            assertTrue(chromosome.name.substringAfter("chr") in genomeQuery)
            assertNotNull(genomeQuery[chromosome.name.substringAfter("chr")])
            assertTrue(chromosome.name in genomeQuery)
            assertNotNull(genomeQuery[chromosome.name])
            assertTrue(chromosome.name.toLowerCase() in genomeQuery)
            assertNotNull(genomeQuery[chromosome.name.toLowerCase()])
        }
    }


    @Test
    fun restrictions() {
        val genome = Genome["to1"]
        val genomeQuery = GenomeQuery(Genome["to1"])
        val chrM = Chromosome(genome, "chrM")
        assertTrue(chrM in genome.chromosomes)
        assertNotNull(genomeQuery[chrM.name])

        val genomeQuery13 = GenomeQuery(genome, "chr1", "chr3")
        assertFalse("chr2" in genomeQuery13)
        assertFalse("CHR2" in genomeQuery13)
        assertFalse("2" in genomeQuery13)
    }

    @Test(expected = IllegalStateException::class)
    fun unknownChromosomeRestriction() {
        GenomeQuery(Genome["to1"], "unknown")
    }


    @Test
    fun testCreateGenomeQuery() {
        withTempFile("foo", ".galaxy.dat") { path ->
            val logContent = ByteArrayOutputStream()
            val appender = WriterAppender(SimpleLayout(), logContent).apply { name = "test appender" }
            GenomeQuery.LOG.addAppender(appender)
            try {
                val genomeQuery = GenomeQuery(path)
                val build = path.fileName.toString().substringBefore(".")
                assertEquals(genomeQuery.build, build)
                assertTrue("Unexpected chrom sizes file name: ${path.fileName}, " +
                        "expected <build>.chrom.sizes. Detected build: $build" in
                        logContent.toString())
            } finally {
                MultitaskProgress.LOG.removeAppender(appender)
            }
        }
    }
}
