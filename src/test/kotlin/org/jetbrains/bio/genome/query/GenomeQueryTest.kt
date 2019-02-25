package org.jetbrains.bio.genome.query

import org.apache.log4j.SimpleLayout
import org.apache.log4j.WriterAppender
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.util.*
import org.junit.Test
import java.io.ByteArrayOutputStream
import kotlin.test.*

class GenomeQueryTest {
    @Test
    fun testUnspecifiedChromosome() {
        val genomeQuery = GenomeQuery("to1")
        val chromosomes = genomeQuery.get()
        assertEquals(Genome["to1"].chromosomes.size - 1,
                chromosomes.size)  // ^^^ minus chrM.
    }

    @Test
    fun testNames() {
        val genomeQuery = GenomeQuery("to1", "chr1", "chr2")
        assertTrue("chr1" in genomeQuery.restriction)
        assertTrue("chr2" in genomeQuery.restriction)
        assertFalse("chr3" in genomeQuery.restriction)
    }

    @Test
    fun testOnly() {
        val genomeQuery = GenomeQuery("to1")
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
                it.write("fooBarBaz\t10\n")
                it.write("chr1\t100000\n")
                it.write("chr10\t100000\n")
                it.write("chr100\t100000\n")
            }
            val g0 = Genome["to1"]
            assertEquals("[chr1, chr2, chr3, chrX, chrM]", g0.chromosomes.map { it.name }.toString())
            val g1 = Genome["to1", chromSizesPath]
            assertEquals("[fooBarBaz, chr1, chr10, chr100]", g1.chromosomes.map { it.name }.toString())
            assertEquals("[chr1, chr10, chr100]", g1.toQuery().get().map { it.name }.toString())
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
            assertNotEquals(genome, GenomeQuery("to1").genome)
        }
    }


    @Test
    fun testNamesAllChromosomes() {
        val genomeQuery = GenomeQuery("to1")
        assertTrue(genomeQuery.restriction.isEmpty())
    }

    @Test
    fun testSize() {
        val genomeQuery = GenomeQuery("to1", "chr1", "chr2")
        assertEquals(2, genomeQuery.get().size)

        val genomeQueryAll = GenomeQuery("to1")
        assertEquals(Genome["to1"].chromosomes.size - 1,
                genomeQueryAll.get().size)
    }

    @Test
    fun testGetShortName() {
        val genomeQuery = GenomeQuery("to1", "chr1", "chr2")
        assertEquals("to1[chr1,chr2]", genomeQuery.id)
    }

    @Test
    fun testParse() {
        assertEquals(Genome["to1"].chromosomes.size - 1,
                "to1".toGenomeQuery().get().size)
        assertEquals(1, "to1[chr1]".toGenomeQuery().get().size)
        assertEquals(2, "to1[chr1,chr3]".toGenomeQuery().get().size)
    }

    @Test
    fun testChromosomeNamesMap() {
        val genomeQuery = GenomeQuery("to1")
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
    fun chrDefaultChoice() {
        val genome = Genome["to1"]
        val genomeQuery = GenomeQuery("to1")
        val chrM = Chromosome(genome, "chrM")
        assertTrue(chrM in genome.chromosomes)
        assertEquals(null, genomeQuery[chrM.name])
        assertFalse(chrM in genomeQuery.get())
        assertFalse(chrM.name in genomeQuery)

        val genomeQuery13 = "to1[chr1,chr3]".toGenomeQuery()
        assertFalse("chr2" in genomeQuery13)
        assertFalse("CHR2" in genomeQuery13)
        assertFalse("2" in genomeQuery13)
    }


    @Test
    fun mappedChromosome() {
        // human
        assertTrue(GenomeQuery.chrDefaultChoice("chr1"))
        assertTrue(GenomeQuery.chrDefaultChoice("chr22"))
        assertTrue(GenomeQuery.chrDefaultChoice("chrX"))
        assertTrue(GenomeQuery.chrDefaultChoice("chrY"))

        assertFalse(GenomeQuery.chrDefaultChoice("chrM"))

        assertFalse(GenomeQuery.chrDefaultChoice("chr6_apd_hap1"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrUn_gl000225"))
        assertFalse(GenomeQuery.chrDefaultChoice("chr4_gl000193_random"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrUn_gl000214"))

        // D. melanogaster
        assertTrue(GenomeQuery.chrDefaultChoice("chr2L"))
        assertTrue(GenomeQuery.chrDefaultChoice("chr2LHet"))

        assertFalse(GenomeQuery.chrDefaultChoice("chrU"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrUextra"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrUn_CP007081v1"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrY_CP007108v1_random"))
        assertFalse(GenomeQuery.chrDefaultChoice("chrUn_DS484396v1"))


        // C. elegance
        assertTrue(GenomeQuery.chrDefaultChoice("chrI"))
        assertTrue(GenomeQuery.chrDefaultChoice("chrII"))
        assertTrue(GenomeQuery.chrDefaultChoice("chrV"))
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
