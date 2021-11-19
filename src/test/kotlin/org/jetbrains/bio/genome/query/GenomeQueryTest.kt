package org.jetbrains.bio.genome.query

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.name
import org.jetbrains.bio.util.withTempDirectory
import org.junit.Test
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
    fun testId() {
        val gq = GenomeQuery(Genome["to1"])
        val same = GenomeQuery(gq.genome, *gq.get().map { it.name }.toTypedArray())

        assertEquals("to1[chr1,chr2,chr3,chrM,chrX]", same.id)
        assertNotEquals(gq.id, same.id)
        assertNotSame(gq, same)

        assertEquals("to1[chr1]", GenomeQuery(gq.genome, "chr1").id)
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
            val genome = Genome[chromSizesPath, "to1.${chromSizesPath.name}"]

            val genomeQuery = GenomeQuery(genome)
            val restricted = GenomeQuery(genomeQuery.genome, "chr1")

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
            assertTrue(chromosome.name.lowercase() in genomeQuery)
            assertNotNull(genomeQuery[chromosome.name.lowercase()])
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


    @Test(expected = IllegalStateException::class)
    fun parseGenomeDefinition_EmptyRestriction() {
        GenomeQuery.parseGenomeQueryId("to1[]")
    }

    @Test
    fun parseGenomeDefinition_Restrictions() {
        checkGenomeDefParsing("to1" to arrayOf("chr1", "chr3"), "to1[chr1,chr3]")
        checkGenomeDefParsing("to1" to arrayOf("chr1"), "to1[chr1]")
    }

    @Test
    fun parseGenomeDefinition_BuildName() {
        checkGenomeDefParsing("to1" to emptyArray(), "to1")
    }

    private fun checkGenomeDefParsing(expected: Pair<String, Array<String>>, genomeStr: String) {
        val (actualBuild, actualRestriction) = GenomeQuery.parseGenomeQueryId(genomeStr)
        assertEquals(expected.first, actualBuild)
        assertEquals(expected.second.joinToString(), actualRestriction.joinToString())
    }
}
