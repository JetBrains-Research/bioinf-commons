package org.jetbrains.bio.methylome

import com.google.common.math.IntMath
import org.apache.commons.math3.distribution.BinomialDistribution
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals

class MethylomeTest {
    private val genomeQuery = GenomeQuery(Genome["to1"])

    @Test fun testLazy() {
        val builder = Methylome.builder(genomeQuery)
        val chromosome = genomeQuery.get().first()
        val strand = Strand.PLUS
        builder.add(chromosome, strand, 42, CytosineContext.CG, 10, 20)
        val methylome0 = builder.build()
        withTempFile("methylome", ".npz") { path ->
            methylome0.save(path)
            val methylome1 = Methylome.lazy(genomeQuery, path)
            assertEquals(methylome0[chromosome, strand], methylome1[chromosome, strand])
        }
    }

    @Test fun testReadWriteEmpty() {
        assertSerializedCorrectly(Methylome.builder(genomeQuery).build())
    }

    @Test fun testReadWriteSingleton() {
        val builder = Methylome.builder(genomeQuery)
        val chromosome = genomeQuery.get().first()
        val strand = Strand.PLUS
        builder.add(chromosome, strand, 42, CytosineContext.CG, 10, 20)
        assertSerializedCorrectly(builder.build())
    }

    @Test fun testReadPartialInternal() {
        val builder = Methylome.builder(genomeQuery)
        val chromosomes = genomeQuery.get()
        val numInfos = IntMath.pow(2, 16)
        val r = Random()
        val n = 42  // coverage.
        val dC = BinomialDistribution(n, .25)
        val dT = BinomialDistribution(n, .5)
        for (i in 0..numInfos - 1) {
            val chromosome = chromosomes.get(r.nextInt(chromosomes.size))
            val strand = if (r.nextBoolean()) Strand.PLUS else Strand.MINUS
            builder.add(chromosome, strand, r.nextInt(100500), CytosineContext.CG,
                        dC.sample(), dT.sample())
        }

        val methylome = builder.build()
        assertEquals(numInfos, methylome.size)

        withTempFile("methylome", ".npz") { path ->
            methylome.save(path)

            for (chromosome in chromosomes) {
                val methylome1 = Methylome.lazy(
                        GenomeQuery(genomeQuery.genome, chromosome.name), path)
                assertEquals(methylome.getCombined(chromosome).size, methylome1.size)
                assertEquals(methylome.getCombined(chromosome),
                             methylome1.getCombined(chromosome))
            }
        }
    }

    @Test
    fun duplicatedOffsets() {
        val builder = Methylome.builder(genomeQuery)
        val chromosome = genomeQuery.get().first()
        val strand = Strand.PLUS
        builder.add(chromosome, strand, 42, CytosineContext.CG, 10, 20)
        builder.add(chromosome, strand, 42, CytosineContext.CG, 11, 18)
        builder.add(chromosome, strand, 42, CytosineContext.CHH, 1, 8)
        builder.add(chromosome, strand, 2, CytosineContext.CHH, 8, 8)
        assertEquals(2, builder.duplicatedOffsets())
        assertEquals(4, builder.build()[chromosome, strand].size)
    }

    @Test
    fun duplicatedOffsetsIgnored() {
        val builder = Methylome.builder(genomeQuery, true)
        val chromosome = genomeQuery.get().first()
        val strand = Strand.PLUS
        builder.add(chromosome, strand, 42, CytosineContext.CG, 10, 20)
        builder.add(chromosome, strand, 42, CytosineContext.CG, 11, 18)
        builder.add(chromosome, strand, 42, CytosineContext.CHH, 1, 8)
        builder.add(chromosome, strand, 2, CytosineContext.CHH, 8, 8)
        assertEquals(2, builder.duplicatedOffsets())
        assertEquals(2, builder.build()[chromosome, strand].size)
    }

    private fun assertSerializedCorrectly(methylome0: Methylome) {
        withTempFile("methylome", ".npz") { path ->
            methylome0.save(path)
            val methylome1 = Methylome.lazy(genomeQuery, path)
            assertEquals(methylome0, methylome1)
        }
    }

    private fun assertEquals(m1: Methylome, m2: Methylome) {
        val genomeQuery = m1.genomeQuery
        assertEquals(genomeQuery, m2.genomeQuery)

        for (chromosome in genomeQuery.get()) {
            assertEquals(m1.getCombined(chromosome), m2.getCombined(chromosome))
        }
    }
}
