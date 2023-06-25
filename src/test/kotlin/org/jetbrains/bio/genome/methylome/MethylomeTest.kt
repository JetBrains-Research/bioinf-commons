package org.jetbrains.bio.genome.methylome

import com.google.common.math.IntMath
import org.apache.commons.math3.distribution.BinomialDistribution
import org.jetbrains.bio.dataframe.dumpHead
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.util.withTempFile
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class MethylomeTest {
    @Rule
    @JvmField
    val thrown = ExpectedException.none()

    private val genomeQuery = GenomeQuery(Genome["to1"])

    @Test
    fun testLazy() {
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

    @Test
    fun testReadWriteEmpty() {
        assertSerializedCorrectly(Methylome.builder(genomeQuery).build())
    }

    @Test
    fun testReadWriteSingleton() {
        val builder = Methylome.builder(genomeQuery)
        val chromosome = genomeQuery.get().first()
        val strand = Strand.PLUS
        builder.add(chromosome, strand, 42, CytosineContext.CG, 10, 20)
        assertSerializedCorrectly(builder.build())
    }

    @Test
    fun testReadPartialInternal() {
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
            builder.add(
                chromosome, strand, r.nextInt(10000), CytosineContext.CG,
                dC.sample(), dT.sample()
            )
        }

        val methylome = builder.build()
        assertEquals(numInfos, methylome.size)

        withTempFile("methylome", ".npz") { path ->
            methylome.save(path)

            for (chromosome in chromosomes) {
                val methylome1 = Methylome.lazy(
                    GenomeQuery(genomeQuery.genome, chromosome.name), path
                )
                assertEquals(methylome.getCombined(chromosome).size, methylome1.size)
                assertEquals(
                    methylome.getCombined(chromosome),
                    methylome1.getCombined(chromosome)
                )
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

    @Test
    fun strandedMethylome() {
        val builder = Methylome.builder(genomeQuery, stranded = false)
        val chromosome = genomeQuery.get().first()
        builder.add(chromosome, Strand.PLUS, 42, CytosineContext.CG, 10, 20)

        val methylome0 = builder.build()
        withTempFile("methylome", ".npz") { path ->
            methylome0.save(path)
            val methylome1 = Methylome.lazy(genomeQuery, path)
            assertFalse(methylome1.stranded)
        }
    }


    @Test
    fun strandIndependentMethylome() {
        val builder = Methylome.builder(genomeQuery, stranded = true)
        val chromosome = genomeQuery.get().first()
        builder.add(chromosome, Strand.PLUS, 42, CytosineContext.CG, 15, 20)
        builder.add(chromosome, Strand.PLUS, 100, CytosineContext.CHH, 1, 20)

        val methylome0 = builder.build()
        withTempFile("methylome", ".npz") { path ->
            methylome0.save(path)
            val methylome1 = Methylome.lazy(genomeQuery, path)

            assertTrue(methylome1.stranded)

            val df = methylome1[chromosome, Strand.PLUS].peel()
            assertEquals(
                "# Integer;\tByte;\tFloat;\tShort;\tShort\n" +
                        "offset\ttag\tlevel\tk\tn\n" +
                        "42\t0\t0.75\t15\t20\n" +
                        "100\t1\t0.05\t1\t20\n",
                df.dumpHead(df.rowsNumber).replace("\r", "")
            )
        }
    }

    @Test
    fun strandedMethylome_AccessMinusStrand() {
        val builder = Methylome.builder(genomeQuery, stranded = true)
        val chromosome = genomeQuery.get().first()
        builder.add(chromosome, Strand.MINUS, 42, CytosineContext.CG, 15, 20)

        assertEquals(1, builder.build()[chromosome, Strand.MINUS].size)
    }

    @Test
    fun strandIndependentMethylome_AccessMinusStrand() {
        val builder = Methylome.builder(genomeQuery, stranded = false)
        val chromosome = genomeQuery.get().first()
        builder.add(chromosome, Strand.PLUS, 42, CytosineContext.CG, 15, 20)

        thrown.expect(IllegalArgumentException::class.java)
        thrown.expectMessage("Cannot access minus strand in strand-independent methylome")

        builder.build()[chromosome, Strand.MINUS]
    }

    @Test
    fun strandIndependentMethylome_AddMinusData() {
        val builder = Methylome.builder(genomeQuery, stranded = false)
        val chromosome = genomeQuery.get().first()

        thrown.expect(IllegalArgumentException::class.java)
        thrown.expectMessage("Cannot add data to minus strand in strand-independent methylome")

        builder.add(chromosome, Strand.MINUS, 43, CytosineContext.CG, 5, 20)
    }

    @Test
    fun strandIndependentMethylome_AddCHHData() {
        val builder = Methylome.builder(genomeQuery, stranded = false)
        val chromosome = genomeQuery.get().first()

        thrown.expect(IllegalArgumentException::class.java)
        thrown.expectMessage("CHH context isn't allowed in strand-independent methylome")

        builder.add(chromosome, Strand.PLUS, 42, CytosineContext.CHH, 10, 20)
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
