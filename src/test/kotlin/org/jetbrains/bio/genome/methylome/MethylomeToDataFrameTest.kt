package org.jetbrains.bio.genome.methylome

import org.jetbrains.bio.dataframe.toBitterSet
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.StrandFilter
import org.junit.Assert.assertArrayEquals
import org.junit.Test
import kotlin.test.assertEquals

class MethylomeToDataFrameTest {
    @Test fun single() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        val chromosome = genomeQuery.get().first()
        val methylome = Methylome.builder(genomeQuery)
                .add(chromosome, Strand.PLUS, 100, CytosineContext.CG, 1, 10)
                .add(chromosome, Strand.MINUS, 200, CytosineContext.CG, 2, 10)
                .add(chromosome, Strand.PLUS, 200, CytosineContext.CHG, 2, 10)
                .build()

        val df = MethylomeToDataFrame.create(StrandFilter.BOTH, CytosineContext.CG)
                .apply(methylome, chromosome)
        assertEquals(2, df.rowsNumber)
        assertArrayEquals(intArrayOf(100, 200), df.sliceAsInt("offset"))
        assertArrayEquals(shortArrayOf(1, 2), df.sliceAsShort("k"))
        assertArrayEquals(shortArrayOf(10, 10), df.sliceAsShort("n"))
        assertArrayEquals(byteArrayOf(CytosineContext.CG.tag,
                                      CytosineContext.CG.tag),
                          df.sliceAsByte("tag"))
        assertArrayEquals(intArrayOf(Int.MAX_VALUE, 100), df.sliceAsInt("d"))
        assertEquals(
                booleanArrayOf(true, false).toBitterSet(),
                df.sliceAsBool("strand")
        )
    }

    @Test fun double() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        val chromosome = genomeQuery.get().first()
        val methylome = Methylome.builder(genomeQuery)
                .add(chromosome, Strand.PLUS, 100, CytosineContext.CG, 1, 10)
                .build()

        val df = MethylomeToDataFrame.DEFAULT
                .apply(listOf(methylome, methylome), chromosome)
        assertEquals(1, df.rowsNumber)
        assertArrayEquals(byteArrayOf(CytosineContext.CG.tag),
                          df.sliceAsByte("tag"))
    }

    @Test fun triple() {
        val genomeQuery = GenomeQuery(Genome["to1"])
        val chromosome = genomeQuery.get().first()
        val methylome = Methylome.builder(genomeQuery)
                .add(chromosome, Strand.PLUS, 100, CytosineContext.CG, 1, 10)
                .build()

        val df = MethylomeToDataFrame.create(StrandFilter.BOTH, CytosineContext.CG)
                .apply(listOf(methylome, methylome, methylome), chromosome)
        assertEquals(1, df.rowsNumber)

        assertEquals(
                booleanArrayOf(true).toBitterSet(),
                df.sliceAsBool("strand")
        )
        assertArrayEquals(byteArrayOf(CytosineContext.CG.tag),
                          df.sliceAsByte("tag"))
    }
}
