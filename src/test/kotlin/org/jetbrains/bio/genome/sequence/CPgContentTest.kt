package org.jetbrains.bio.genome.sequence

import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.toQuery
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CpGContentTest {
    @Test
    fun testComputeCpG() {
        assertEquals(1, CpGContent.computeCpG("cacgg".toCharArray()))
    }

    @Test
    fun testComputeCG() {
        assertEquals(2, CpGContent.computeCG("acg".toCharArray()))
    }

    @Test
    fun testClassify() {
        assertEquals(CpGContent.LCP, CpGContent.classify("AAAAAA".asNucleotideSequence(), 3))
        assertEquals(CpGContent.LCP, CpGContent.classify("CCCCCC".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 5))
        assertEquals(CpGContent.ICP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 8))
        assertEquals(CpGContent.ICP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 8))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCGCG".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCGCG".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CGCGCG".asNucleotideSequence(), 5))
    }

    @Test
    fun testBinnedMeanCG() {
        val chr = Genome["to1"].toQuery()["chr1"]!!
        val binned100 = CpGContent.binnedMeanCG(chr, 100)
        val binned200 = CpGContent.binnedMeanCG(chr, 200)
        assertEquals((chr.length - 1) / 100 + 1, binned100.size)
        assertEquals((chr.length - 1) / 200 + 1, binned200.size)
        binned200.forEachIndexed { index, value ->
            assertTrue(value >= 0, "Mean CG content was negative ($value at bin $index)")
            assertTrue(value <= 1, "Mean CG content was greater than 1 ($value at bin $index)")
            if (index < binned200.size - 1) {
                Tests.assertEquals(
                    (binned100[2 * index] + binned100[2 * index + 1]) / 2,
                    value,
                    1E-10, // we want to ignore double arithmetic precision issues
                    "Expected mean CG $value at 200 bps bin $index to be equal to mean of 100bps bins " +
                            "${2 * index} (${binned100[2 * index]}) and " +
                            "${2 * index + 1} (${binned100[2 * index + 1]})."
                )
            }
        }
    }
}
