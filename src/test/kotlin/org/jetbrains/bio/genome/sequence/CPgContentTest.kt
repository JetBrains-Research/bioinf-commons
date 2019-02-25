package org.jetbrains.bio.genome.sequence

import org.junit.Test
import kotlin.test.assertEquals

class CpGContentTest {
    @Test fun testComputeCpG() {
        assertEquals(1, CpGContent.computeCpG("cacgg".toCharArray()))
    }

    @Test fun testComputeCG() {
        assertEquals(2, CpGContent.computeCG("acg".toCharArray()))
    }

    @Test fun testClassify() {
        assertEquals(CpGContent.LCP, CpGContent.classify("AAAAAA".asNucleotideSequence(), 3))
        assertEquals(CpGContent.LCP, CpGContent.classify("CCCCCC".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 5))
        assertEquals(CpGContent.ICP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 8))
        assertEquals(CpGContent.ICP, CpGContent.classify("CCCCGCCC".asNucleotideSequence(), 8))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCGCG".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CCCGCG".asNucleotideSequence(), 5))
        assertEquals(CpGContent.HCP, CpGContent.classify("CGCGCG".asNucleotideSequence(), 5))
    }
}
