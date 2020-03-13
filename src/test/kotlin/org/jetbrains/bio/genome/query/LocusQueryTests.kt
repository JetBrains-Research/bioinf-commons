package org.jetbrains.bio.genome.query

import org.junit.Assert.assertNotSame
import org.junit.Test
import kotlin.test.assertEquals

class LocusQueryTest {
    @Test fun testEquals() {
        assertEquals(TranscriptQuery(), TranscriptQuery())
        assertEquals(TssQuery(), TssQuery(-2000, 2000))
        assertNotSame(TssQuery(), TssQuery(-2000, 2001))
        assertEquals(TesQuery(), TesQuery(-2000, 2000))
        assertNotSame(TesQuery(), TesQuery(-2000, 2001))
    }

    @Test fun testParseLoci() {
        for (locus in arrayOf(
                TssQuery(-5000, -2500), // [1] Distal enhancer
                TssQuery(), // [3, 4, 7] Default ±2kbp
                UTR5Query(), // [4]
                UTR3Query(), // [4]
                CDSQuery(), // [4]
                IntronsQuery(), // [4]
                ExonsQuery(), // [4]
                TranscriptQuery(), // Important for H3K36me3
                TesQuery(), // Default ±2kbp
                TesQuery(2500, 5000) // Distal TES
        )) {
            assertEquals(locus.id, toString(locus.id))
        }
    }

    @Test fun testParseTssWidth() {
        assertEquals("tss[-500..500]", toString("tss500"))
        assertEquals("tes[-500..500]", toString("tes500"))
        assertEquals("tss[-500..500]", toString("tss[500]"))
    }

    @Test fun testParseTssWithBoundaries() {
        assertEquals("tss[-500..500]", toString("tss[-500..500]"))
        assertEquals("tes[-500..500]", toString("tes[-500..500]"))
    }

    private fun toString(locus: String) = TranscriptLocusQuery.parse(locus)!!.id
}