package org.jetbrains.bio.methylome

import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.sequence.asNucleotideSequence
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull

class CytosineContextTest {
    @Test fun testDetermine() {
        assertEquals(CytosineContext.CG,
                     CytosineContext.determine("atcgta".asNucleotideSequence(), 2, Strand.PLUS))
        assertNull(CytosineContext.determine("atcgta".asNucleotideSequence(), 0, Strand.PLUS))
        assertNull(CytosineContext.determine("atcgta".asNucleotideSequence(), 3, Strand.PLUS))
        assertNull(CytosineContext.determine("atcgta".asNucleotideSequence(), 2, Strand.MINUS))
        assertNull(CytosineContext.determine("atcgtat".asNucleotideSequence(), 0, Strand.MINUS))
        assertEquals(CytosineContext.CG,
                     CytosineContext.determine("atcgta".asNucleotideSequence(), 3, Strand.MINUS))

        assertEquals(CytosineContext.CHG,
                     CytosineContext.determine("atcagta".asNucleotideSequence(), 2, Strand.PLUS))
        assertEquals(CytosineContext.CHG,
                     CytosineContext.determine("atctgta".asNucleotideSequence(), 2, Strand.PLUS))
        assertEquals(CytosineContext.CHG,
                     CytosineContext.determine("atccgta".asNucleotideSequence(), 2, Strand.PLUS))
        assertEquals(CytosineContext.CHG,
                     CytosineContext.determine("tcgggta".asNucleotideSequence(), 3, Strand.MINUS))

        assertEquals(CytosineContext.CHH,
                     CytosineContext.determine("atcaata".asNucleotideSequence(), 2, Strand.PLUS))
        assertEquals(CytosineContext.CHH,
                     CytosineContext.determine("atcacta".asNucleotideSequence(), 2, Strand.PLUS))
        assertEquals(CytosineContext.CHH,
                     CytosineContext.determine("atgggta".asNucleotideSequence(), 3, Strand.MINUS))
    }
}
