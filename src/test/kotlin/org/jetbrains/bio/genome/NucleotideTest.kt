package org.jetbrains.bio.genome

import org.jetbrains.bio.genome.sequence.Nucleotide
import org.junit.Test
import kotlin.test.assertEquals

class NucleotideTest {
    @Test fun testAlphabet() {
        for (n in Nucleotide.values()) {
            assertEquals(Nucleotide.ALPHABET[n.byte.toInt()], n.char)
        }
    }

    @Test fun testComplementId() {
        for (n in Nucleotide.values()) {
            val ch = n.char
            assertEquals(ch, Nucleotide.complement(Nucleotide.complement(ch)))
        }
    }

    @Test fun testComplement() {
        assertEquals(Nucleotide.T.char, Nucleotide.complement(Nucleotide.A.char))
        assertEquals(Nucleotide.C.char, Nucleotide.complement(Nucleotide.G.char))
    }
}
