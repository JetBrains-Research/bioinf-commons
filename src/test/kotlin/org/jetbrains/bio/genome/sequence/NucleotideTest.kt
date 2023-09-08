package org.jetbrains.bio.genome.sequence

import org.junit.Test
import kotlin.test.assertEquals

class NucleotideTest {
    @Test
    fun testAlphabet() {
        for (n in Nucleotide.values()) {
            assertEquals(Nucleotide.ALPHABET[n.byte.toInt()], n.char)
        }
    }

    @Test
    fun testComplementId() {
        for (n in Nucleotide.values()) {
            val ch = n.char
            assertEquals(ch, Nucleotide.complement(Nucleotide.complement(ch)))
        }
    }

    @Test
    fun testChar() {
        assertEquals('t', Nucleotide.T.char)
        assertEquals('c', Nucleotide.C.char)
        assertEquals('g', Nucleotide.G.char)
        assertEquals('a', Nucleotide.A.char)
    }

    @Test
    fun testComplement() {
        assertEquals('t', Nucleotide.complement(Nucleotide.A.char))
        assertEquals('c', Nucleotide.complement(Nucleotide.G.char))
        assertEquals(Nucleotide.C.char, Nucleotide.complement(Nucleotide.G.char))
    }

    @Test
    fun testGetChar() {
        assertEquals('t', Nucleotide.getChar(Nucleotide.T.byte))
        assertEquals('g', Nucleotide.getChar(Nucleotide.G.byte))
        assertEquals('a', Nucleotide.getChar(Nucleotide.A.byte))
        assertEquals('c', Nucleotide.getChar(Nucleotide.C.byte))
        assertEquals('n', Nucleotide.getChar(Nucleotide.ANY_NUCLEOTIDE_BYTE))
    }

    @Test
    fun testFromChar() {
        assertEquals(Nucleotide.T, Nucleotide.fromChar('t'))
        assertEquals(Nucleotide.A, Nucleotide.fromChar('a'))
        assertEquals(Nucleotide.G, Nucleotide.fromChar('g'))
        assertEquals(Nucleotide.T, Nucleotide.fromChar('T'))
        assertEquals(Nucleotide.C, Nucleotide.fromChar('C'))
        assertEquals('n', Nucleotide.getChar(Nucleotide.ANY_NUCLEOTIDE_BYTE))
    }
}
