package org.jetbrains.bio.genome.sequence

import org.junit.Test
import java.lang.IllegalArgumentException
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull

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
    fun testFromByte() {
        assertEquals(Nucleotide.T, Nucleotide.fromByte(Nucleotide.T.byte))
        assertEquals(Nucleotide.A, Nucleotide.fromByte(Nucleotide.A.byte))
        assertEquals(Nucleotide.C, Nucleotide.fromByte(Nucleotide.C.byte))
        assertEquals(Nucleotide.G, Nucleotide.fromByte(Nucleotide.G.byte))

        assertNull(Nucleotide.fromByte(Nucleotide.ANY_NUCLEOTIDE_BYTE))
    }

    @Test
    fun testFromInvalidByte(){
        assertFailsWith(IndexOutOfBoundsException::class, message = "invalid nucleotide byte = 6: (6) must be less than size (4)") {
            assertEquals(Nucleotide.G, Nucleotide.fromByte(6.toByte()))
        }
        assertFailsWith(IndexOutOfBoundsException::class, message = "invalid nucleotide byte = 128: (128) must be less than size (4)") {
            assertEquals(Nucleotide.G, Nucleotide.fromByte(128.toByte()))
        }
    }

    @Test
    fun testFromChar() {
        assertEquals(Nucleotide.T, Nucleotide.fromChar('t'))
        assertEquals(Nucleotide.A, Nucleotide.fromChar('a'))
        assertEquals(Nucleotide.G, Nucleotide.fromChar('g'))
        assertEquals(Nucleotide.T, Nucleotide.fromChar('T'))
        assertEquals(Nucleotide.C, Nucleotide.fromChar('C'))

        assertNull( Nucleotide.fromChar('N'))
        assertNull( Nucleotide.fromChar('n'))
        assertNull( Nucleotide.fromChar('?'))
    }
}
