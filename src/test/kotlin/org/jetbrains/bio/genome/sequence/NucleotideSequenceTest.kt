package org.jetbrains.bio.genome.sequence

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Strand
import org.junit.Test
import kotlin.test.assertEquals

class NucleotideSequenceTest {
    @Test fun testSubstringBasic() {
        val sequence = "aaccggtt".asNucleotideSequence()
        assertEquals("ccggt", sequence.substring(2, 7))
        assertEquals("ccggt", sequence.substring(2, 7, Strand.PLUS))
        assertEquals("accgg", sequence.substring(2, 7, Strand.MINUS))
    }

    @Test fun testSubstringTestOrganism() {
        val chromosome = Chromosome("to1", "chr1")
        val sequence = chromosome.sequence
        val from = sequence.length / 3
        val to = sequence.length - from
        testSubstring(sequence, from, to, Strand.PLUS)
        testSubstring(sequence, from, to, Strand.MINUS)
    }

    private fun testSubstring(sequence: NucleotideSequence, from: Int,
                              to: Int, strand: Strand) {
        val sub = sequence.substring(from, to, strand)
        for (offset in from..to - 1) {
            val expected = sequence.charAt(offset)
            val actual = sub[if (strand.isMinus()) to - offset - 1 else offset - from]
            if (expected == Nucleotide.N && actual == Nucleotide.N) {
                continue
            }

            if (strand.isMinus()) {
                // Now it is obvious how ugly this is.
                assertEquals(expected, Nucleotide.complement(actual))
            } else {
                assertEquals(expected, actual)
            }
        }
    }

    @Test fun byteAtCasing() {
        // XXX currently 'NucleotideSequence' doesn't distinguish between
        // character casing.
        val sequence = "aacCgGtT".asNucleotideSequence()
        assertEquals('a', sequence.charAt(0))
        assertEquals('c', sequence.charAt(2))
        assertEquals('c', sequence.charAt(3))

        assertEquals('c', sequence.charAt(2, Strand.PLUS))
        assertEquals('c', sequence.charAt(3, Strand.PLUS))

        assertEquals('t', sequence.charAt(0, Strand.MINUS))
        assertEquals('g', sequence.charAt(3, Strand.MINUS))
        assertEquals('c', sequence.charAt(4, Strand.MINUS))
        assertEquals('a', sequence.charAt(6, Strand.MINUS))
    }
}
