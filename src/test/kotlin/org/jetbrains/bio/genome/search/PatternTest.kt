package org.jetbrains.bio.genome.search

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.junit.Test
import java.util.*
import java.util.stream.Collectors
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PatternTest {
    @Test fun testNoMismatches() {
        assertEquals("acgtacgt", seeds("acgtacgt", 0))
    }

    @Test fun testMismatches1() {
        assertEquals("acgt, acgt", seeds("acgtacgt", 1))
    }

    @Test fun testMismatches2() {
        assertEquals("ac, gta, cgt", seeds("acgtacgt", 2))
    }

    @Test fun testMismatchesAlternatives() {
        assertEquals("ata, aca, aaa, aga, ctc, ccc, cac, cgc", seeds("aNacNc", 1))
    }

    @Test fun testMismatchesTO() {
        val chromosome = Chromosome(Genome["to1"], "chr1")
        val sequence = chromosome.sequence
        for (mismatches in intArrayOf(1, 2, 3)) {
            for (length in intArrayOf(5, 10, 20)) {
                RANDOM.ints(0, chromosome.length - length).limit(1000).forEach { i ->
                    val pattern = sequence.substring(i, i + length)
                    if ('n' !in pattern) {
                        print("Length $length, mismatches $mismatches, position $i, pattern $pattern")
                        Pattern.getSeeds(pattern, mismatches).forEach { s ->
                            assertTrue(Pattern.matches(sequence, s, i, 0))
                        }
                    }
                }
            }
        }
    }

    private fun seeds(text: String, mismatches: Int): String {
        return Pattern.getSeeds(text, mismatches).flatMap { it.getAlternatives() }
            .map { it.substring(0, it.length) }
            .collect(Collectors.joining(", "))
    }

    companion object {
        private val RANDOM = Random()
    }
}
