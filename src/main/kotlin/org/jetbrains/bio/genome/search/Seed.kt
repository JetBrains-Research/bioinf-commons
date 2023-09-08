package org.jetbrains.bio.genome.search

import org.jetbrains.bio.genome.sequence.Nucleotide
import org.jetbrains.bio.genome.sequence.NucleotideAlternative
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream
import java.util.stream.Stream

/**
 * Seed for search with mismatches in extended alphabet.
 * @author Oleg Shpynov
 */
data class Seed(val text: String, val length: Int, val start: Int) {

    val nucleotideAlternatives = Array(text.length) { i ->
        NucleotideAlternative.fromChar(text[i])
    }

    /**
     * Encodes all the possible alternatives of seed [org.jetbrains.bio.genome.NucleotideAlternative] as
     * as [NucleotideSequence]
     */
    fun getAlternatives(): Stream<NucleotideSequence> {
        val alternatives = IntStream.range(start, start + length)
            .map { nucleotideAlternatives[it].alternativesCount }
            .reduce(Int::times).asInt
        return IntStream.range(0, alternatives).mapToObj<NucleotideSequence> {
            val bytes = ByteArray(length)
            var acc = it
            for (i in 0 until length) {
                val n = nucleotideAlternatives[start + i]
                val count = n.alternativesCount
                bytes[i] = n.alternatives[acc % count]
                acc /= count
            }

            object : NucleotideSequence {
                override fun charAt(pos: Int): Char = Nucleotide.ALPHABET[bytes[pos].toInt()]

                override val length: Int
                    get() = bytes.size

            }
        }
    }
}
