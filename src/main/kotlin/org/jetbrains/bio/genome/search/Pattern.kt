package org.jetbrains.bio.genome.search

import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream
import java.util.stream.Stream

/**
 * Pattern, supporting extended nucleotide alphabet.
 *
 * @see org.jetbrains.bio.genome.sequence.NucleotideAlternative
 */
object Pattern {
    @JvmStatic fun matches(sequence: NucleotideSequence, seed: Seed, position: Int,
                           mismatches: Int): Boolean {
        val length = seed.nucleotideAlternatives.size
        if (position < 0 || position + length >= sequence.length) {
            return false
        }
        var mm = 0
        // Before seed check
        if (seed.start > 0) {
            for (i in 0 until seed.start) {
                if (!seed.nucleotideAlternatives[i].match(sequence.charAt(position + i))) {
                    if (++mm > mismatches) {
                        return false;
                    }
                }
            }
        }
        // After seed check
        val seedEnd = seed.start + seed.length
        if (seedEnd < length) {
            for (i in seedEnd until length) {
                if (!seed.nucleotideAlternatives[i].match(sequence.charAt(position + i))) {
                    if (++mm > mismatches) {
                        return false;
                    }
                }
            }
        }
        return true
    }

    /**
     * Calculates and returns the seeds for the given pattern.
     *
     * The seeds should be non-intersecting and their number should be
     * strictly greater than the maximum mismatch count allowed.
     */
    @JvmStatic fun getSeeds(text: String, mismatches: Int): Stream<Seed> {
        val seed = 1f * text.length / (mismatches + 1)
        return IntStream.range(0, mismatches + 1).mapToObj { i ->
            val seedStart = (i * seed).toInt();
            val seedEnd = ((i + 1) * seed).toInt();
            Seed(text, seedEnd - seedStart, seedStart);
        }
    }
}
