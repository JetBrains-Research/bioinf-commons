package org.jetbrains.bio.genome.search

import htsjdk.samtools.util.SequenceUtil
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream
import java.util.stream.Stream

/**
 * The chromosome searcher which allows searching with a fixed number of mismatches allowed.
 * @author: Alexey Dievsky
 */
class ChromosomeSearcher(private val searcher: Searcher,
                         private val mismatches: Int,
                         private val chromosome: Chromosome) {

    companion object {
        const val UNIQUE = -3
        const val NOT_UNIQUE = -2
        const val NOT_FOUND = -1
    }

    private val sequence: NucleotideSequence = chromosome.sequence

    /**
     * Searches for the occurrences of the pattern on the chromosome.
     * @param text The search pattern.
     */
    fun find(text: String): Stream<Location> {
        return Stream.concat(
            findOnPositiveStrand(text).mapToObj<Location> {
                Location(it, it + text.length, chromosome, Strand.PLUS)
            },
            findOnPositiveStrand(SequenceUtil.reverseComplement(text)).mapToObj<Location> {
                Location(it, it + text.length, chromosome, Strand.MINUS)
            })
    }

    private fun findOnPositiveStrand(text: String): IntStream {
        return Pattern.getSeeds(text, mismatches).flatMapToInt { seed ->
            seed.getAlternatives().flatMapToInt {
                searcher.find(it)
                    .map { it - seed.start }
                    .filter { Pattern.matches(sequence, seed, it, mismatches) }
            }
        }.distinct()
    }

    fun findUnique(text: String): Location? {
        val positive = findOnPositiveStrandUnique(text)
        if (positive == NOT_UNIQUE) {
            return null
        }
        val negative = findOnPositiveStrandUnique(SequenceUtil.reverseComplement(text))
        if (negative == NOT_UNIQUE) {
            return null
        }
        if (positive != NOT_FOUND && negative == NOT_FOUND) {
            return Location(positive, positive + text.length, chromosome, Strand.PLUS)
        }
        if (negative != NOT_FOUND && positive == NOT_FOUND) {
            return Location(negative, negative + text.length, chromosome, Strand.MINUS)
        }
        return null
    }

    /**
     * Returns one of [UNIQUE], [NOT_UNIQUE], [NOT_FOUND]
     */
    fun test(text: String): Int {
        val positive = findOnPositiveStrandUnique(text)
        if (positive == NOT_UNIQUE) {
            return NOT_UNIQUE
        }
        val negative = findOnPositiveStrandUnique(SequenceUtil.reverseComplement(text))
        if (negative == NOT_UNIQUE) {
            return NOT_UNIQUE
        }
        if (positive != NOT_FOUND && negative == NOT_FOUND ||
            negative != NOT_FOUND && positive == NOT_FOUND) {
            return UNIQUE
        }
        return NOT_FOUND
    }

    /**
     * Searches for the occurrences of the pattern on the positive strand.
     * @param text The search pattern.
     * @return int value:
     *  [NOT_FOUND] in case when nothing found,
     *  [NOT_UNIQUE] in case of multiple
     *  index in case of unique
     */
    private fun findOnPositiveStrandUnique(text: String): Int {
        val iterator = Pattern.getSeeds(text, mismatches)
            .flatMapToInt { seed ->
                seed.getAlternatives()
                    .flatMapToInt {
                        searcher.find(it)
                            .map { it - seed.start }
                            .filter { Pattern.matches(sequence, seed, it, mismatches) }
                    }
            }
            .distinct().iterator()
        if (!iterator.hasNext()) {
            return NOT_FOUND
        }
        val index = iterator.nextInt()
        return if (iterator.hasNext()) NOT_UNIQUE else index
    }
}
