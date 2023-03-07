package org.jetbrains.bio.genome.search

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.search.sa.SuffixArraySearcher
import org.jetbrains.bio.genome.sequence.asNucleotideSequence
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.statistics.distribution.CategoricalDistribution
import org.junit.Assert.assertArrayEquals
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals


class SimpleSearcherTest {
    @Test
    fun testFind() {
        val pattern = "ab".asNucleotideSequence()
        val matches = SimpleSearcher("abracadabra").find(pattern)
        assertEquals("[0, 7]", Arrays.toString(matches.toArray()))
    }
}

class SuffixArraySearcherTest {

    @Test
    fun testTTTT() {
        for (chromosome in Genome["to1"].toQuery().get().reversed()) {
            val saSearcher = SuffixArraySearcher(chromosome)
            val ssSearcher = SimpleSearcher(chromosome.sequence.toString())
            val pattern = "tttt".asNucleotideSequence()
            val ssResult = ssSearcher.find(pattern).toArray()
            val saResult = saSearcher.find(pattern).sorted().toArray()
            assertArrayEquals(chromosome.toString(), ssResult, saResult)
        }
    }

    @Test
    fun testRandom() {
        for (chromosome in Genome["to1"].toQuery().get()) {
            val saSearcher = SuffixArraySearcher(chromosome)
            val ssSearcher = SimpleSearcher(chromosome.sequence.toString())
            for (k in KMER_MIN_LENGTH until KMER_MAX_LENGTH) {
                for (searchIteration in 0 until NUMBER_OF_SEARCHES) {
                    val sampleString = sampleString(charArrayOf('a', 't', 'g', 'c'), k)
                    val pattern = sampleString.asNucleotideSequence()
                    val ssResult = ssSearcher.find(pattern).toArray()
                    val saResult = saSearcher.find(pattern).sorted().toArray()
                    assertArrayEquals(
                        "Wrong answer for $k[$searchIteration]: $pattern ${chromosome.name}",
                        saResult, ssResult
                    )
                }
            }
        }
    }

    companion object {

        /**
         * Defaults `probabilities` to be a uniform distribution
         * over characters.
         */
        fun sampleString(alphabet: CharArray, length: Int): String {
            val probabilities = DoubleArray(alphabet.size)
            probabilities.fill(1.0 / alphabet.size)
            return sampleString(alphabet, length, *probabilities)
        }

        /**
         * Samples a random string in a given alphabet.
         *
         * @param alphabet letters, e.g. `"ATCG"`.
         * @param length desired length.
         * @param probabilities probabilities of all (or all but last) letters.
         * @return random string.
         */
        @JvmStatic
        fun sampleString(
            alphabet: CharArray, length: Int,
            vararg probabilities: Double
        ): String {
            val n = probabilities.size
            val d = when (alphabet.size) {
                n -> CategoricalDistribution(probabilities)
                n - 1 -> {
                    val copy = probabilities.copyOf(n - 1)
                    copy[n - 1] = 1 - probabilities.sum()
                    CategoricalDistribution(copy)
                }

                else -> throw IllegalArgumentException("missing or redundant probabilities")
            }

            val sb = StringBuilder(length)
            for (i in 0 until length) {
                sb.append(alphabet[d.sample()])
            }

            return sb.toString()
        }

        const val KMER_MIN_LENGTH: Int = 4
        const val KMER_MAX_LENGTH: Int = 20
        const val NUMBER_OF_SEARCHES: Int = 10
    }
}
