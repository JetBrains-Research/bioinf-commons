package org.jetbrains.bio.genome.sequence

import com.google.common.annotations.VisibleForTesting
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Location
import kotlin.math.max

/**
 * Computes CpG content (High, Intermediate or Low) for given [Location].
 *
 * Based on the paper: Genome-wide maps of chromatin state in pluripotent and lineage-committed cells.
 * http://www.ncbi.nlm.nih.gov/pubmed/17603471
 *
 * HCPs contain a  500-bp interval within -0.5 kb to +2 kb with
 * a (G+C)-fraction >=0.55 and a CpG observed to expected ratio (O/E) >=0.6.
 * LCPs contain no 500-bp interval with CpG O/E >=0.4.

 * @author Oleg Shpynov
 * @since 11/5/14
 */
enum class CpGContent {
    HCP, ICP, LCP;

    internal data class Counters(var cg: Int, var cpg: Int)

    companion object {
        const val MIN_LENGTH = 500

        fun classify(location: Location): CpGContent {
            return classify(location.sequence.asNucleotideSequence(), MIN_LENGTH)
        }

        fun classify(sequence: NucleotideSequence, l: Int): CpGContent {
            check(sequence.length >= l) { "Cannot classify sequences < $l" }

            val buffer = CharArray(l)
            for (i in buffer.indices) {
                buffer[i] = sequence.charAt(i)
            }

            // We are going to update these values incrementally
            var maxOE = java.lang.Double.MIN_VALUE
            var cgcpg: Counters? = null
            // We use shift at first so that we start with previous byte,
            // location.startOffset = 0 is quite impossible
            for (i in l - 1 until sequence.length) {
                if (cgcpg == null) {
                    cgcpg = Counters(computeCG(buffer), computeCpG(buffer))
                } else {
                    val byteOut = buffer[0]
                    // Shift i -> i - 1
                    System.arraycopy(buffer, 1, buffer, 0, buffer.size - 1)
                    buffer[buffer.size - 1] = sequence.charAt(i)
                    // Update cg and cpg numbers
                    update(buffer, cgcpg, byteOut, buffer[buffer.size - 1])
                }

                val observedToExpectedCpG = observedToExpected(cgcpg.cg, cgcpg.cpg)
                val cgFraction = cgFraction(l, cgcpg.cg)
                if (cgFraction >= 0.55 && observedToExpectedCpG >= 0.6) {
                    return HCP
                }
                maxOE = max(maxOE, observedToExpectedCpG)
            }
            return if (maxOE < 0.4) LCP else ICP
        }

        private fun cgFraction(length: Int, cg: Int): Double {
            return 1.0 * cg / length
        }

        /**
         * See wikipedia for more details: http://en.wikipedia.org/wiki/CpG_site
         */
        private fun observedToExpected(cg: Int, cpg: Int): Double {
            return if (cg != 0) 2.0 * cpg / cg else 0.0
        }

        internal fun update(buffer: CharArray, cgcpg: Counters, charOut: Char, charIn: Char) {
            var cg = cgcpg.cg
            var cpg = cgcpg.cpg
            // Out update
            if (charOut == 'c') {
                cg--
                if ('g' == buffer[0]) {
                    cpg--
                }
            }
            if (charOut == 'g') {
                cg--
                if ('c' == buffer[0]) {
                    cpg--
                }
            }
            // In update
            if (charIn == 'c') {
                cg++
                if ('g' == buffer[buffer.size - 2]) {
                    cpg++
                }
            }
            if (charIn == 'g') {
                cg++
                if ('c' == buffer[buffer.size - 2]) {
                    cpg++
                }
            }
            cgcpg.cg = cg
            cgcpg.cpg = cpg
        }

        @VisibleForTesting
        internal fun computeCpG(buffer: CharArray): Int {
            var cpg = 0
            for (i in 0 until buffer.size - 1) {
                val b = buffer[i]
                val nextB = buffer[i + 1]
                if ('c' == b && 'g' == nextB || 'g' == b && 'c' == nextB) {
                    cpg++
                }
            }
            return cpg
        }

        @VisibleForTesting
        internal fun computeCG(buffer: CharArray): Int {
            var cg = 0
            for (nByte in buffer) {
                if ('c' == nByte || 'g' == nByte) {
                    cg++
                }
            }
            return cg
        }

        /**
         * Slice the given chromosome into bins and return an array containing the mean GC content for each bin.
         */
        fun binnedMeanCG(chromosome: Chromosome, binSize: Int): DoubleArray {
            val sequence = chromosome.sequence
            return chromosome.range.slice(binSize).mapToDouble { bin ->
                (bin.startOffset until bin.endOffset).count { pos ->
                    sequence.charAt(pos).let { it == 'c' || it == 'g' }
                }.toDouble() / bin.length()
            }.toArray()
        }
    }
}

