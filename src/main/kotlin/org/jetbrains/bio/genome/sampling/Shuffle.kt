package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import java.util.concurrent.ThreadLocalRandom


/**
 * Fast reimplementation of `bedtools shuffle`.
 * Uniformly shuffles regions in such way that they starts are in background. Resulting region will be not intersecting.
 * Such behaviour allows fast implementations. Also it is consistent with `bedtools shuffle`.
 * @param genomeQuery genome to use
 * @param regions regions to shuffle. Only length of regions are used.
 * @param background background regions, if null shuffle from whole genome.
 * @param maxRetries number of attempts to generate regions satisfying all conditions.
 */

fun shuffleChromosomeRanges(genomeQuery: GenomeQuery,
                            regions: List<ChromosomeRange>,
                            background: List<ChromosomeRange>? = null,
                            maxRetries: Int = 100): List<ChromosomeRange> {
    val lengths = regions.map { it.length() }.toIntArray()

    val backgroundRegions = background ?: genomeQuery.get().map {
        ChromosomeRange(0, it.length, it)
    }

    val prefixSum = LongArray(backgroundRegions.size)

    var s = 0L

    for (i in backgroundRegions.indices) {
        s += backgroundRegions[i].length()
        prefixSum[i] = s
    }

    for (i in 1..maxRetries) {
        val result = tryShuffle(genomeQuery, backgroundRegions, lengths, prefixSum, maxRetries)
        if (result != null) {
            return result
        }
    }
    throw RuntimeException("Too many shuffle attempts")
}

private fun hasIntersection(regionsMap: GenomeMap<MutableList<Range>>, chromosomeRange: ChromosomeRange): Boolean {
    val list = regionsMap[chromosomeRange.chromosome]
    val range = chromosomeRange.toRange()
    for (r in list) {
        if (r intersects range) {
            return true
        }
    }
    return false
}

private fun tryShuffle(genomeQuery: GenomeQuery,
                       background: List<ChromosomeRange>,
                       lengths: IntArray,
                       prefixSum: LongArray,
                       maximalReTries: Int = 100): List<ChromosomeRange>? {
    val regionMap = genomeMap<MutableList<Range>>(genomeQuery) {
        ArrayList()
    }
    val result = ArrayList<ChromosomeRange>(lengths.size)

    val r = ThreadLocalRandom.current()

    val sum = prefixSum.last()

    for (i in lengths.indices) {
        var nextRange: ChromosomeRange? = null

        for (rangeTry in 1..maximalReTries) {
            val rVal = r.nextLong(sum)

            val index = prefixSum.binarySearch(rVal)

            val j = if (index >= 0) {
                index + 1
            } else {
                -index - 1
            }

            val l = background[j]
            val offset = l.length() + (rVal - prefixSum[j])

            val startOffset = (l.startOffset + offset).toInt()

            val candidate = ChromosomeRange(startOffset, startOffset + lengths[i], l.chromosome)
            if (candidate.endOffset <= candidate.chromosome.length && !hasIntersection(regionMap, candidate)) {
                // success
                nextRange = candidate
                break
            }
            // else retry
        }

        if (nextRange == null) {
            // Cannot sample next range in given attempts number
            return null
        }
        regionMap[nextRange.chromosome].add(nextRange.toRange())
        result.add(nextRange)
    }

    return result
}
