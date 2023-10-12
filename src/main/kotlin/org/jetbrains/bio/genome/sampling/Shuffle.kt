package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.containers.intersectMap
import org.jetbrains.bio.gsea.IntHistogram
import kotlin.math.max


/**
 * Fast reimplementation of `bedtools shuffle`.
 * Uniformly shuffles regions in such way that they starts are in background. Resulting region will be not intersecting.
 * Such behaviour allows fast implementations. Also it is consistent with `bedtools shuffle`.
 * @param genomeQuery genome to use
 * @param regions regions to shuffle. Only length of regions are used.
 * @param background background regions, if null shuffle from whole genome.
 * @param maxRetries number of attempts to generate regions satisfying all conditions.
 * @param withReplacement If true then resulted regions could intersect.
 */

fun shuffleChromosomeRanges(
    genomeQuery: GenomeQuery,
    regions: List<ChromosomeRange>,
    background: List<ChromosomeRange>? = null,
    regionSetMaxRetries: Int = 100,
    singleRegionMaxRetries: Int = 100,
    withReplacement: Boolean,
    genomeMaskedAreaFilter: LocationsMergingList?,
    randomGen: kotlin.random.Random  // Need Random.nextLong(bound) method that is available in Java 17 or in Kotlin
): Pair<List<ChromosomeRange>, IntHistogram> {
    val lengths = regions.map { it.length() }.toIntArray()

    val backgroundRegions = background ?: genomeQuery.get().map {
        ChromosomeRange(0, it.length, it)
    }

    val prefixSum = createPrefixSumFor(backgroundRegions)

    for (i in 1..regionSetMaxRetries) {
        val (sampledRegions, attemptsHistogram) = tryShuffle(
            genomeQuery, backgroundRegions, lengths, prefixSum, singleRegionMaxRetries, withReplacement,
            genomeMaskedAreaFilter, randomGen
        )
        if (sampledRegions != null) {
            // TODO: report max retries!!!!!!!!!!!!!!!!!!!!!!!!!!!
            return sampledRegions to attemptsHistogram
        }
            // TODO: report max retries!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // TODO: report max retries!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // TODO: report max retries!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // TODO: report max retries!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    throw RuntimeException("Too many shuffle attempts. Max limit is: $regionSetMaxRetries")
}

fun createPrefixSumFor(backgroundRegions: List<ChromosomeRange>): LongArray {
    val prefixSum = LongArray(backgroundRegions.size)

    var s = 0L

    for (i in backgroundRegions.indices) {
        s += backgroundRegions[i].length()
        prefixSum[i] = s
    }
    return prefixSum
}

fun tryShuffle(
    gq: GenomeQuery,
    background: List<ChromosomeRange>,
    lengths: IntArray,
    prefixSum: LongArray,
    singleRegionMaxRetries: Int = 100,
    withReplacement: Boolean,
    genomeMaskedAreaFilter: LocationsMergingList?,
    randomGen: kotlin.random.Random  // Need Random.nextLong(bound) method that is available in Java 17 or in Kotlin
): Pair<List<ChromosomeRange>?, IntHistogram> {
    // Trick: used one storage for masked & seen intervals that is forbidden to intersect
    val maskedArea: GenomeMap<MutableList<Range>>? = createMaskedArea(genomeMaskedAreaFilter, withReplacement, gq)

    val result = ArrayList<ChromosomeRange>(lengths.size)
    val attemptsHistogram = IntHistogram()

    val sum = prefixSum.last()

    for (expectedLen in lengths) {
        require(expectedLen > 0) { "Empty region not supported" }
        var nextRange: ChromosomeRange? = null

        var lastRangeTry = 0
        for (rangeTry in 1..singleRegionMaxRetries) {
            lastRangeTry = rangeTry
            val rVal = randomGen.nextLong(sum)

            val index = prefixSum.binarySearch(rVal)

            val j = if (index >= 0) {
                index + 1
            } else {
                -index - 1
            }

            val l = background[j]
            val offset = l.length() + (rVal - prefixSum[j])  // random offset in selected random interval
            val anchorStartOffset = (l.startOffset + offset).toInt()
            // [Not BED Tools]: just improvement so not to start a sampled region always directly in
            // a background position, but just intersect background position, called 'anchor'
            // Example: exp_len = 20
            //  * anchor=5 =>
            //  * anchor=50 => [31,51)..[50, 70): shift=0..19=randInt(exp_len) => [anchor - shift, anchor - shift + exp_len]
            //  * anchor=5 => [0,20) .. [5, 25): shift=0..5=0..min(anchor, exp_len)
            // [anchor=5, exp_len=20] => [0, 20) ..
            // [5, 25) : shift=0..5=0..min(anchor, exp_len)

            // XXX: [option 1] like default BED Tools behavior
            // val startOffset = anchorStartOffset

            // XXX: [option 2] like default BED Tools behavior

            val startOffset = max(0,  anchorStartOffset - randomGen.nextInt(expectedLen))

            // Make Range
            val candidate = ChromosomeRange(startOffset, startOffset + expectedLen, l.chromosome)
            if (candidate.endOffset <= candidate.chromosome.length) {
                val isCandidateValid = if (maskedArea == null) {
                    // All candidates are valid, i.e., sampling with replacement and not masked genome
                    true
                } else {
                    // Doesn't intersect already seen regions or initially masked genome
                    !candidate.intersectMap(maskedArea)
                }

                if (isCandidateValid) {
                    // success
                    nextRange = candidate
                    break
                }
            }
            // else retry
        }
        attemptsHistogram.increment(lastRangeTry)

        if (nextRange == null) {
            // Cannot sample next range in given attempts number
            return null to attemptsHistogram
        }
        if (!withReplacement) {
            maskedArea!![nextRange.chromosome].add(nextRange.toRange())
        }
        result.add(nextRange)
    }

    return result to attemptsHistogram
}

fun createMaskedArea(
    genomeMaskedAreaFilter: LocationsMergingList?,
    withReplacement: Boolean,
    gq: GenomeQuery
): GenomeMap<MutableList<Range>>? {
    // Trick: used one storage for masked & seen intervals that is forbidden to intersect

    val maskedArea: GenomeMap<MutableList<Range>>? = if (genomeMaskedAreaFilter == null && withReplacement) {
        null
    } else {
        val filter = genomeMaskedAreaFilter ?: LocationsMergingList.create(gq, emptyList())

        genomeMap(gq) {
            require(filter[it, Strand.MINUS].isEmpty()) { "Sampling ignores strand. Minus strand isn't supported." }
            filter[it, Strand.PLUS].map { it.toRange() }.toMutableList()
        }
    }
    return maskedArea
}

