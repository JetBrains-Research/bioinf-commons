package org.jetbrains.bio.genome.coverage

import gnu.trove.list.TIntList
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeStrandMap
import org.slf4j.LoggerFactory
import kotlin.math.ceil
import kotlin.math.sqrt

/**
 * Estimates fragment size by cross-correlation approach.
 *
 * Reference:
 *  Kharchenko, et al. "Design and analysis of ChIP-seq experiments for DNA-binding proteins."
 *  Nature biotechnology, 2008.
 */
object FragmentSize {

    private val LOG = LoggerFactory.getLogger(FragmentSize::class.java)

    /**
     * Kharchenko et al. use this value as an upper bound, according to the figures in the article.
     * Also, "chipseq" R package uses 500 as the default upper bound.
     */
    private const val MAX_FRAGMENT_SIZE = 500

    /**
     * Detects the fragment size using the cross correlation approach.
     * We ignore candidate fragment sizes less than [averageReadLength] according to
     * Ramachandran et al., "MaSC: mappability-sensitive cross-correlation for estimating mean fragment
     * length of single-end short-read sequencing data", Bioinformatics, 2013.
     *
     * Fragment size is computed as arg max 1/N ∑Nc P[tags(c, d, +), tags(c, 0, -)]
     * where N is the total number of tags, c is the chromosome,
     * Nc is the number of tags on c, tags(c, d, s) is the boolean vector
     * of tags on chromosome c and strand s, shifted by d, P[] is the Pearson
     * correlation coefficient.
     *
     * Important: we cannot use 3rd party Pearson correlation method, because
     * total size of genome may exceed max int size, so cannot be stored in array.
     */
    fun detectFragmentSize(
        data: GenomeStrandMap<TIntList>,
        averageReadLength: Double
    ): Int {
        if (averageReadLength.isNaN()) {
            // empty data, return a placeholder value
            return 0
        }
        // Ignore phantom peaks <= read length * 1.5
        val candidateRange = ceil(averageReadLength * 1.5).toInt()..MAX_FRAGMENT_SIZE
        if (candidateRange.isEmpty()) {
            return averageReadLength.toInt()
        }
        val ccs = computePearsonCorrelationTransform(candidateRange.toList(), data)
        val fragment = ccs.maxByOrNull { it.pearsonTransform }!!.fragment
        LOG.debug("Detected fragment: $fragment")
        LOG.debug("All non-scaled cross-correlations: " +
                ccs.filter { it.pearsonTransform > 0.01 }.joinToString(",") {
                    "${it.fragment}:${String.format("%.2f", it.pearsonTransform)}"
                }
        )
        return fragment
    }

    /**
     * Stores a certain monotone substitute for Pearson correlation, see extensive comments
     * to [computePearsonCorrelationTransform].
     */
    private data class CrossCorrelation(val fragment: Int, var pearsonTransform: Double = 0.0)


    /**
     * Find the index of the next differing element in a sorted [TIntList].
     */
    private fun nextIndex(tags: TIntList, tagsSize: Int, currentTagIndex: Int): Int {
        var res = currentTagIndex + 1
        // don't assume uniqueness, just emulate it
        while (res < tagsSize && tags[res] == tags[res - 1]) {
            res++
        }
        return res
    }


    /**
     * Compute Pearson correlation transform (see below) for a given range
     * of candidate fragment sizes.
     * Pearson correlation formula
     * https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Mathematical_properties
     *
     * Linear algorithm to compute ∑xi*yi, ∑xi, ∑xi^2, ∑yi, ∑yi^2
     * for Pearson coefficient computation
     * In case of boolean unique alignment we deal with boolean vectors, where:
     *      matchedTags = ∑xi*yi
     *      positiveTags = ∑xi^2 = ∑xi
     *      negativeTags = ∑yi^2 = ∑yi
     *      Nc = ∑xi + ∑yi
     *      P[xi, yi] = (Lc∑xi*yi - ∑xi*∑yi) / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi))
     *      where Lc is the length of c.
     * The only value that depends on d is ∑xi*yi, so we can simplify the arg max.
     *      arg max 1/N ∑Nc P[tags(c, d, +), tags(c, 0, -)] =
     *      arg max ∑Nc P[tags(c, d, +), tags(c, 0, -)] =
     *      arg max ∑(∑xi + ∑yi) (Lc∑xi*yi - ∑xi*∑yi) / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi)) =
     *      arg max ∑(∑xi + ∑yi) Lc ∑xi*yi / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi)) =
     *      arg max ∑ ∑xi*yi * (∑xi + ∑yi) Lc / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi))
     *      We denote the value under arg max as "Pearson correlation transform".
     */
    private fun computePearsonCorrelationTransform(
        fragments: List<Int>,
        data: GenomeStrandMap<TIntList>
    ): List<CrossCorrelation> {
        val ccs = fragments.map { CrossCorrelation(it) }
        for (chr in data.genomeQuery.get()) {
            val positive = data[chr, Strand.PLUS]
            val positiveSize = positive.size()
            val negative = data[chr, Strand.MINUS]
            val negativeSize = negative.size()
            if (positiveSize == 0 || negativeSize == 0) continue
            val chrLength = chr.length.toDouble()
            ccs.parallelStream().forEach { cc ->
                cc.updatePearsonTransform(positive, negative, chrLength)
            }
        }
        return ccs
    }

    /**
     * Calculates a Pearson correlation transform summand for a given chromosome.
     * See [computePearsonCorrelationTransform] comments
     * for the definition of the transform and its calculation details.
     */
    private fun CrossCorrelation.updatePearsonTransform(
        positive: TIntList, negative: TIntList,
        chrLength: Double
    ) {
        val positiveSize = positive.size()
        val negativeSize = negative.size()
        var matchedTags = 0
        var positiveTags = 0
        var negativeTags = 0
        var positiveIndex = 0
        var negativeIndex = 0
        while (positiveIndex < positiveSize && negativeIndex < negativeSize) {
            val t1 = positive[positiveIndex] + fragment
            val t2 = negative[negativeIndex]
            when {
                t1 < t2 -> {
                    positiveIndex = nextIndex(positive, positiveSize, positiveIndex)
                    positiveTags++
                }

                t1 > t2 -> {
                    negativeIndex = nextIndex(negative, negativeSize, negativeIndex)
                    negativeTags++
                }

                else -> {
                    matchedTags += 1
                    positiveTags++
                    negativeTags++
                    positiveIndex = nextIndex(positive, positiveSize, positiveIndex)
                    negativeIndex = nextIndex(negative, negativeSize, negativeIndex)
                }
            }
        }
        val coefficient = (positiveSize + negativeSize) * chrLength /
                sqrt(
                    positiveSize * (chrLength - positiveSize) *
                            negativeSize * (chrLength - negativeSize)
                )
        pearsonTransform += matchedTags * coefficient
    }

}