package org.jetbrains.bio.genome.coverage

import com.google.common.base.MoreObjects
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import gnu.trove.set.hash.TIntHashSet
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeStrandMap
import org.jetbrains.bio.genome.containers.genomeStrandMap
import org.jetbrains.bio.npy.NpzFile
import java.io.IOException
import java.nio.file.Path
import kotlin.math.sqrt

/**
 * The container maintains a sorted list of tag offsets for each chromosome and strand.
 * Immutable. Saves data in [NpzFile] format.
 *
 * [detectedFragment] is estimated from the data using cross-correlation and saved
 * in the cache file.
 * [actualFragment] is used when computing location coverage and can be altered
 * by the user using [withFragment] method. By default it's equal to [detectedFragment].
 *
 * @author Oleg Shpynov
 * @author Aleksei Dievskii
 */
class SingleEndCoverage private constructor(
        override val genomeQuery: GenomeQuery,
        val detectedFragment: Int,
        val actualFragment: Int = detectedFragment,
        internal val data: GenomeStrandMap<TIntList>
): Coverage {

    override fun getCoverage(location: Location) = getTags(location).size

    override val depth = genomeQuery.get().flatMap { chr ->
        Strand.values().map { strand ->
            data[chr, strand].size().toLong()
        }
    }.sum()

    /**
     * Returns a sorted array of reads covered by a given [location].
     */
    internal fun getTags(location: Location): IntArray {
        val data = data[location.chromosome, location.strand]
        /* we don't really care if the offsets are outside of chromosome range,
           since this won't lead to incorrect results */
        val locationShift = location.strand.choose(
                (-actualFragment) / 2,
                actualFragment / 2
        )
        val startOffset = location.startOffset + locationShift
        val endOffset = location.endOffset + locationShift
        val index = data.binarySearchLeft(startOffset)
        var size = 0
        while (index + size < data.size() &&
                data[index + size] < endOffset) {
            size++
        }

        return data.toArray(index, size)
    }

    @Throws(IOException::class)
    internal fun save(outputPath: Path) {
        NpzFile.write(outputPath).use { writer ->
            writer.write(Coverage.VERSION_FIELD, intArrayOf(Coverage.VERSION))
            writer.write(Coverage.PAIRED_FIELD, booleanArrayOf(false))
            writer.write(FRAGMENT_FIELD, intArrayOf(detectedFragment))

            for (chromosome in genomeQuery.get()) {
                for (strand in Strand.values()) {
                    val key = chromosome.name + '/' + strand
                    writer.write(key, data[chromosome, strand].toArray())
                }
            }
        }
    }

    override fun toString() = MoreObjects.toStringHelper(this)
            .addValue(genomeQuery).toString()

    /**
     * Sets fragment size to specified value or reverts to a detected one if "null" is provided.
     * Returns a copy of the immutable [SingleEndCoverage] object.
     */
    fun withFragment(fragment: Fragment): SingleEndCoverage = if (fragment is FixedFragment) {
        SingleEndCoverage(genomeQuery, detectedFragment, fragment.size, data)
    } else {
        SingleEndCoverage(genomeQuery, detectedFragment, data = data)
    }

    /**
     * Stores a certain monotone substitute for Pearson correlation, see extensive comments
     * to [computePearsonCorrelationTransform].
     */
    private data class CrossCorrelation(val fragment: Int, var pearsonTransform: Double = 0.0)

    class Builder(val genomeQuery: GenomeQuery) {

        val data: GenomeStrandMap<TIntList> = genomeStrandMap(genomeQuery) {
            _, _ -> TIntArrayList()
        }

        private var readLengthSum = 0L
        private var readCount = 0L

        /**
         * Add a tag to the coverage being built. Only the 5' end of [read] is relevant.
         */
        fun process(read: Location): Builder {
            data[read.chromosome, read.strand].add(read.get5Bound())
            readLengthSum += read.length()
            readCount++
            return this
        }

        /**
         * Generate a [SingleEndCoverage] object.
         * [unique] controls whether duplicate tags should be preserved ([unique] == false)
         * or squished into one tag ([unique] == true).
         * Only tags at the exact same offset on the exact same strand
         * are considered duplicate.
         * [detectedFragment] is imputed at this point.
         */
        fun build(unique: Boolean): SingleEndCoverage {
            for (chromosome in data.genomeQuery.get()) {
                for (strand in Strand.values()) {
                    if (unique) {
                        // XXX we can do linear time de-duplication on
                        //     the sorted sequence.
                        val tags = data[chromosome, strand]
                        data[chromosome, strand] = TIntArrayList(TIntHashSet(tags))
                    }
                    data[chromosome, strand].sort()
                }
            }

            val detectedFragment = detectFragmentSize(
                    data,
                    readLengthSum * 1.0 / readCount
            )

            return SingleEndCoverage(genomeQuery, detectedFragment, data = data)
        }
    }

    companion object {

        const val FRAGMENT_FIELD = "fragment"

        /**
         * Kharchenko et al. (see link below) use this value as an upper bound,
         * judging by the figures in the article.
         * Also, "chipseq" R package uses 500 as the default upper bound.
         */
        private const val MAX_FRAGMENT_SIZE = 500

        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        internal fun load(
                npzReader: NpzFile.Reader,
                genomeQuery: GenomeQuery
        ): SingleEndCoverage {
            check(!npzReader[Coverage.PAIRED_FIELD].asBooleanArray().single()) {
                "${npzReader.path} attempting to read single-end coverage from paired-end cache file"
            }
            val detectedFragment = npzReader[FRAGMENT_FIELD].asIntArray().single()
            val data: GenomeStrandMap<TIntList> = genomeStrandMap(genomeQuery) { _, _ ->
                TIntArrayList()
            }
            for (chromosome in genomeQuery.get()) {
                for (strand in Strand.values()) {
                    val key = chromosome.name + '/' + strand
                    try {
                        val npyArray = npzReader[key]
                        data[chromosome, strand] = TIntArrayList.wrap(npyArray.asIntArray())
                    } catch (e: IllegalStateException) {
                        throw IllegalStateException(
                                "Cache file ${npzReader.path} doesn't contain $key.\n" +
                                        "It's likely that chrom.sizes file used for its creation differs " +
                                        "from the one being used to read it now.\n" +
                                        "If problem persists, delete the cache file ${npzReader.path} " +
                                        "and Span will recreate it with correct settings.",
                                e
                        )
                    }
                }
            }
            return SingleEndCoverage(genomeQuery, detectedFragment, data = data)
        }

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
         */
        private fun computePearsonCorrelationTransform(
                fragments: List<Int>,
                data: GenomeStrandMap<TIntList>
        ): List<CrossCorrelation> {
            /*
                Using the approach from Kharchenko et. al, 2008
                    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2597701/
                we compute
                    arg max 1/N ∑Nc P[tags(c, d, +), tags(c, 0, -)]
                where N is the total number of tags, c is the chromosome,
                Nc is the number of tags on c, tags(c, d, s) is the boolean vector
                of tags on chromosome c and strand s, shifted by d, P[] is the Pearson
                correlation coefficient.

                Pearson correlation formula
                https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Mathematical_properties

                Linear algorithm to compute ∑xi*yi, ∑xi, ∑xi^2, ∑yi, ∑yi^2
                for Pearson coefficient computation
                In case of boolean unique alignment we deal with boolean vectors, where:
                    matchedTags = ∑xi*yi
                    positiveTags = ∑xi^2 = ∑xi
                    negativeTags = ∑yi^2 = ∑yi

                    Nc = ∑xi + ∑yi
                    P[xi, yi] = (Lc∑xi*yi - ∑xi*∑yi) / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi))
                where Lc is the length of c.

                The only value that depends on d is ∑xi*yi, so we can simplify the arg max.
                    arg max 1/N ∑Nc P[tags(c, d, +), tags(c, 0, -)] =
                    arg max ∑Nc P[tags(c, d, +), tags(c, 0, -)] =
                    arg max ∑(∑xi + ∑yi) (Lc∑xi*yi - ∑xi*∑yi) / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi)) =
                    arg max ∑(∑xi + ∑yi) Lc ∑xi*yi / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi)) =
                    arg max ∑ ∑xi*yi * (∑xi + ∑yi) Lc / √(∑xi * (Lc-∑xi) * ∑yi * (Lc-∑yi))
                We denote the value under arg max as "Pearson correlation transform".
            */
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
                    t1 == t2 -> {
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

        /**
         * Imputes the fragment size using the approach from Kharchenko et al., 2008
         * (see link above).
         * We ignore candidate fragment sizes less than [averageReadLength], following
         * the advice of Ramachandran et al., 2013:
         *      https://academic.oup.com/bioinformatics/article/29/4/444/200320
         * Marked internal for testing.
         */
        internal fun detectFragmentSize(
                data: GenomeStrandMap<TIntList>,
                averageReadLength: Double
        ): Int {
            if (averageReadLength.isNaN()) {
                // empty data, return a placeholder value
                return 0
            }
            // Ignore phantom peaks <= read length
            val candidateRange = averageReadLength.toInt()..MAX_FRAGMENT_SIZE
            if (candidateRange.isEmpty()) {
                return averageReadLength.toInt()
            }
            val ccs = computePearsonCorrelationTransform(candidateRange.toList(), data)
            return ccs.maxBy { it.pearsonTransform }!!.fragment
        }
    }
}
