package org.jetbrains.bio.genome.coverage

import com.google.common.base.MoreObjects
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import gnu.trove.set.hash.TIntHashSet
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.npy.NpzFile
import org.slf4j.LoggerFactory
import java.io.IOException
import java.nio.file.Path

/**
 * The container stores paired-end coverage information.
 * Immutable. Saves data in [NpzFile] format.
 *
 * For paired-end coverage, we don't need any fragment size or shift information.
 * We can just reduce each read pair to a tag at the middle of the pair's 5' bounds
 * on the plus strand.
 * This coverage will always be perfectly strand-asymmetric, residing on plus strand entirely.
 */
class PairedEndCoverage private constructor(
    override val genomeQuery: GenomeQuery,
    val averageFragmentSize: Int,
    val data: GenomeMap<TIntList>
) : Coverage {

    /**
     * Returns the number of tags covered by a given [location].
     */
    override fun getCoverage(location: Location) = location.strand.choose(
        ifPlus = { getTags(location.toChromosomeRange()).size },
        ifMinus = { 0 }
    )

    override fun getBothStrandsCoverage(chromosomeRange: ChromosomeRange): Int =
        getTags(chromosomeRange).size

    /**
     * Returns a sorted array of tags covered by a given [chromosomeRange].
     */
    fun getTags(chromosomeRange: ChromosomeRange): IntArray {
        val data = data[chromosomeRange.chromosome]
        val index = data.binarySearchLeft(chromosomeRange.startOffset)
        var size = 0
        while (index + size < data.size() &&
            data[index + size] < chromosomeRange.endOffset
        ) {
            size++
        }

        return data.toArray(index, size)
    }

    override val depth = genomeQuery.get().map { chr -> data[chr].size().toLong() }.sum()

    @Throws(IOException::class)
    fun save(outputPath: Path) {
        NpzFile.write(outputPath).use { writer ->
            writer.write(Coverage.VERSION_FIELD, intArrayOf(Coverage.VERSION))
            writer.write(Coverage.PAIRED_FIELD, booleanArrayOf(true))
            writer.write(PAIRED_VERSION_FIELD, intArrayOf(PAIRED_VERSION))
            writer.write(AVERAGE_FRAGMENT_SIZE_FIELD, intArrayOf(averageFragmentSize))

            for (chromosome in genomeQuery.get()) {
                val key = chromosome.name
                writer.write(key, data[chromosome].toArray())
            }
        }
    }

    override fun toString() = MoreObjects.toStringHelper(this)
        .addValue(genomeQuery).toString()

    class Builder(
        val genomeQuery: GenomeQuery,
        val data: GenomeMap<TIntList> = genomeMap(genomeQuery) { TIntArrayList() }
    ) {

        private var readPairsCount = 0L
        private var totalFragmentSize = 0L

        /**
         * We expect this method to be called for only one read of each read pair: the one
         * on the negative strand. We thus expect that no same-strand pairs are provided.
         * We also expect only same-reference pairs.
         *
         * [chromosome] is the mapping chromosome.
         * [pos] and [pnext] are POS and PNEXT, the standard SAM fields.
         * [length] is the length of the read.
         */
        fun process(
            chromosome: Chromosome,
            pos: Int, pnext: Int, length: Int
        ): Builder {
            // Process case when mate pair begins before
            if (pos > pnext) {
                return process(chromosome, pnext, pos, length)
            }
            val fragmentSize = pnext + length - pos
            check(fragmentSize >= 0) {
                "Negative fragment size chromosome=$chromosome, pos=$pos, pnext=$pnext, length=$length"
            }
            if (fragmentSize > FragmentSize.MAX_FRAGMENT_SIZE) {
                return this  // Ignore too long fragment sizes
            }
            readPairsCount++
            data[chromosome].add(pos + fragmentSize / 2)
            totalFragmentSize += fragmentSize
            return this
        }

        /**
         * Generate a [PairedEndCoverage] object.
         * [unique] controls whether duplicate tags should be preserved ([unique] == false)
         * or squished into one tag ([unique] == true).
         * Only tags at the exact same offset are considered duplicate.
         */
        fun build(unique: Boolean): PairedEndCoverage {
            for (chromosome in data.genomeQuery.get()) {
                if (unique) {
                    // XXX we can do linear time de-duplication on
                    //     the sorted sequence.
                    val tags = data[chromosome]
                    data[chromosome] = TIntArrayList(TIntHashSet(tags))
                }
                data[chromosome].sort()
            }

            val averageFragmentSize = if (readPairsCount != 0L) {
                (totalFragmentSize / readPairsCount).toInt()
            } else {
                0
            }
            return PairedEndCoverage(
                genomeQuery,
                averageFragmentSize = averageFragmentSize,
                data = data
            )
        }
    }

    companion object {
        const val PAIRED_VERSION = 3
        const val PAIRED_VERSION_FIELD = "paired_version"
        const val AVERAGE_FRAGMENT_SIZE_FIELD = "average_fragment_size"

        private val LOG = LoggerFactory.getLogger(PairedEndCoverage::class.java)

        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        internal fun load(
            npzReader: NpzFile.Reader,
            path: Path,
            genomeQuery: GenomeQuery,
            failOnMissingChromosomes: Boolean
        ): PairedEndCoverage {
            check(npzReader[Coverage.PAIRED_FIELD].asBooleanArray().single()) {
                "$path attempting to read paired-end coverage from single-end cache file"
            }
            try {
                val version = npzReader[PAIRED_VERSION_FIELD].asIntArray().single()
                check(version == PAIRED_VERSION) {
                    "$path paired-end coverage version is $version instead of $PAIRED_VERSION"
                }
            } catch (e: IllegalStateException) {
                throw IllegalStateException("$path paired-end coverage version is missing", e)
            }
            val averageFragmentSize = npzReader[AVERAGE_FRAGMENT_SIZE_FIELD].asIntArray().single()
            val data: GenomeMap<TIntList> = genomeMap(genomeQuery) { TIntArrayList() }
            for (chromosome in genomeQuery.get()) {
                try {
                    val npyArray = npzReader[chromosome.name]
                    data[chromosome] = TIntArrayList.wrap(npyArray.asIntArray())
                } catch (e: NullPointerException) { // JDK11 doesn't throw ISE in case of missing zip entry
                    val msg = "File $path doesn't contain data for ${chromosome.name}."
                    LOG.trace(msg)
                    if (failOnMissingChromosomes) {
                        throw java.lang.IllegalStateException(msg, e)
                    }
                } catch (e: IllegalStateException) {
                    val msg = "File $path doesn't contain data for ${chromosome.name}."
                    LOG.trace(msg)
                    if (failOnMissingChromosomes) {
                        throw java.lang.IllegalStateException(msg, e)
                    }
                }
            }
            return PairedEndCoverage(genomeQuery, averageFragmentSize, data = data)
        }
    }

}
