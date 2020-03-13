package org.jetbrains.bio.genome.coverage

import com.google.common.base.MoreObjects
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import gnu.trove.set.hash.TIntHashSet
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.npy.NpzFile
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
        val averageInsertSize: Int,
        internal val data: GenomeMap<TIntList>
): Coverage {

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
    internal fun getTags(chromosomeRange: ChromosomeRange): IntArray {
        val data = data[chromosomeRange.chromosome]
        val index = data.binarySearchLeft(chromosomeRange.startOffset)
        var size = 0
        while (index + size < data.size() &&
                data[index + size] < chromosomeRange.endOffset) {
            size++
        }

        return data.toArray(index, size)
    }

    override val depth = genomeQuery.get().map { chr -> data[chr].size().toLong() }.sum()

    @Throws(IOException::class)
    internal fun save(outputPath: Path) {
        NpzFile.write(outputPath).use { writer ->
            writer.write(Coverage.VERSION_FIELD, intArrayOf(Coverage.VERSION))
            writer.write(Coverage.PAIRED_FIELD, booleanArrayOf(true))
            writer.write(PAIRED_VERSION_FIELD, intArrayOf(PAIRED_VERSION))
            writer.write(AVERAGE_INSERT_SIZE_FIELD, intArrayOf(averageInsertSize))

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
        private var totalInsertLength = 0L

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
            val insertSize = pos + length - pnext
            data[chromosome].add(pnext + insertSize / 2)
            readPairsCount++
            totalInsertLength += insertSize
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

            val averageInsertSize = if (readPairsCount != 0L) {
                (totalInsertLength / readPairsCount).toInt()
            } else {
                0
            }
            return PairedEndCoverage(
                    genomeQuery,
                    averageInsertSize = averageInsertSize,
                    data = data
            )
        }
    }

    companion object {

        const val PAIRED_VERSION = 2
        const val PAIRED_VERSION_FIELD = "paired_version"
        const val AVERAGE_INSERT_SIZE_FIELD = "average_insert_size"

        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        internal fun load(
                npzReader: NpzFile.Reader,
                genomeQuery: GenomeQuery
        ): PairedEndCoverage {
            check(npzReader[Coverage.PAIRED_FIELD].asBooleanArray().single()) {
                "${npzReader.path} attempting to read paired-end coverage from single-end cache file"
            }
            try {
                val version = npzReader[PAIRED_VERSION_FIELD].asIntArray().single()
                check(version == PAIRED_VERSION) {
                    "${npzReader.path} paired-end coverage version is $version " +
                            "instead of $PAIRED_VERSION"
                }
            } catch (e: IllegalStateException) {
                throw IllegalStateException(
                        "${npzReader.path} paired-end coverage version is missing",
                        e
                )
            }
            val averageInsertSize = npzReader[AVERAGE_INSERT_SIZE_FIELD].asIntArray().single()
            val data: GenomeMap<TIntList> = genomeMap(genomeQuery) { TIntArrayList() }
            for (chromosome in genomeQuery.get()) {
                try {
                    val npyArray = npzReader[chromosome.name]
                    data[chromosome] = TIntArrayList.wrap(npyArray.asIntArray())
                } catch (e: IllegalStateException) {
                    throw IllegalStateException(
                            "Cache file ${npzReader.path} doesn't contain ${chromosome.name}.\n" +
                                    "It's likely that chrom.sizes file used for its creation differs " +
                                    "from the one being used to read it now.\n" +
                                    "If problem persists, delete the cache file ${npzReader.path} " +
                                    "and Span will recreate it with correct settings.",
                            e
                    )
                }
            }
            return PairedEndCoverage(genomeQuery, averageInsertSize, data = data)
        }
    }

}
