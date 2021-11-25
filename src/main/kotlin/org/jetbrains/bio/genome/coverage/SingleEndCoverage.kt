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
import org.jetbrains.bio.genome.coverage.FragmentSize.detectFragmentSize
import org.jetbrains.bio.npy.NpzFile
import java.io.IOException
import java.nio.file.Path

/**
 * The container maintains a sorted list of tag offsets for each chromosome and strand.
 * Immutable. Saves data in [NpzFile] format.
 *
 * [detectedFragment] is estimated from the data using cross-correlation and saved
 * in the cache file.
 * [actualFragment] is used when computing location coverage and can be altered
 * by the user using [withFragment] method. By default, it's equal to [detectedFragment].
 *
 * @author Oleg Shpynov
 * @author Aleksei Dievskii
 */
class SingleEndCoverage private constructor(
    override val genomeQuery: GenomeQuery,
    val detectedFragment: Int,
    val actualFragment: Int = detectedFragment,
    internal val data: GenomeStrandMap<TIntList>
) : Coverage {

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
            data[index + size] < endOffset
        ) {
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


    class Builder(val genomeQuery: GenomeQuery) {

        val data: GenomeStrandMap<TIntList> = genomeStrandMap(genomeQuery) { _, _ ->
            TIntArrayList()
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

    }
}
