package org.jetbrains.bio.methylome

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.argSort
import org.jetbrains.bio.dataframe.div
import org.jetbrains.bio.dataframe.reorder
import org.jetbrains.bio.npy.NpzFile
import java.util.*

/**
 * A frame hold WGBS counts for a specific chromosome and strand.
 *
 * @author Sergei Lebedev
 */
internal class MethylomeFrame {
    private var offsets = IntArray(0)
    private var cytosineContextTags = ByteArray(0)
    private var methylatedCounts = ShortArray(0)
    private var totalCounts = ShortArray(0)
    var size = 0
        private set

    fun add(offset: Int, contextTag: Byte, methylatedCount: Short, totalCount: Short) {
        ensureCapacity(size + 1)

        offsets[size] = offset
        cytosineContextTags[size] = contextTag
        methylatedCounts[size] = methylatedCount
        totalCounts[size] = totalCount
        size++
    }

    fun sort() {
        trimToSize()

        val indices = offsets.argSort()
        offsets.reorder(indices)
        cytosineContextTags.reorder(indices)
        methylatedCounts.reorder(indices)
        totalCounts.reorder(indices)
    }

    /**
     * Loads a methylome frame in _sorted_ order from a given HDF5 file.
     */
    fun load(key: String, reader: NpzFile.Reader) {
        offsets = reader[key + "/offsets"].asIntArray()
        cytosineContextTags = reader[key + "/tags"].asByteArray()
        methylatedCounts = reader[key + "/methylated"].asShortArray()
        totalCounts = reader[key + "/total"].asShortArray()
        size = reader[key + "/size"].asIntArray().single()
    }

    /**
     * Saves a methylome frame in _sorted_ order into a given HDF5 file.
     */
    fun save(key: String, writer: NpzFile.Writer) {
        trimToSize()

        with(writer) {
            write("$key/offsets", offsets)
            write("$key/tags", cytosineContextTags)
            write("$key/methylated", methylatedCounts)
            write("$key/total", totalCounts)
            write("$key/size", intArrayOf(size))
        }
    }

    /**
     * Every methylome is also a data frame. Peel it and see for yourself.
     */
    internal fun peel(): DataFrame = DataFrame()
            .with("offset", offsets)
            .with("tag", cytosineContextTags)
            .with("level", methylatedCounts / totalCounts)
            .with("k", methylatedCounts)
            .with("n", totalCounts)

    private fun ensureCapacity(capacity: Int) {
        if (capacity > offsets.size) {
            val newCapacity = capacity + (capacity shr 1)  // 1.5x
            offsets = offsets.copyOf(newCapacity)
            cytosineContextTags = cytosineContextTags.copyOf(newCapacity)
            methylatedCounts = methylatedCounts.copyOf(newCapacity)
            totalCounts = totalCounts.copyOf(newCapacity)
        }
    }

    private fun trimToSize() {
        if (size < offsets.size) {
            offsets = offsets.copyOf(size)
            cytosineContextTags = cytosineContextTags.copyOf(size)
            methylatedCounts = methylatedCounts.copyOf(size)
            totalCounts = totalCounts.copyOf(size)
        }
    }

    override fun equals(other: Any?): Boolean = when {
        this === other -> true
        other == null || other !is MethylomeFrame -> false
        else -> {
            Arrays.equals(offsets, other.offsets) &&
            Arrays.equals(cytosineContextTags, other.cytosineContextTags) &&
            Arrays.equals(methylatedCounts, other.methylatedCounts) &&
            Arrays.equals(totalCounts, other.totalCounts)
        }
    }

    override fun hashCode() = Objects.hash(
            Arrays.hashCode(offsets),
            Arrays.hashCode(cytosineContextTags),
            Arrays.hashCode(methylatedCounts),
            Arrays.hashCode(totalCounts))
}
