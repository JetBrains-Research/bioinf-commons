package org.jetbrains.bio.genome.containers

import com.google.common.collect.Iterables
import com.google.common.collect.UnmodifiableIterator
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import java.io.IOException
import java.nio.file.Path
import java.util.*

/**
 * A container for possibly overlapping ranges.
 *
 * @see [org.jetbrains.bio.genome.Range] for details.
 * @author Sergei Lebedev
 */
class RangesMergingList internal constructor(
        private val startOffsets: TIntList,
        private val endOffsets: TIntList) : Iterable<Range> {

    init {
        require(startOffsets.size() == endOffsets.size())
    }

    /**
     * Performs element-wise union on ranges in the two lists.
     */
    infix fun or(other: RangesMergingList) = Iterables.concat(this, other).toRangeList()

    /**
     * Performs element-wise intersection on ranges in the two lists.
     */
    infix fun and(other: RangesMergingList): RangesMergingList {
        val acc = ArrayList<Range>()
        for ((startOffset, endOffset) in other) {
            var i = Math.max(0, lookup(startOffset))
            while (i < size &&  // vvv inlined Range#intersects.
                   endOffset > startOffsets[i]) {
                if (endOffsets[i] > startOffset) {
                    val intersection = Range(
                            Math.max(startOffsets[i], startOffset),
                            Math.min(endOffsets[i], endOffset))
                    acc.add(intersection)
                }

                i++
            }
        }

        return acc.toRangeList()
    }

    /**
     * Leaves only ranges intersecting some range form other
     *
     * THIS            |----------|
     * ****************************************
     * OTHER 1      |----|    |----|        : +
     * OTHER 2     |-----------------|      : +
     * OTHER 3          |--|                : +
     * OTHER 4  |--|                        : -
     * OTHER 5                       |--|   : -
     * OTHER 6  |-|     |----|              : +
     */
    infix fun overlap(other: RangesMergingList): RangesMergingList {
        val acc = ArrayList<Range>()

        for ((startOffset, endOffset) in this) {
            var i = Math.max(0, other.lookup(startOffset))
            while (i < other.size && other.startOffsets[i] < endOffset) {
                if (other.endOffsets[i] > startOffset) {
                    // doesn't intersect
                    acc.add(Range(startOffset, endOffset))
                    break
                }
                i++
            }
        }
        return acc.toRangeList()
    }

    operator fun contains(range: Range) = contains(range.startOffset, range.endOffset)
    operator fun contains(offset: Int) = contains(offset, offset + 1)

    fun contains(startOffset: Int, endOffset: Int): Boolean {
        val i = lookup(startOffset)
        return i >= 0 && startOffset >= startOffsets[i] && endOffset <= endOffsets[i]
    }

    fun intersectionLength(range: Range): Int {
        return intersect(range).map(Range::length).sum()
    }

    fun intersect(range: Range): List<Range> {
        var i = Math.max(0, lookup(range.startOffset))
        val result = arrayListOf<Range>()
        while (i < size) {
            // Iterate over nearby ranges.
            if (startOffsets[i] >= range.endOffset) {
                break
            }

            val start = Math.max(range.startOffset, startOffsets[i])
            val end = Math.min(range.endOffset, endOffsets[i])
            if (start < end) {
                result.add(Range(start, Math.max(0, end)))
            }
            i++
        }
        return result
    }

    /**
     * Returns the index of the first range starting before the given
     * [offset].
     */
    private fun lookup(offset: Int): Int {
        val i = startOffsets.binarySearch(offset)
        return if (i < 0) i.inv() - 1 else i
    }

    /** The number of ranges in this list. */
    val size: Int get() = startOffsets.size()


    override fun toString() = "[${joinToString(", ")}]"

    @Throws(IOException::class)
    fun save(path: Path) = CSVFormat.TDF.print(path.bufferedWriter()).use {
        for ((startOffset, endOffset) in this) {
            it.printRecord(startOffset, endOffset)
        }
    }

    override fun iterator() = object : UnmodifiableIterator<Range>() {
        private var current = 0

        override fun hasNext() = current < size

        override fun next(): Range {
            val range = Range(startOffsets[current], endOffsets[current])
            current++
            return range
        }
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is RangesMergingList -> false
        else -> startOffsets == other.startOffsets &&
                endOffsets == other.endOffsets
    }

    override fun hashCode() = Objects.hash(startOffsets, endOffsets)

    companion object {
        @Throws(IOException::class)
        fun load(path: Path) = CSVFormat.TDF.parse(path.bufferedReader()).use { csvParser ->
            csvParser.map { Range(it[0].toInt(), it[1].toInt()) }.toRangeList()
        }
    }
}

/**
 * Constructs a list of complementary ranges.
 *
 *   |---------------------|  this
 *     |--|  |-----|     |-|  ranges
 *
 *   |-|  |--|     |-----|    result
 *
 * Input ranges may be in any order and are allowed to overlap.
 */
operator fun Range.minus(ranges: List<Range>): List<Range> {
    if (ranges.isEmpty()) {
        return listOf(this)
    }

    val acc = ArrayList<Range>()
    var current = startOffset
    for ((startOffset, endOffset) in ranges.toRangeList()) {
        if (startOffset > current) {
            acc.add(Range(current, startOffset))
        }

        current = endOffset
    }

    if (current < endOffset) {
        acc.add(Range(current, endOffset))
    }

    return acc
}

fun rangeList(vararg ranges: Range) = ranges.asIterable().toRangeList()

fun Sequence<Range>.toRangeList() = asIterable().toRangeList()

fun Iterable<Range>.toRangeList(): RangesMergingList {
    val copy = sorted()

    // Try to compress overlapping ranges.
    val startOffsets = TIntArrayList()
    val endOffsets = TIntArrayList()
    var end = 0
    for (range in copy) {
        val start = range.startOffset
        if (startOffsets.isEmpty) {
            startOffsets.add(start)
        } else if (start > end) {
            endOffsets.add(end)
            startOffsets.add(start)
        }

        end = Math.max(end, range.endOffset)
    }

    if (!startOffsets.isEmpty) {
        endOffsets.add(end)
    }

    return RangesMergingList(startOffsets, endOffsets)
}