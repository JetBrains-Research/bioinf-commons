package org.jetbrains.bio.genome.containers

import com.google.common.collect.Iterables
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import org.apache.commons.csv.CSVFormat
import org.jetbrains.annotations.TestOnly
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.util.bufferedReader
import java.io.IOException
import java.nio.file.Path
import java.util.*
import kotlin.math.max
import kotlin.math.min

/**
 * A container for possibly overlapping ranges which merges overlapping ranges.
 *
 * @see [org.jetbrains.bio.genome.Range] for details.
 * @author Sergei Lebedev
 */
class RangesMergingList internal constructor(
    startOffsets: TIntList,
    endOffsets: TIntList
) : BaseRangesList(startOffsets, endOffsets) {

    /**
     * Performs element-wise union on ranges in the two lists.
     */
    infix fun or(other: RangesMergingList) = Iterables.concat(this, other).toRangeMergingList()

    override fun overlapRanges(startOffset: Int, endOffset: Int): Boolean {
        val rangesNumber = size
        var i = max(0, lookup(startOffset))
        // check if this 'idx' range intersects at least one other region:
        while (i < rangesNumber && startOffsets[i] < endOffset) {
            if (endOffsets[i] > startOffset) {
                return true
            } // else doesn't intersect
            i++
        }
        return false
    }

    override fun includesRange(startOffset: Int, endOffset: Int): Boolean {
        val i = lookup(startOffset)
        return i >= 0 && startOffset >= startOffsets[i] && endOffset <= endOffsets[i]
    }

    fun intersectionLength(range: Range) = intersectRanges(range).map(Range::length).sum()

    override fun intersectRanges(startOffset: Int, endOffset: Int): List<Range> {
        var i = max(0, lookup(startOffset))
        val result = arrayListOf<Range>()
        val len = size
        while (i < len) {
            // Iterate over nearby ranges.
            if (startOffsets[i] >= endOffset) {
                break  // break if start offset became out of range
            }

            val start = max(startOffset, startOffsets[i])
            val end = min(endOffset, endOffsets[i])
            if (start < end) {
                // add only valid ranges, invalid range could appear e.g.:
                // list = [10, 30] [50, 100] , rage=[45, 60]
                // candidates: [45, 30], [45, 60]
                result.add(Range(start, max(0, end)))
            }
            i++
        }
        return result
    }

    /**
     * Constructs a list of complementary ranges.
     *
     *     |----------| |---------|  |---|  this
     *  |--|          |-|         |--|      result
     *
     */
    fun complementaryRanges(chrLen: Int): RangesMergingList {
        val len = size

        val complStarts = TIntArrayList(len + 1)
        val complEnds = TIntArrayList(len + 1)

        if (len == 0) {
            complStarts.add(0)
            complEnds.add(chrLen)
        } else {
            var i = 0
            while (i < len) {
                val start = startOffsets[i]
                val end = endOffsets[i]

                if (i == 0) {
                    if (start > 0) {
                        complStarts.add(0)
                        complEnds.add(start)
                    }
                } else {
                    complEnds.add(start)
                }
                complStarts.add(end)

                i++
            }

            if (endOffsets[len - 1] != chrLen) {
                complEnds.add(chrLen)
            } else {
                complStarts.removeAt(complStarts.size() - 1)
            }
        }
        return RangesMergingList(complStarts, complEnds)
    }

    /**
     * Returns the index of the insertion point. So the [startOffsets] at the index is >= [startOffset] and
     * [startOffsets] at the index-1 is < [startOffset]. The index could be `-1`.
     */
    private fun lookup(startOffset: Int): Int {
        val i = startOffsets.binarySearch(startOffset)
        return if (i < 0) i.inv() - 1 else i
    }

    @TestOnly
    fun internalLookup(startOffset: Int) = lookup(startOffset)

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
            csvParser.map { Range(it[0].toInt(), it[1].toInt()) }.toRangeMergingList()
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
    for ((startOffset, endOffset) in ranges.toRangeMergingList()) {
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

fun rangeMergingList(vararg ranges: Range) = ranges.asIterable().toRangeMergingList()

fun Sequence<Range>.toRangeMergingList() = asIterable().toRangeMergingList()

fun Iterable<Range>.toRangeMergingList(): RangesMergingList {
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

        end = max(end, range.endOffset)
    }

    if (!startOffsets.isEmpty) {
        endOffsets.add(end)
    }

    // make compact:
    startOffsets.trimToSize()
    endOffsets.trimToSize()

    return RangesMergingList(startOffsets, endOffsets)
}