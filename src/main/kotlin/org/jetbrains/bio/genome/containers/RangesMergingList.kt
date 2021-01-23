package org.jetbrains.bio.genome.containers

import com.google.common.collect.Iterables
import com.google.common.collect.UnmodifiableIterator
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import org.apache.commons.csv.CSVFormat
import org.jetbrains.annotations.TestOnly
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import java.io.IOException
import java.nio.file.Path
import java.util.*
import kotlin.math.max
import kotlin.math.min

/**
 * A container for possibly overlapping ranges.
 *
 * @see [org.jetbrains.bio.genome.Range] for details.
 * @author Sergei Lebedev
 */
class RangesMergingList internal constructor(
    private val startOffsets: TIntList,
    private val endOffsets: TIntList
) : Iterable<Range> {

    init {
        require(startOffsets.size() == endOffsets.size())
    }

    /**
     * Performs element-wise union on ranges in the two lists.
     */
    infix fun or(other: RangesMergingList) = Iterables.concat(this, other).toRangeMergingList()

    /**
     * Performs element-wise intersection on ranges in the two lists.
     */
    infix fun and(other: RangesMergingList): RangesMergingList {
        val acc = ArrayList<Range>()

        (0 until other.size).forEach { idx ->
            val startOffset = other.startOffsets[idx]
            val endOffset = other.endOffsets[idx]

            var i = max(0, lookup(startOffset))
            while (i < size &&  // vvv inlined Range#intersects.
                endOffset > startOffsets[i]
            ) {
                if (endOffsets[i] > startOffset) {
                    val intersection = Range(
                        max(startOffsets[i], startOffset),
                        min(endOffsets[i], endOffset)
                    )
                    acc.add(intersection)
                }

                i++
            }
        }

        return acc.toRangeMergingList()
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
     *
     * @param other Other list. Methods finds ranges from target list with intersects any range from [other] list
     * @param flankBothSides Flank each range from target (this) list with [flankBothSides] positions at start and end
     */
    fun overlap(other: RangesMergingList, flankBothSides: Int = 0): RangesMergingList {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        val acc = ArrayList<Range>()

        (0 until size).forEach { idx ->
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]
            val startOffsetFlnk = max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides

            var i = max(0, other.lookup(startOffsetFlnk))
            // check if this 'idx' range intersects at least one other region:
            while (i < other.size && other.startOffsets[i] < endOffsetFlnk) {
                if (other.endOffsets[i] > startOffsetFlnk) {
                    acc.add(Range(startOffset, endOffset))
                    break
                } // else doesn't intersect
                i++
            }
        }

        return acc.toRangeMergingList()
    }

    private fun overlap(idx: Int, other: RangesMergingList, flankBothSides: Int): Range? {

        return null
    }

    operator fun contains(range: Range) = contains(range.startOffset, range.endOffset)
    operator fun contains(offset: Int) = contains(offset, offset + 1)

    fun contains(startOffset: Int, endOffset: Int): Boolean {
        val i = lookup(startOffset)
        return i >= 0 && startOffset >= startOffsets[i] && endOffset <= endOffsets[i]
    }

    fun intersectionLength(range: Range) = intersect(range).map(Range::length).sum()

    /**
     * Intersect each list interval with requested [range], empty intervals not reported
     */
    fun intersect(range: Range): List<Range> {
        var i = max(0, lookup(range.startOffset))
        val result = arrayListOf<Range>()
        val len = size;
        while (i < len) {
            // Iterate over nearby ranges.
            if (startOffsets[i] >= range.endOffset) {
                break  // break if start offset became out of range
            }

            val start = max(range.startOffset, startOffsets[i])
            val end = min(range.endOffset, endOffsets[i])
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
     * Returns the index of the insertion point. So the [startOffsets] at the index is >= [startOffset] and
     * [startOffsets] at the index-1 is < [startOffset]. The index could be `-1`.
     */
    private fun lookup(startOffset: Int): Int {
        val i = startOffsets.binarySearch(startOffset)
        return if (i < 0) i.inv() - 1 else i
    }

    @TestOnly
    fun internalLookup(startOffset: Int) = lookup(startOffset)

    /** The number of ranges in this list. */
    val size: Int get() = startOffsets.size()


    override fun toString() = "[ranges=${joinToString(", ")}]"
    fun dump() = "${toString()}\nstarts=$startOffsets\nends=$endOffsets"

    @Throws(IOException::class)
    fun save(path: Path) = CSVFormat.TDF.print(path.bufferedWriter()).use {
        (0 until size).forEach { i ->
            it.printRecord(startOffsets[i], endOffsets[i])
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

    return RangesMergingList(startOffsets, endOffsets)
}