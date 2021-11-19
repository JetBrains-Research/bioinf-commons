package org.jetbrains.bio.genome.containers

import com.google.common.collect.UnmodifiableIterator
import gnu.trove.list.TIntList
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.util.bufferedWriter
import java.io.IOException
import java.nio.file.Path

interface RangesList : Iterable<Range> {
    val size: Int

    fun overlapRanges(startOffset: Int, endOffset: Int): Boolean

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
    fun overlapRanges(other: RangesList, flankBothSides: Int = 0): List<Range>

    /**
     * If overlap ranges not needed use this function to do less GS
     */
    fun overlapRangesNumber(other: RangesList, flankBothSides: Int = 0): Int

    fun includesRange(startOffset: Int, endOffset: Int): Boolean
    fun includesRange(range: Range) = includesRange(range.startOffset, range.endOffset)
    fun includesRange(offset: Int) = includesRange(offset, offset + 1)

    /**
     * Intersect each list interval with requested [startOffset], [endOffset] range,
     *  empty intervals not reported
     */
    fun intersectRanges(startOffset: Int, endOffset: Int): List<Range>

    /**
     * Intersect each list interval with requested [range], empty intervals not reported
     */
    fun intersectRanges(range: Range) = intersectRanges(range.startOffset, range.endOffset)
    infix fun intersectRanges(other: RangesList): List<Range>

    /**
     * If intersected ranges not needed use this function to do less GS
     */
    fun intersectRangesNumber(other: RangesList, flankBothSides: Int = 0): Int
}

abstract class BaseRangesList(
    protected val startOffsets: TIntList,
    protected val endOffsets: TIntList
) : RangesList {

    init {
        require(startOffsets.size() == endOffsets.size())
    }

    override fun overlapRanges(other: RangesList, flankBothSides: Int): List<Range> {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        val acc = ArrayList<Range>()

        for (idx in 0 until size) {
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]

            val startOffsetFlnk = kotlin.math.max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides
            if (other.overlapRanges(startOffsetFlnk, endOffsetFlnk)) {
                acc.add(Range(startOffset, endOffset))
            }
        }

        return acc
    }

    override infix fun intersectRanges(other: RangesList): List<Range> {
        val acc = ArrayList<Range>()

        for (idx in 0 until size) {
            acc.addAll(other.intersectRanges(startOffsets[idx], endOffsets[idx]))
        }

        return acc
    }

    override fun intersectRangesNumber(other: RangesList, flankBothSides: Int): Int {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        var acc = 0

        for (idx in 0 until size) {
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]

            val startOffsetFlnk = kotlin.math.max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides
            acc += other.intersectRanges(startOffsetFlnk, endOffsetFlnk).size
        }

        return acc
    }

    override fun overlapRangesNumber(other: RangesList, flankBothSides: Int): Int {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        var acc = 0

        for (idx in 0 until size) {
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]

            val startOffsetFlnk = kotlin.math.max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides
            if (other.overlapRanges(startOffsetFlnk, endOffsetFlnk)) {
                acc++
            }
        }

        return acc
    }

    /** The number of ranges in this list. */
    override val size: Int get() = startOffsets.size()

    override fun iterator() = object : UnmodifiableIterator<Range>() {
        private var current = 0

        override fun hasNext() = current < size

        override fun next(): Range {
            val range = Range(startOffsets[current], endOffsets[current])
            current++
            return range
        }
    }

    override fun toString() = "[ranges=${joinToString(", ")}]"
    fun dump() = "${toString()}\nstarts=$startOffsets\nends=$endOffsets]"

    @Throws(IOException::class)
    fun save(path: Path) = CSVFormat.TDF.print(path.bufferedWriter()).use {
        for (i in 0 until size) {
            it.printRecord(startOffsets[i], endOffsets[i])
        }
    }
}