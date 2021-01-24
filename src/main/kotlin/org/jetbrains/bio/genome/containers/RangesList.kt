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

    fun overlap(startOffset: Int, endOffset: Int): Boolean

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
    fun overlap(other: RangesList, flankBothSides: Int = 0): List<Range>

    /**
     * If overlap ranges not needed use this function to do less GS
     */
    fun overlapNumber(other: RangesList, flankBothSides: Int = 0): Int

    // TODO: rename includes?
    fun contains(startOffset: Int, endOffset: Int): Boolean
    operator fun contains(range: Range) = contains(range.startOffset, range.endOffset)
    operator fun contains(offset: Int) = contains(offset, offset + 1)

    /**
      * Intersect each list interval with requested [startOffset], [endOffset] range,
      *  empty intervals not reported
      */
    fun intersect(startOffset: Int, endOffset: Int): List<Range>
    /**
      * Intersect each list interval with requested [range], empty intervals not reported
      */
    fun intersect(range: Range) = intersect(range.startOffset, range.endOffset)
    infix fun and(other: RangesList): List<Range>
}

abstract class BaseRangesList(
    protected val startOffsets: TIntList,
    protected val endOffsets: TIntList
) : RangesList {

    init {
        require(startOffsets.size() == endOffsets.size())
    }

    override fun overlap(other: RangesList, flankBothSides: Int): List<Range> {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        val acc = ArrayList<Range>()

        (0 until size).forEach { idx ->
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]

            val startOffsetFlnk = kotlin.math.max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides
            if (other.overlap(startOffsetFlnk, endOffsetFlnk)) {
                acc.add(Range(startOffset, endOffset))
            }

        }

        return acc
    }

    override infix fun and(other: RangesList): List<Range> {
        val acc = ArrayList<Range>()

        (0 until size).forEach { idx ->
            acc.addAll(other.intersect(startOffsets[idx], endOffsets[idx]))
        }

        return acc
    }

    override fun overlapNumber(other: RangesList, flankBothSides: Int): Int {
        require(flankBothSides >= 0) { "Expected to be non-negative, but was $flankBothSides" }

        var acc = 0

        (0 until size).forEach { idx ->
            val startOffset = startOffsets[idx]
            val endOffset = endOffsets[idx]

            val startOffsetFlnk = kotlin.math.max(0, startOffset - flankBothSides)
            val endOffsetFlnk = endOffset + flankBothSides
            if (other.overlap(startOffsetFlnk, endOffsetFlnk)) {
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
        (0 until size).forEach { i ->
            it.printRecord(startOffsets[i], endOffsets[i])
        }
    }
}