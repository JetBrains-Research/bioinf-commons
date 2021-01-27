package org.jetbrains.bio.genome.containers

import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import org.jetbrains.bio.genome.Range


/**
 * A container for possibly overlapping ranges. Which doesn't merge overlapping ranges.
 *
 * // TODO: optional 'uniq' mode
 *
 * @see [org.jetbrains.bio.genome.Range] for details.
 * @author Roman Chernyatchik
 */
class RangesSortedList internal constructor(
    startOffsets: TIntList,
    endOffsets: TIntList
) : BaseRangesList(startOffsets, endOffsets) {

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
    override fun overlapRanges(startOffset: Int, endOffset: Int): Boolean {
        val rangesNumber = size
        for (idx in 0 until rangesNumber) {
            if (startOffsets[idx] >= endOffset) {
                // cannot overlap here and further, start offsets are sorted
                return false
            }

            if (endOffsets[idx] > startOffset) {
                return true
            }

        }
        return false
    }

    override fun includesRange(startOffset: Int, endOffset: Int): Boolean {
        val rangesNumber = size
        for (idx in 0 until rangesNumber) {
            if (startOffset < startOffsets[idx]) {
                // cannot contain here and further, start offsets are sorted
                return false
            }

            if (endOffset <= endOffsets[idx]) {
                return true
            }

        }
        return false
    }

    override fun intersectRanges(startOffset: Int, endOffset: Int): List<Range> {
        val rangesNumber = size

        val result = arrayListOf<Range>()

        for (idx in 0 until rangesNumber) {
            if (startOffsets[idx] >= endOffset) {
                // cannot intersect here and further, start offsets are sorted
                break
            }

            if (endOffsets[idx] > startOffset) {
                val start = kotlin.math.max(startOffset, startOffsets[idx])
                val end = kotlin.math.min(endOffset, endOffsets[idx])
                if(start < end) {
                    result.add(Range(start, end))
                }
            }
        }

        return result;
    }
}

fun Iterable<Range>.toRangeSortedList(): RangesSortedList {
    // sort by start,end
    val copy = sorted()

    // store ranges as: start, end offsets
    val startOffsets = TIntArrayList(copy.size)
    val endOffsets = TIntArrayList(copy.size)

    for (range in copy) {
        startOffsets.add(range.startOffset)
        endOffsets.add(range.endOffset)
    }
    return RangesSortedList(startOffsets, endOffsets)
}

fun rangeSortedList(vararg ranges: Range) = ranges.asIterable().toRangeSortedList()