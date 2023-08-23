package org.jetbrains.bio.genome.containers

import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import org.jetbrains.annotations.TestOnly
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Range
import kotlin.math.max
import kotlin.math.min


/**
 * A container for possibly overlapping ranges. Which doesn't merge overlapping ranges.
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
        var i = max(0, lookupLeft(startOffset))
        // check if this 'idx' range intersects at least one other region:
        while (i < rangesNumber && startOffsets[i] < endOffset) {
            if (endOffsets[i] > startOffset) {
                return true
            } // else doesn't intersect
            i++
        }
        return false
    }

    /**
     * Returns the index of the insertion point. So the [startOffsets] at the index is >= [startOffset] and
     * [startOffsets] at the index-1 is < [startOffset]. The index could be `-1`.
     */
    private fun lookupLeft(startOffset: Int): Int {
        val i = startOffsets.binarySearch(startOffset)
        var idx = (if (i < 0) i.inv() - 1 else i)

        // optional scroll left, to make binary search more consistent:
        while (idx > 0 && (startOffsets[idx] == startOffsets[idx-1])) {
            idx-=1;
        }

        return idx
    }
    @TestOnly
    fun internalLookupLeft(startOffset: Int) = lookupLeft(startOffset)

    override fun includesRange(startOffset: Int, endOffset: Int): Boolean {
        // XXX: cannot use binary search here due to intersected ranges (see tests)

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
        // XXX: cannot use binary search here due to intersected ranges (see tests)

        val rangesNumber = size

        val result = arrayListOf<Range>()

        for (idx in 0 until rangesNumber) {
            if (startOffsets[idx] >= endOffset) {
                // cannot intersect here and further, start offsets are sorted
                break
            }

            if (endOffsets[idx] > startOffset) {
                val start = max(startOffset, startOffsets[idx])
                val end = min(endOffset, endOffsets[idx])
                if (start < end) {
                    result.add(Range(start, end))
                }
            }
        }

        return result
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

fun Location.toRangeSortedList(): RangesSortedList {
    // store ranges as: start, end offsets
    val startOffsets = TIntArrayList(1)
    val endOffsets = TIntArrayList(1)

    startOffsets.add(startOffset)
    endOffsets.add(endOffset)

    return RangesSortedList(startOffsets, endOffsets)
}

fun rangeSortedList(vararg ranges: Range) = ranges.asIterable().toRangeSortedList()