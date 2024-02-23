package org.jetbrains.bio.genome.containers

import com.google.common.collect.Iterables
import org.jetbrains.bio.genome.Range
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class RangesMergingListTest {
    @Test
    fun empty() {
        assertFalse(rangeMergingList().includesRange(Range(0, 12)))
        assertEquals(0, rangeMergingList().intersectionLength(Range(0, 12)))
    }

    @Test
    fun single() {
        assertTrue(rangeMergingList(Range(0, 12)).includesRange(Range(0, 12)))
        assertFalse(rangeMergingList(Range(0, 12)).includesRange(Range(0, 20)))
        assertEquals(12, rangeMergingList(Range(0, 12)).intersectionLength(Range(0, 12)))
    }

    @Test
    fun singleBigSegment() {
        assertTrue(rangeMergingList(Range(0, 100)).includesRange(Range(10, 20)))
    }

    @Test
    fun singleExceedingRange() {
        assertFalse(rangeMergingList(Range(0, 10)).includesRange(Range(-100, 100)))
    }

    @Test
    fun intersectionLength() {
        val rl = rangeMergingList(Range(0, 10))
        assertFalse(rl.includesRange(Range(5, 15)))
        assertEquals(5, rl.intersectionLength(Range(5, 15)))
        assertEquals(10, rl.intersectionLength(Range(0, 15)))
        assertEquals(5, rl.intersectionLength(Range(5, 10)))
        assertEquals(5, rl.intersectionLength(Range(-10, 5)))
        assertEquals(10, rl.intersectionLength(Range(-10, 10)))
        assertEquals(0, rl.intersectionLength(Range(20, 30)))
        assertEquals(0, rl.intersectionLength(Range(10, 15)))
    }

    @Test
    fun doubleIntersectionLength() {
        assertTrue(rangeMergingList(Range(0, 10), Range(10, 20)).includesRange(Range(5, 15)))
        assertFalse(rangeMergingList(Range(0, 10), Range(11, 20)).includesRange(Range(5, 15)))
        assertEquals(10, rangeMergingList(Range(0, 10), Range(10, 20)).intersectionLength(Range(5, 15)))
        assertEquals(5, rangeMergingList(Range(0, 10), Range(20, 30)).intersectionLength(Range(15, 25)))
        assertEquals(0, rangeMergingList(Range(0, 10), Range(10, 20)).intersectionLength(Range(-10, -5)))

        val rl = rangeMergingList(Range(0, 10), Range(20, 30), Range(40, 50))
        assertEquals(1, rl.intersectionLength(Range(35, 41)))
        assertEquals(0, rl.intersectionLength(Range(30, 40)))
        assertEquals(1, rl.intersectionLength(Range(29, 40)))
        assertEquals(1, rl.intersectionLength(Range(30, 41)))
        assertEquals(2, rl.intersectionLength(Range(29, 41)))
        assertEquals(20, rl.intersectionLength(Range(5, 45)))
    }

    @Test
    fun nonSorted() {
        assertTrue(rangeMergingList(Range(10, 20), Range(0, 5)).includesRange(Range(11, 12)))
    }

    @Test
    fun intersectingRanges() {
        assertEquals(30, rangeMergingList(Range(0, 20), Range(10, 30)).intersectionLength(Range(0, 30)))
    }

    @Test
    fun iteration() {
        assertEquals(0, Iterables.size(rangeMergingList()))

        val ranges = arrayOf(Range(0, 20), Range(25, 30))
        assertEquals(ranges.toList(), rangeMergingList(*ranges).toList())
    }

    @Test
    fun orWithEmptyOrSelf() {
        val rl = rangeMergingList(Range(0, 10), Range(20, 30))
        assertEquals(rl.toList(), (rl or rangeMergingList()).toList())
        assertEquals(rl.toList(), (rl or rl).toList())
    }

    @Test
    fun orWithoutIntersections() {
        val rl1 = rangeMergingList(Range(0, 10))
        val rl2 = rangeMergingList(Range(20, 30))
        assertEquals(
                rangeMergingList(Range(0, 10), Range(20, 30)).toList(),
                (rl1 or rl2).toList()
        )
    }

    @Test
    fun orWithIntersections() {
        assertEquals(
                rangeMergingList(Range(0, 30)).toList(),
                (rangeMergingList(Range(0, 10)) or rangeMergingList(Range(5, 30))).toList()
        )
    }

    @Test
    fun intersectRangesWithEmpty() {
        val rl = rangeMergingList(Range(0, 10), Range(20, 30))
        doCheckIntersectRanges(emptyList<Range>(), rl, rangeMergingList())
    }

    @Test
    fun intersectRangesWithSelf() {
        val rl = rangeMergingList(Range(0, 10), Range(20, 30))
        doCheckIntersectRanges(rl.toList(), rl, rl)
    }

    @Test
    fun intersectRanges() {
        val rl1 = rangeMergingList(Range(0, 10))
        val rl2 = rangeMergingList(Range(5, 30))
        doCheckIntersectRanges(rangeMergingList(Range(5, 10)).toList(), rl1, rl2)
        doCheckIntersectRanges(
                rangeMergingList(Range(5, 10), Range(25, 30)).toList(),
                rl1 or rangeMergingList(Range(25, 30)),
                rl2
        )
        doCheckIntersectRanges(
                rangeMergingList(Range(5, 10), Range(25, 30)).toList(),
                rl2,
                rl1 or rangeMergingList(Range(25, 30))
        )

        doCheckIntersectRanges(rangeMergingList(Range(6, 8)).toList(), rangeMergingList(Range(0, 2), Range(6, 8)), rl2)
    }

    @Test
    fun intersectRangesFlnk() {
        val rl1 = rangeMergingList(Range(0, 10), Range(20, 25))
        val rl2 = rangeMergingList(Range(0, 30))
        assertEquals(listOf(Range(0, 12), Range(18, 27)), (rl1.intersectRanges(rl2, 2)).toList())
        assertEquals(listOf(Range(11, 12), Range(18, 20)), (rl1.intersectRanges(rangeMergingList(Range(11, 20)), 2)))
    }

    @Test
    fun intersectRangesFlnkOverlapAfterFlanking() {
        val rl1 = rangeMergingList(Range(0, 11), Range(13, 25))
        val rl2 = rangeMergingList(Range(0, 30))
        assertEquals(listOf(Range(0, 14), Range(10, 28)), (rl1.intersectRanges(rl2, 3)).toList())
    }

    private fun doCheckIntersectRanges(
            expected: List<Range>,
            rl1: RangesMergingList,
            rl2: RangesMergingList
    ) {
        assertEquals(expected, (rl1.intersectRanges(rl2)).toList())
        assertEquals(expected.size, (rl1.intersectRangesNumber(rl2)))
    }

    @Test
    fun overlapRegion() {
        assertFalse(rangeMergingList().overlapRanges(7, 12))

        val rl1 = rangeMergingList(Range(10, 30))
        assertTrue(rl1.overlapRanges(7, 12))
        assertTrue(rl1.overlapRanges(20, 40))

        assertFalse(rl1.overlapRanges(2, 7))
        assertFalse(rl1.overlapRanges(7, 10))
        assertFalse(rl1.overlapRanges(30, 40))
        assertFalse(rl1.overlapRanges(40, 50))

        val rl2 = rangeMergingList(Range(0, 30))
        assertTrue(rl2.overlapRanges(0, 4))
        assertTrue(rl2.overlapRanges(29, 30))
    }

    @Test
    fun overlap() {
        val rl1 = rangeMergingList(Range(10, 30))

        checkOverlap(listOf(), rl1, rangeMergingList())
        checkOverlap(rl1.toList(), rl1, rangeMergingList(Range(7, 12)))
        checkOverlap(rl1.toList(), rl1, rangeMergingList(Range(2, 7), Range(20, 40)))
        checkOverlap(listOf(), rl1, rangeMergingList(Range(2, 7)))
    }


    @Test
    fun overlapFlnk() {
        val rl1 = rangeMergingList(Range(10, 30), Range(50, 100))

        checkOverlap(listOf(Range(10, 30)), rl1, rangeMergingList(Range(7, 12)), 10)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeMergingList(Range(35, 40)), 10)
        checkOverlap(rl1.toList(), rl1, rangeMergingList(Range(35, 40)), 11)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeMergingList(Range(2, 7), Range(20, 40)), 2)
        checkOverlap(listOf(), rl1, rangeMergingList(Range(2, 7)), 3)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeMergingList(Range(2, 7)), 4)
    }

    private fun checkOverlap(
            expected: List<Range>,
            a: RangesMergingList,
            b: RangesMergingList,
            flankBothSides: Int = 0
    ) {
        assertEquals(expected, (a.overlapRanges(b, flankBothSides)).toList())
        assertEquals(expected.size, a.overlapRangesNumber(b, flankBothSides))
    }


    @Test
    fun rangeMinus() {
        val universe = Range(0, 100)
        assertEquals(listOf(universe), universe - emptyList())
        assertEquals(
                listOf(Range(0, 20), Range(40, 100)),
                universe - listOf(Range(20, 40))
        )
        assertEquals(listOf(Range(40, 100)), universe - listOf(Range(0, 40)))
        assertEquals(listOf(Range(0, 40)), universe - listOf(Range(40, 100)))
        assertEquals(
                emptyList<Range>(),
                universe - listOf(Range(0, 40), Range(40, 100))
        )
    }

    @Test
    fun intersect() {
        assertEquals(listOf(Range(0, 100)), rangeMergingList(Range(0, 100)).intersectRanges(Range(0, 100)))
        assertEquals(listOf(Range(10, 90)), rangeMergingList(Range(0, 100)).intersectRanges(Range(10, 90)))
        assertEquals(listOf(Range(10, 90)), rangeMergingList(Range(10, 90)).intersectRanges(Range(0, 100)))
        assertEquals(
                listOf(Range(20, 40), Range(50, 70)),
                rangeMergingList(Range(0, 40), Range(50, 100)).intersectRanges(Range(20, 70))
        )
    }

    @Test
    fun lookup() {
        val data = listOf(
                Range(1, 2),
                Range(4, 6),
                Range(4, 9),
                Range(4, 7),
                Range(6, 10),
                Range(20, 30),
                Range(40, 50)
        ).toRangeMergingList()
        assertEquals(-1, data.internalLookup(0))
        assertEquals(0, data.internalLookup(1))
        assertEquals(0, data.internalLookup(3))
        assertEquals(1, data.internalLookup(4))
        assertEquals(1, data.internalLookup(5))
        assertEquals(1, data.internalLookup(10))
        assertEquals(2, data.internalLookup(30))
        assertEquals(3, data.internalLookup(100))
    }

    @Test
    fun complementaryRanges() {
        checkComplement(
                listOf(),
                listOf(Range(0, 100))
        )
        checkComplement(
                listOf(Range(4, 6)),
                listOf(Range(0, 4), Range(6, 100))
        )
        checkComplement(
                listOf(Range(3, 6), Range(14, 16), Range(40, 90)),
                listOf(Range(0, 3), Range(6, 14), Range(16, 40), Range(90, 100))
        )
        checkComplement(
                listOf(Range(1, 2), Range(4, 6), Range(10, 99)),
                listOf(Range(0, 1), Range(2, 4), Range(6, 10), Range(99, 100))
        )
        checkComplement(
                listOf(Range(0, 2), Range(4, 6), Range(10, 100)),
                listOf(Range(2, 4), Range(6, 10))
        )
        checkComplement(
                listOf(Range(0, 100)),
                listOf()
        )
    }

    private fun checkComplement(
            data: List<Range>,
            expected: List<Range>
    ) {
        val actual = data.toRangeMergingList().complementaryRanges(100).toList()
        assertEquals(expected, actual)
    }
}