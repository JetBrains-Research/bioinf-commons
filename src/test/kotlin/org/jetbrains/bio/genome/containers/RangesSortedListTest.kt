package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.genome.Range
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class RangesSortedListTest {
    @Test
    fun overlapRegion() {
        assertFalse(rangeSortedList().overlapRanges(7, 12));

        val rl1 = rangeSortedList(Range(10, 30))
        assertTrue(rl1.overlapRanges(7, 12));
        assertTrue(rl1.overlapRanges(20, 40));

        assertFalse(rl1.overlapRanges(2, 7));
        assertFalse(rl1.overlapRanges(7, 10));
        assertFalse(rl1.overlapRanges(30, 40));
        assertFalse(rl1.overlapRanges(40, 50));

        val rl2 = rangeSortedList(
            Range(0, 40),
            Range(2, 6),
            Range(20, 40),
            Range(30, 70),
            Range(50, 100),
        )
        assertTrue(rl2.overlapRanges(20, 70));
        assertTrue(rl2.overlapRanges(0, 5));
        assertTrue(rl2.overlapRanges(99, 100));
    }

    @Test
    fun overlap() {
        val rl1 = rangeSortedList(Range(10, 30))
        checkOverlap(rl1.toList(), rl1, rangeSortedList(Range(7, 12)))
        checkOverlap(rl1.toList(), rl1, rangeSortedList(Range(2, 7), Range(20, 40)))
        checkOverlap(listOf(), rl1, rangeSortedList(Range(2, 7)))

        val rl2 = rangeSortedList(
            Range(0, 40),
            Range(2, 6),
            Range(9, 15),
            Range(10, 30),
            Range(10, 50),
            Range(20, 40),
            Range(20, 50),
            Range(30, 70),
            Range(50, 100),
            Range(60, 70),
            Range(90, 140)
        )
        checkOverlap(
            listOf(
                Range(0, 40),
                Range(10, 30),
                Range(10, 50),
                Range(20, 40),
                Range(20, 50),
                Range(30, 70),
                Range(50, 100),
                Range(60, 70),
            ), rl2, rangeSortedList(Range(20, 70))
        )
    }

    @Test
    fun overlapFlnk() {
        val rl1 = rangeSortedList(Range(10, 30), Range(50, 100))

        checkOverlap(listOf(Range(10, 30)), rl1, rangeSortedList(Range(7, 12)), 10)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeSortedList(Range(35, 40)), 10)
        checkOverlap(rl1.toList(), rl1, rangeSortedList(Range(35, 40)), 11)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeSortedList(Range(2, 7), Range(20, 40)), 2)
        checkOverlap(listOf(), rl1, rangeSortedList(Range(2, 7)), 3)
        checkOverlap(listOf(Range(10, 30)), rl1, rangeSortedList(Range(2, 7)), 4)

        val rl2 = rangeSortedList(
            Range(0, 40),
            Range(2, 6),
            Range(9, 15),
            Range(10, 30),
            Range(10, 50),
            Range(20, 40),
            Range(20, 50),
            Range(30, 70),
            Range(50, 100),
            Range(60, 70),
            Range(190, 240),
            Range(90, 140)
        )

        checkOverlap(
            listOf(
                Range(0, 40), Range(10, 30), Range(10, 50), Range(20, 40),
                Range(20, 50), Range(30, 70), Range(50, 100), Range(60, 70),
            ), rl2, rangeSortedList(Range(20, 70)), 5
        )

        checkOverlap(
            listOf(
                Range(0, 40), Range(9, 15), Range(10, 30), Range(10, 50), Range(20, 40),
                Range(20, 50), Range(30, 70), Range(50, 100), Range(60, 70),
            ), rl2, rangeSortedList(Range(20, 70)), 6
        )

        checkOverlap(
            listOf(
                Range(0, 40), Range(2, 6), Range(9, 15), Range(10, 30), Range(10, 50), Range(20, 40),
                Range(20, 50), Range(30, 70), Range(50, 100), Range(60, 70),
            ), rl2, rangeSortedList(Range(20, 70)), 15
        )

        checkOverlap(
            listOf(
                Range(0, 40), Range(2, 6), Range(9, 15), Range(10, 30), Range(10, 50), Range(20, 40),
                Range(20, 50), Range(30, 70), Range(50, 100), Range(60, 70), Range(90, 140)
            ), rl2, rangeSortedList(Range(20, 70)), 29
        )
    }

    @Test
    fun contains() {
        assertFalse(rangeSortedList().includesRange(Range(0, 12)))

        // single
        assertTrue(rangeSortedList(Range(0, 12)).includesRange(Range(0, 12)))
        assertFalse(rangeSortedList(Range(0, 12)).includesRange(Range(0, 20)))

        // two
        assertTrue(rangeSortedList(Range(0, 10), Range(5, 20)).includesRange(Range(5, 15)))
        assertFalse(rangeSortedList(Range(0, 10), Range(10, 20)).includesRange(Range(5, 15)))
        assertFalse(rangeSortedList(Range(0, 10), Range(11, 20)).includesRange(Range(5, 15)))

        // multiple
        assertTrue(rangeSortedList(Range(0, 100), Range(5, 8), Range(10, 20), Range(20, 30)).includesRange(Range(5, 15)))
        assertFalse(rangeSortedList(Range(0, 10), Range(5, 8), Range(10, 20), Range(20, 30)).includesRange(Range(5, 15)))
        assertTrue(
            rangeSortedList(
                Range(0, 10),
                Range(2, 20),
                Range(5, 8),
                Range(10, 20),
                Range(20, 30)
            ).includesRange(Range(5, 15))
        )
    }

    @Test
    fun containsSingleBigSegment() {
        assertTrue(rangeSortedList(Range(0, 100)).includesRange(Range(10, 20)))
    }

    @Test
    fun containsSingleExceedingRange() {
        assertFalse(rangeSortedList(Range(0, 10)).includesRange(Range(-100, 100)))
    }

    @Test
    fun intersect() {
        assertEquals(listOf(Range(0, 100)), rangeSortedList(Range(0, 100)).intersectRanges(Range(0, 100)))
        assertEquals(listOf(Range(10, 90)), rangeSortedList(Range(0, 100)).intersectRanges(Range(10, 90)))
        assertEquals(listOf(Range(10, 90)), rangeSortedList(Range(10, 90)).intersectRanges(Range(0, 100)))
        assertEquals(
            listOf(Range(20, 40), Range(50, 70)),
            rangeSortedList(Range(0, 40), Range(50, 100)).intersectRanges(Range(20, 70))
        )

        val rl2 = rangeSortedList(
            Range(0, 40),
            Range(2, 6),
            Range(9, 15),
            Range(10, 30),
            Range(10, 50),
            Range(20, 40),
            Range(20, 50),
            Range(30, 70),
            Range(50, 100),
            Range(60, 70),
            Range(70, 80),
            Range(90, 140)
        )
        assertEquals(
            listOf(
                Range(20, 40),
                Range(20, 30),
                Range(20, 50),
                Range(20, 40),
                Range(20, 50),
                Range(30, 70),
                Range(50, 70),
                Range(60, 70),
            ),
            rl2.intersectRanges(Range(20, 70))
        )
    }

    @Test
    fun intersectRangesWithEmpty() {
        val rl = rangeSortedList(Range(0, 10), Range(20, 30))
        doCheckIntersectedRanges(emptyList(), rl, rangeSortedList())
    }

    @Test
    fun intersectRangesWithSelf() {
        val rl = rangeSortedList(Range(0, 10), Range(20, 30))
        doCheckIntersectedRanges(rl.toList(), rl, rl)
    }

    @Test
    fun intersectRanges() {
        val rl1 = rangeSortedList(Range(0, 10))
        val rl2 = rangeSortedList(Range(5, 30))

        doCheckIntersectedRanges(rangeSortedList(Range(5, 10)).toList(), rl1, rl2)
        doCheckIntersectedRanges(rangeSortedList(Range(6, 8)).toList(), rangeSortedList(Range(0, 2), Range(6, 8)), rl2)
    }

    @Test
    fun intersectRangesComplex1() {
        val rl1 = rangeSortedList(
                    Range(0, 40),
                    Range(2, 6),
                    Range(9, 15),
                    Range(10, 30),
                    Range(10, 50),
                    Range(20, 40),
                    Range(20, 50),
                    Range(30, 70),
                    Range(50, 100),
                    Range(60, 70),
                    Range(70, 80),
                    Range(90, 140)
                )
        val rl2 = rangeSortedList(Range(0, 5), Range(0, 90), Range(10, 15))
        // List with duplicated!
        doCheckIntersectedRanges(
            listOf(
                Range(0, 5),
                Range(0, 40),
                Range(10, 15),
                Range(2, 5),
                Range(2, 6),
                Range(9, 15),
                Range(10, 15),
                Range(10, 30),
                Range(10, 15),
                Range(10, 50),
                Range(10, 15),
                Range(20, 40),
                Range(20, 50),
                Range(30, 70),
                Range(50, 90),
                Range(60, 70),
                Range(70, 80)
            ), rl1, rl2
        )


        val expected = listOf(
            Range(0, 5),
            Range(0, 40),
            Range(2, 5),
            Range(2, 6),
            Range(9, 15),
            Range(10, 15),
            Range(10, 15),
            Range(10, 15),
            Range(10, 15),
            Range(10, 30),
            Range(10, 50),
            Range(20, 40),
            Range(20, 50),
            Range(30, 70),
            Range(50, 90),
            Range(60, 70),
            Range(70, 80)
        )
        doCheckIntersectedRanges(expected, rl1, rl2, true)
    }

    private fun doCheckIntersectedRanges(
        expected: List<Range>,
        rl1: RangesSortedList,
        rl2: RangesSortedList,
        sortActual: Boolean = false
    ) {
        assertEquals(
            expected,
            rl1.intersectRanges(rl2).let {
                if (sortActual) {
                    it.toRangeSortedList().toList()
                }  else {
                    it
                }
            }
        )
        assertEquals(
            expected.size, rl1.intersectRangesNumber(rl2)
        )
    }

    private fun checkOverlap(
        expected: List<Range>,
        a: RangesSortedList,
        b: RangesSortedList,
        flankBothSides: Int = 0
    ) {
        assertEquals(expected, (a.overlapRanges(b, flankBothSides)).toList())
        assertEquals(expected.size, a.overlapRangesNumber(b, flankBothSides))
    }
}