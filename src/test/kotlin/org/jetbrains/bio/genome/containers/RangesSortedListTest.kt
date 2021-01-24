package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.genome.Range
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class RangesSortedListTest {
    @Test
    fun overlapRegion() {
        assertFalse(rangeSortedList().overlap(7, 12));

        val rl1 = rangeSortedList(Range(10, 30))
        assertTrue(rl1.overlap(7, 12));
        assertTrue(rl1.overlap(20, 40));

        assertFalse(rl1.overlap(2, 7));
        assertFalse(rl1.overlap(7, 10));
        assertFalse(rl1.overlap(30, 40));
        assertFalse(rl1.overlap(40, 50));

        val rl2 = rangeSortedList(
            Range(0, 40),
            Range(2, 6),
            Range(20, 40),
            Range(30, 70),
            Range(50, 100),
        )
        assertTrue(rl2.overlap(20, 70));
        assertTrue(rl2.overlap(0, 5));
        assertTrue(rl2.overlap(99, 100));
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
        assertFalse(Range(0, 12) in rangeSortedList())

        // single
        assertTrue(Range(0, 12) in rangeSortedList(Range(0, 12)))
        assertFalse(Range(0, 20) in rangeSortedList(Range(0, 12)))

        // two
        assertTrue(Range(5, 15) in rangeSortedList(Range(0, 10), Range(5, 20)))
        assertFalse(Range(5, 15) in rangeSortedList(Range(0, 10), Range(10, 20)))
        assertFalse(Range(5, 15) in rangeSortedList(Range(0, 10), Range(11, 20)))

        // multiple
        assertTrue(Range(5, 15) in rangeSortedList(Range(0, 100), Range(5, 8), Range(10, 20), Range(20, 30)))
        assertFalse(Range(5, 15) in rangeSortedList(Range(0, 10), Range(5, 8), Range(10, 20), Range(20, 30)))
        assertTrue(
            Range(5, 15) in rangeSortedList(
                Range(0, 10),
                Range(2, 20),
                Range(5, 8),
                Range(10, 20),
                Range(20, 30)
            )
        )
    }

    @Test
    fun containsSingleBigSegment() {
        assertTrue(Range(10, 20) in rangeSortedList(Range(0, 100)))
    }

    @Test
    fun containsSingleExceedingRange() {
        assertFalse(Range(-100, 100) in rangeSortedList(Range(0, 10)))
    }

    @Test
    fun intersect() {
        assertEquals(listOf(Range(0, 100)), rangeSortedList(Range(0, 100)).intersect(Range(0, 100)))
        assertEquals(listOf(Range(10, 90)), rangeSortedList(Range(0, 100)).intersect(Range(10, 90)))
        assertEquals(listOf(Range(10, 90)), rangeSortedList(Range(10, 90)).intersect(Range(0, 100)))
        assertEquals(
            listOf(Range(20, 40), Range(50, 70)),
            rangeSortedList(Range(0, 40), Range(50, 100)).intersect(Range(20, 70))
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
            rl2.intersect(Range(20, 70))
        )
    }

    @Test
    fun andWithEmpty() {
        val rl = rangeSortedList(Range(0, 10), Range(20, 30))
        assertEquals(emptyList<Range>(), (rl and rangeSortedList()).toList())
    }

    @Test
    fun andWithSelf() {
        val rl = rangeSortedList(Range(0, 10), Range(20, 30))
        assertEquals(rl.toList(), (rl and rl).toList())
    }

    @Test
    fun and() {
        val rl1 = rangeSortedList(Range(0, 10))
        val rl2 = rangeSortedList(Range(5, 30))
        assertEquals(rangeSortedList(Range(5, 10)).toList(), (rl1 and rl2).toList())

        assertEquals(
            rangeSortedList(Range(6, 8)).toList(), (
                    rangeSortedList(Range(0, 2), Range(6, 8)) and rl2).toList()
        )
    }

    @Test
    fun andComplex1() {
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
        assertEquals(
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
            ),
            rl1.and(rl2)
        )

        assertEquals(
            listOf(
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
            ),
            rl1.and(rl2).toRangeSortedList().toList()
        )
    }

    private fun checkOverlap(
        expected: List<Range>,
        a: RangesSortedList,
        b: RangesSortedList,
        flankBothSides: Int = 0
    ) {
        assertEquals(expected, (a.overlap(b, flankBothSides)).toList())
        assertEquals(expected.size, a.overlapNumber(b, flankBothSides))
    }
}