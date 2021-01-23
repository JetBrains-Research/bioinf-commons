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
        assertFalse(Range(0, 12) in rangeMergingList())
        assertEquals(0, rangeMergingList().intersectionLength(Range(0, 12)))
    }

    @Test
    fun single() {
        assertTrue(Range(0, 12) in rangeMergingList(Range(0, 12)))
        assertFalse(Range(0, 20) in rangeMergingList(Range(0, 12)))
        assertEquals(12, rangeMergingList(Range(0, 12)).intersectionLength(Range(0, 12)))
    }

    @Test
    fun singleBigSegment() {
        assertTrue(Range(10, 20) in rangeMergingList(Range(0, 100)))
    }

    @Test
    fun singleExceedingRange() {
        assertFalse(Range(-100, 100) in rangeMergingList(Range(0, 10)))
    }

    @Test
    fun intersectionLength() {
        val rl = rangeMergingList(Range(0, 10))
        assertFalse(Range(5, 15) in rl)
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
        assertTrue(Range(5, 15) in rangeMergingList(Range(0, 10), Range(10, 20)))
        assertFalse(Range(5, 15) in rangeMergingList(Range(0, 10), Range(11, 20)))
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
        assertTrue(Range(11, 12) in rangeMergingList(Range(10, 20), Range(0, 5)))
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
    fun andWithEmpty() {
        val rl = rangeMergingList(Range(0, 10), Range(20, 30))
        assertEquals(emptyList<Range>(), (rl and rangeMergingList()).toList())
    }

    @Test
    fun andWithSelf() {
        val rl = rangeMergingList(Range(0, 10), Range(20, 30))
        assertEquals(rl.toList(), (rl and rl).toList())
    }

    @Test
    fun and() {
        val rl1 = rangeMergingList(Range(0, 10))
        val rl2 = rangeMergingList(Range(5, 30))
        assertEquals(rangeMergingList(Range(5, 10)).toList(), (rl1 and rl2).toList())
        assertEquals(
            rangeMergingList(Range(5, 10), Range(25, 30)).toList(),
            ((rl1 or rangeMergingList(Range(25, 30))) and rl2).toList()
        )
        assertEquals(
            rangeMergingList(Range(5, 10), Range(25, 30)).toList(),
            (rl2 and (rl1 or rangeMergingList(Range(25, 30)))).toList()
        )

        assertEquals(
            rangeMergingList(Range(6, 8)).toList(), (
                    rangeMergingList(Range(0, 2), Range(6, 8)) and rl2).toList()
        )
    }

    @Test
    fun overlap() {
        val rl1 = rangeMergingList(Range(10, 30))

        assertEquals(rl1.toList(), (rl1.overlap(rangeMergingList(Range(7, 12)))).toList())
        assertEquals(rl1.toList(), (rl1.overlap(rangeMergingList(Range(2, 7), Range(20, 40)))).toList())
        assertEquals(listOf(), (rl1.overlap(rangeMergingList(Range(2, 7)))).toList())
    }


    @Test
    fun overlapFlnk() {
        val rl1 = rangeMergingList(Range(10, 30), Range(50, 100))

        assertEquals(listOf(Range(10, 30)), (rl1.overlap(rangeMergingList(Range(7, 12)), 10)).toList())
        assertEquals(listOf(Range(10, 30)), (rl1.overlap(rangeMergingList(Range(35, 40)), 10)).toList())
        assertEquals(rl1.toList(), (rl1.overlap(rangeMergingList(Range(35, 40)), 11)).toList())
        assertEquals(listOf(Range(10, 30)), (rl1.overlap(rangeMergingList(Range(2, 7), Range(20, 40)), 2)).toList())
        assertEquals(listOf(), (rl1.overlap(rangeMergingList(Range(2, 7)), 3)).toList())
        assertEquals(listOf(Range(10, 30)), (rl1.overlap(rangeMergingList(Range(2, 7)), 4)).toList())
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
        assertEquals(listOf(Range(0, 100)), rangeMergingList(Range(0, 100)).intersect(Range(0, 100)))
        assertEquals(listOf(Range(10, 90)), rangeMergingList(Range(0, 100)).intersect(Range(10, 90)))
        assertEquals(listOf(Range(10, 90)), rangeMergingList(Range(10, 90)).intersect(Range(0, 100)))
        assertEquals(
            listOf(Range(20, 40), Range(50, 70)),
            rangeMergingList(Range(0, 40), Range(50, 100)).intersect(Range(20, 70))
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
}