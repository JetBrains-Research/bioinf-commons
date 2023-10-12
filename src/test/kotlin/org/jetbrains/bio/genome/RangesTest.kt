package org.jetbrains.bio.genome

import com.google.common.math.IntMath
import org.jetbrains.bio.genome.Strand.MINUS
import org.jetbrains.bio.genome.Strand.PLUS
import org.jetbrains.bio.genome.sequence.asNucleotideSequence
import org.junit.Test
import java.math.RoundingMode
import kotlin.test.*

class RangeTest {
    @Test
    fun testContainsOffset() {
        val range = Range(0, 10)

        assertTrue(0 in range)
        assertTrue(9 in range)
        assertFalse(10 in range)
        assertFalse(-1 in range)
    }

    @Test
    fun testIntersects() {
        assertFalse(Range(10, 11) intersects Range(11, 12))
        assertFalse(Range(11, 12) intersects Range(10, 11))
        assertTrue(Range(0, 10) intersects Range(0, 10))
        assertTrue(Range(0, 10) intersects Range(5, 15))
        assertTrue(Range(5, 15) intersects Range(0, 10))
        assertTrue(Range(0, 10) intersects Range(4, 5))
        assertTrue(Range(4, 5) intersects Range(0, 10))
    }

    @Test
    fun testDistanceTo() {
        assertEquals(10, Range(0, 10) distanceTo Range(20, 30))
        assertEquals(10, Range(20, 30) distanceTo Range(0, 10))
        assertEquals(0, Range(0, 10) distanceTo Range(10, 20))
        assertEquals(0, Range(0, 10) distanceTo Range(5, 15))
    }


    @Test
    fun testIntersection() {
        assertNull(Range(10, 11) intersection Range(11, 12))
        assertEquals(Range(11, 12), Range(10, 12) intersection Range(11, 12))
    }

    @Test
    fun testSlice() {
        var bins = 0
        var prevBound = 0
        for (r in Range(0, 103).slice(10)) {
            assertEquals(prevBound, r.startOffset)
            prevBound = r.endOffset
            assertTrue(r.length() <= 10)
            bins++
        }

        assertEquals(11, bins)
        assertEquals(103, prevBound)
        assertEquals(IntMath.divide(103, 10, RoundingMode.CEILING), bins)
        assertEquals(0, Range(0, 0).slice(150).count())
        assertEquals(1, Range(0, 103).slice(150).count())
        assertEquals(1, Range(0, 103).slice(103).count())
        assertEquals(2, Range(0, 103).slice(53).count())
        assertEquals(3, Range(0, 103).slice(50).count())
    }
}

class LocationTest {
    private var chromosome = Chromosome(Genome["to1"], "chr1")

    @Test
    fun testContainsOffset() {
        val loc = Location(0, 10, chromosome, PLUS)

        assertTrue(0 in loc)
        assertTrue(9 in loc)
        assertFalse(10 in loc)
        assertFalse(-1 in loc)
    }

    @Test
    fun testGetLength_EmptyLocation() {
        assertEquals(0, Location(0, 0, chromosome, PLUS).length())
        assertEquals(0, Location(1, 1, chromosome, PLUS).length())
    }

    @Test
    fun testGetLength() {
        assertEquals(2, Location(0, 2, chromosome, PLUS).length())
        assertEquals(2, Location(1, 3, chromosome, PLUS).length())
    }

    @Test
    fun testConvert5BoundRelativeToAbsoluteOffset() {
        assertEquals(0, Location(1, 5, chromosome, PLUS).get5Bound(-1))
        assertEquals(1, Location(1, 5, chromosome, PLUS).get5Bound(0))
        assertEquals(2, Location(1, 5, chromosome, PLUS).get5Bound(1))

        assertEquals(5, Location(1, 5, chromosome, Strand.MINUS).get5Bound(-1))
        assertEquals(4, Location(1, 5, chromosome, Strand.MINUS).get5Bound(0))
        assertEquals(3, Location(1, 5, chromosome, Strand.MINUS).get5Bound(1))
    }

    @Test
    fun testConvert53BoundRelativeToAbsoluteOffset() {
        assertEquals(3, Location(1, 5, chromosome, PLUS).get3Bound(-1))
        assertEquals(4, Location(1, 5, chromosome, PLUS).get3Bound(0))
        assertEquals(5, Location(1, 5, chromosome, PLUS).get3Bound(1))

        assertEquals(2, Location(1, 5, chromosome, Strand.MINUS).get3Bound(-1))
        assertEquals(1, Location(1, 5, chromosome, Strand.MINUS).get3Bound(0))
        assertEquals(0, Location(1, 5, chromosome, Strand.MINUS).get3Bound(1))
    }

    @Test
    fun testGetSequence() {
        val plusSequence = Location(1, 10, chromosome, PLUS).sequence
        val rcSequence = plusSequence.asNucleotideSequence()
            .substring(0, plusSequence.length, Strand.MINUS)
        assertEquals(rcSequence, Location(1, 10, chromosome, Strand.MINUS).sequence)
    }

    @Test
    fun testAroundStart_Plus() {
        val location = Location(100, 200, chromosome, PLUS)
        assertEquals(
            "chr1:+[100, 106)",
            RelativePosition.FIVE_PRIME.of(location, 0, 6).toString()
        )
        assertEquals(
            "chr1:+[100, 200)",
            RelativePosition.FIVE_PRIME.of(location, 0, 100).toString()
        )
        assertEquals(
            "chr1:+[0, 2100)", // Edge case.
            RelativePosition.FIVE_PRIME.of(location, -2000, 2000).toString()
        )
    }

    @Test
    fun testAroundStart_Minus() {
        val location = Location(100, 200, chromosome, Strand.MINUS)
        assertEquals(
            "chr1:-[194, 200)",
            RelativePosition.FIVE_PRIME.of(location, 0, 6).toString()
        )
        assertEquals(
            "chr1:-[100, 200)",
            RelativePosition.FIVE_PRIME.of(location, 0, 100).toString()
        )
    }

    @Test
    fun testAroundEnd_Plus() {
        val location = Location(100, 200, chromosome, PLUS)
        assertEquals(
            "chr1:+[199, 205)",
            RelativePosition.THREE_PRIME.of(location, 0, 6).toString()
        )
        assertEquals(
            "chr1:+[100, 200)",
            RelativePosition.THREE_PRIME.of(location, -99, 1).toString()
        )
    }

    @Test
    fun testAroundEnd_Minus() {
        val location = Location(100, 200, chromosome, Strand.MINUS)
        assertEquals(
            "chr1:-[95, 101)",
            RelativePosition.THREE_PRIME.of(location, 0, 6).toString()
        )
        assertEquals(
            "chr1:-[100, 200)",
            RelativePosition.THREE_PRIME.of(location, -99, 1).toString()
        )
    }

    @Test
    fun testAroundWhole_Plus() {
        val location = Location(100, 200, chromosome, PLUS)
        assertEquals(
            "chr1:+[100, 200)",
            RelativePosition.ALL.of(location, 0, 1).toString()
        )
    }

    @Test
    fun testAroundWhole_End() {
        val location = Location(100, 200, chromosome, Strand.MINUS)
        assertEquals(
            "chr1:-[100, 200)",
            RelativePosition.ALL.of(location, 0, 1).toString()
        )
    }

    @Test
    fun testComparator() {
        val location1 = Location(0, 100, chromosome, PLUS)
        val location2 = Location(0, 100, chromosome, Strand.MINUS)
        assertNotEquals(0, location1.compareTo(location2))

        val location3 = Location(0, 100, chromosome, PLUS)
        assertEquals(0, location1.compareTo(location3))

        val location4 = Location(0, 100, Chromosome(Genome["to1"], "chr2"), PLUS)
        assertEquals(-1, location1.compareTo(location4))
        assertEquals(1, location4.compareTo(location1))
    }

    @Test
    fun testIntersectsWithNull() {
        assertFalse(Location.intersects(null, null))

        assertFalse(Location.intersects(null, Location(1, 5, chromosome, PLUS)))
        assertFalse(Location.intersects(Location(1, 5, chromosome, PLUS), null))
    }
    @Test
    fun testIntersectsWithOtherChr() {
        assertFalse(
            Location.intersects(
                Location(1, 5, chromosome, PLUS),
                Location(1, 5, Chromosome(Genome["to1"], "chr2"), PLUS)
            )
        )
    }

    @Test
    fun testIntersectsWithOtherStand() {
        assertTrue(Location.intersects(
            Location(11, 15, chromosome, PLUS),
            Location(11, 15, chromosome, MINUS)
        ))
        assertTrue(Location.intersects(
            Location(11, 15, chromosome, PLUS),
            Location(14, 18, chromosome, MINUS)
        ))
        assertTrue(Location.intersects(
            Location(4, 12, chromosome, MINUS),
            Location(11, 15, chromosome, PLUS),
        ))
        assertFalse(Location.intersects(
            Location(11, 15, chromosome, PLUS),
            Location(15, 18, chromosome, MINUS)
        ))
        assertFalse(Location.intersects(
            Location(11, 15, chromosome, PLUS),
            Location(4, 11, chromosome, MINUS)
        ))
    }
    @Test
    fun testIntersects() {
        val l1 = Location(11, 15, chromosome, PLUS)
        assertTrue(Location.intersects(l1, l1))
        assertTrue(Location.intersects(l1, Location(14, 16, chromosome, PLUS)))
        assertFalse(Location.intersects(l1, Location(15, 16, chromosome, PLUS)))
        assertTrue(Location.intersects(l1, Location(4, 12, chromosome, PLUS)))
        assertFalse(Location.intersects(l1, Location(4, 11, chromosome, PLUS)))

        val l2 = Location(11, 15, chromosome, MINUS)
        assertTrue(Location.intersects(l2, l2))
        assertTrue(Location.intersects(l2, Location(14, 16, chromosome, MINUS)))
        assertFalse(Location.intersects(l2, Location(15, 16, chromosome, MINUS)))
        assertTrue(Location.intersects(l2, Location(4, 12, chromosome, MINUS)))
        assertFalse(Location.intersects(l2, Location(4, 11, chromosome, MINUS)))
    }

    @Test
    fun testChromosomeRangeIntersectsWithNull() {
        assertFalse(ChromosomeRange.intersects(null, null))
        assertFalse(ChromosomeRange.intersects(null, ChromosomeRange(1, 5, chromosome)))
        assertFalse(ChromosomeRange.intersects(ChromosomeRange(1, 5, chromosome), null))
    }
    @Test
    fun testChromosomeRangeIntersectsWithOtherChr() {
        assertFalse(
            ChromosomeRange.intersects(
                ChromosomeRange(1, 5, chromosome),
                ChromosomeRange(1, 5, Chromosome(Genome["to1"], "chr2"))
            )
        )
    }

    @Test
    fun testChromosomeRangeIntersects() {
        val range1 = ChromosomeRange(11, 15, chromosome)

        assertTrue(ChromosomeRange.intersects(range1, range1))
        assertTrue(ChromosomeRange.intersects(range1, ChromosomeRange(14, 16, chromosome)))
        assertFalse(ChromosomeRange.intersects(range1, ChromosomeRange(15, 16, chromosome)))
        assertTrue(ChromosomeRange.intersects(range1, ChromosomeRange(4, 12, chromosome)))
        assertFalse(ChromosomeRange.intersects(range1, ChromosomeRange(4, 11, chromosome)))
    }

}
