package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class LocationsMergingListTest {
    @Test fun containsSingleLocus() {
        val locationList = testList(Location(10, 100, chromosome),
                                    Location(300, 400, chromosome))

        assertTrue(Location(20, 90, chromosome) in locationList)
    }

    @Test fun containsOppositeStrand() {
        val locationList = testList(Location(10, 100, chromosome))
        assertFalse(Location(20, 90, chromosome, Strand.MINUS) in locationList)
    }

    @Test fun containsIntersectedLocus() {
        val locationList = testList(Location(10, 100, chromosome))
        assertTrue(Location(20, 100, chromosome) in locationList)
    }

    @Test fun containsNonIntersectedLocus() {
        val locationList = testList(Location(10, 100, chromosome))
        assertFalse(Location(100, 200, chromosome) in locationList)
        assertFalse(Location(110, 200, chromosome) in locationList)
    }

    @Test fun containsSameLocus() {
        val location = Location(10, 100, chromosome)
        assertTrue(location in testList(location))
    }

    @Test fun containsSeveralLoci() {
        val locationList = testList(Location(10, 100, chromosome),
                                    Location(300, 400, chromosome))

        assertTrue(Location(20, 90, chromosome) in locationList)
        assertTrue(Location(320, 390, chromosome) in locationList)
    }

    @Test fun containsSeveralLociNonSorted() {
        val locationList = testList(Location(300, 400, chromosome),
                                    Location(10, 100, chromosome))

        assertTrue(Location(20, 90, chromosome) in locationList)
        assertTrue(Location(320, 390, chromosome) in locationList)
        assertFalse(Location(200, 250, chromosome) in locationList)
        assertFalse(Location(0, 5, chromosome) in locationList)
        assertFalse(Location(410, 420, chromosome) in locationList)
    }

    @Test fun testCollapse() {
        val locationList = testList(Location(10, 20, chromosome),
                                    Location(20, 30, chromosome))
        val get = locationList.get(chromosome, Strand.PLUS)
        assertTrue(Location(10, 30, chromosome) in get)
    }

    @Test fun containsEmpty() {
        assertFalse(Location(0, 100, chromosome) in testList())
    }

    @Test fun testLocationMinus() {
        val universe = Location(0, 100, chromosome, Strand.PLUS)
        assertEquals(listOf(universe), universe - emptyList())
        assertEquals(listOf(universe),
                     universe - listOf(Location(0, 40, chromosome, Strand.MINUS)))
        assertEquals(listOf(universe),
                     universe - listOf(Location(0, 40, Chromosome("to1", "chr2"), Strand.PLUS)))

        assertEquals(listOf(Location(0, 20, chromosome, Strand.PLUS),
                            Location(40, 100, chromosome, Strand.PLUS)),
                     universe - listOf(Location(20, 40, chromosome, Strand.PLUS)))
    }

    @Test fun testLocationBothStrand() {
        val locationList = testList(Location(10, 100, chromosome), Location(300, 400, chromosome))
        assertTrue(locationList.intersectsBothStrands(Location(20, 30, chromosome, Strand.MINUS)))
        assertEquals(2, locationList.intersectBothStrands(Location(90, 310, chromosome, Strand.MINUS)).size)
    }

    companion object {
        private val chromosome = Chromosome("to1", "chr1")
    }
}

private fun testList(vararg locations: Location): LocationsMergingList {
    return locationList(GenomeQuery("to1"), *locations)
}
