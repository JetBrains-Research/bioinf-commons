package org.jetbrains.bio.statistics

import junit.framework.TestCase.assertTrue
import org.junit.Assert.assertEquals
import org.junit.Test

class MoreMathTest {
    @Test
    fun testIntDegenerateCases() {
        assertTrue(IntArray(0).average().isNaN())
        assertTrue(IntArray(0).standardDeviation().isNaN())
    }

    @Test
    fun testIntBasicCases() {
        val data = intArrayOf(9, 17, 11, 17, 11, 18, 17, 13, 10, 12)
        assertEquals("Mean is wrong", 13.5, data.average(), 1E-5)
        assertEquals("SD is wrong", 3.407508, data.standardDeviation(), 1E-6)
    }

    @Test
    fun testLongBasicCases() {
        val data = longArrayOf(9, 17, 11, 17, 11, 18, 17, 13, 10, 12)
        assertEquals("Mean is wrong", 13.5, data.average(), 1E-5)
        assertEquals("SD is wrong", 3.407508, data.standardDeviation(), 1E-6)
    }

    @Test
    fun testShortDegenerateCases() {
        assertTrue(ShortArray(0).average().isNaN())
        assertTrue(ShortArray(0).standardDeviation().isNaN())
    }

    @Test
    fun testShortBasicCases() {
        val data = shortArrayOf(9, 17, 11, 17, 11, 18, 17, 13, 10, 12)
        assertEquals("Mean is wrong", 13.5, data.average(), 1E-5)
        assertEquals("SD is wrong", 3.407508, data.standardDeviation(), 1E-6)
    }

}
