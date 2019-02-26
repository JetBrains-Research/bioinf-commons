package org.jetbrains.bio.statistics

import junit.framework.TestCase.assertTrue
import org.apache.commons.math3.util.Precision
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

class PrefixSumTableTest {
    private val values = doubleArrayOf(1e-5, -42.0, 2.4, 1.25, 5.0, -1.0, 2.0, 7.0, -9.0)

    @Test
    fun testPrefixSum() {
        val pst = PrefixSumTable(values)
        for (i in values.indices) {
            val sub = values.copyOfRange(0, i + 1)
            kotlin.test.assertTrue(Precision.equals(sub.sum(), pst[i], 5))
        }
    }

    @Test
    fun testIntervalSum() {
        val pst = PrefixSumTable(values)
        for (i in values.indices) {
            for (j in i until values.size) {
                val sub = values.copyOfRange(i, j + 1)
                kotlin.test.assertTrue(Precision.equals(sub.sum(), pst[i, j], 5))
            }
        }
    }
}

