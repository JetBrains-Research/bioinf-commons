package org.jetbrains.bio.statistics.hypothesis

import org.junit.Assert
import org.junit.Test

class StofferLiptakTestTest {
    @Test
    fun testCorrelations() {
        val  correlations = StofferLiptakTest.computeCorrelations(doubleArrayOf(0.0, 0.1, 0.2, 0.3, 0.4), 20)
        Assert.assertEquals(3, correlations.size)
        Assert.assertTrue(doubleArrayOf(0.0, -0.48507125007266605, -0.7298004491997617) contentEquals correlations)
    }

    @Test
    fun testZeroesCorrelations() {
        val correlations = StofferLiptakTest.computeCorrelations(doubleArrayOf(0.0, 0.0, 0.0, 0.0), 2)
        Assert.assertEquals(3, correlations.size)
        Assert.assertTrue(doubleArrayOf(0.0, 0.0, 0.0) contentEquals correlations)
    }

    @Test
    fun testZScore() {
        val stofferLiptakTest = StofferLiptakTest(doubleArrayOf(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.2, 0.3))
        Assert.assertEquals(7.650730905155645, stofferLiptakTest.zscore(0.0), 1e-6)
        Assert.assertEquals(-7.650730905155645, stofferLiptakTest.zscore(1.0), 1e-6)
        Assert.assertEquals(1.6448536269514724, stofferLiptakTest.zscore(0.05), 1e-6)
    }

    @Test
    fun testCombine() {
        val stofferLiptakTest = StofferLiptakTest(doubleArrayOf(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.2, 0.3))
        // Check that results are monotonous by pvalues and size
        Assert.assertEquals(1e-6, stofferLiptakTest.combine(doubleArrayOf(1e-6)), 1e-6)
        Assert.assertEquals(0.004272487313437656, stofferLiptakTest.combine(doubleArrayOf(1e-6, 1e-4)), 1e-6)
        Assert.assertEquals(0.1824166413929258, stofferLiptakTest.combine(doubleArrayOf(1e-2, 0.1)), 1e-6)
        Assert.assertEquals(4.988280790652055E-7, stofferLiptakTest.combine(doubleArrayOf(1e-8, 1e-6, 1e-4)), 1e-6)
        Assert.assertEquals(
            3.207676690930583E-8, stofferLiptakTest.combine(doubleArrayOf(1e-10, 1e-8, 1e-6, 1e-4)), 1e-6
        )
        Assert.assertEquals(
            1.1102230246251565E-16, stofferLiptakTest.combine(doubleArrayOf(1e-10, 1e-8, 1e-10, 1e-4)), 1e-12
        )
        Assert.assertEquals(
            1.0E-14, stofferLiptakTest.combine(doubleArrayOf(1e-12, 1e-10, 1e-8, 1e-10, 1e-4)), 1e-12
        )
    }

}