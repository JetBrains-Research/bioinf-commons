package org.jetbrains.bio.gsea

import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.util.Precision
import org.junit.Assert
import org.junit.Test
import kotlin.math.sqrt
import kotlin.test.assertTrue

class IntHistogramTest {
    @Test
    fun countValues() {
        Assert.assertEquals(0, IntHistogram.create(intArrayOf()).countValues())
        Assert.assertEquals(1, IntHistogram.create(intArrayOf(2)).countValues())
        Assert.assertEquals(8, IntHistogram.create(intArrayOf(3, 4, 5, 5, 5, 6, 7, 3)).countValues())
    }

    @Test
    fun mean() {
        doCheckMean(intArrayOf())
        doCheckMean(intArrayOf(1))
        doCheckMean(intArrayOf(1, 2))
        doCheckMean(intArrayOf(3, 4, 5, 5, 5, 6, 7, 101))
    }

    @Test
    fun variance() {
        doCheckStdev(intArrayOf())
        doCheckStdev(intArrayOf(1))
        doCheckStdev(intArrayOf(1, 2))
        doCheckStdev(intArrayOf(3, 4, 5, 5, 5, 6, 7, 101))
    }

    @Test
    fun median() {
        doCheckMedian(intArrayOf(), -1)
        doCheckMedian(intArrayOf(1), 1)
        doCheckMedian(intArrayOf(1, 2), 2)
        doCheckMedian(intArrayOf(1, 2, 3), 2)
        doCheckMedian(intArrayOf(1, 2, 2, 3), 2)
        doCheckMedian(intArrayOf(0, 0, 0, 2, 2, 2, 2, 2, 3, 4, 4), 2)
        doCheckMedian(intArrayOf(0, 0, 0, 2, 2, 2, 3, 4, 4, 4, 4), 2)
        doCheckMedian(intArrayOf(0, 0, 0, 2, 2, 3, 4, 4, 4, 4, 4), 3)
        doCheckMedian(intArrayOf(3, 4, 5, 5, 5, 6, 7, 101), 5)
    }

    @Test
    fun mergeHists() {
        checkMergeHists(intArrayOf(), intArrayOf(1, 1, 4, 3, 3), intArrayOf(1, 1, 4, 3, 3))
        checkMergeHists(intArrayOf(1, 1, 4, 3, 3), intArrayOf(), intArrayOf(1, 1, 4, 3, 3))
        checkMergeHists(intArrayOf(1, 1, 2, 4), intArrayOf(1, 1, 4, 3, 3), intArrayOf(1, 1, 2, 4, 1, 1, 4, 3, 3))
    }

    private fun checkMergeHists(src: IntArray, target: IntArray, expected: IntArray) {
        val srcHist = IntHistogram.create(src)
        val targetHist = IntHistogram.create(target)
        targetHist += srcHist
        Assert.assertEquals(
            IntHistogram.create(expected), targetHist
        )
    }

    fun assertEquals(expected: Double, actual: Double, eps: Double = 0.0001) {
        assertTrue(
            Precision.equals(expected, actual, eps),
            "Expected: $expected, but was $actual (precision: $eps)"
        )
    }

    private fun doCheckMean(data: IntArray) {
        val metricHist = IntHistogram.create(data)
        val expected = StatUtils.mean(data.map { it.toDouble() }.toDoubleArray())
        val n = metricHist.countValues()
        val actual = metricHist.mean(n)
        if (expected.isNaN()) {
            assertTrue(actual.isNaN())
        } else {
            assertEquals(expected, actual)
        }
    }

    private fun doCheckMedian(data: IntArray, expected: Int) {
        val metricHist = IntHistogram.create(data)
        val n = metricHist.countValues()
        val actual = metricHist.median(n)
        Assert.assertEquals(expected, actual)
        Assert.assertEquals(actual, metricHist.median())
    }

    private fun doCheckStdev(data: IntArray) {
        val metricHist = IntHistogram.create(data)

        val expected = sqrt(StatUtils.variance(data.map { it.toDouble() }.toDoubleArray()))

        val n = metricHist.countValues()
        val actualMean = metricHist.mean(n)
        val actual = metricHist.stdev(n, actualMean)

        if (expected.isNaN()) {
            assertTrue(actual.isNaN())
            assertTrue(metricHist.stdev().isNaN())
            assertTrue(metricHist.stdev(n).isNaN())
        } else {
            assertEquals(expected, actual)
            assertEquals(actual, metricHist.stdev())
            assertEquals(actual, metricHist.stdev(n))
        }
    }
}