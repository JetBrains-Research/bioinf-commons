package org.jetbrains.bio.statistics.distribution

import com.google.common.math.IntMath
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.junit.Assert.assertArrayEquals
import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CategoricalDistributionTest {
    @Before fun setUp() {
        Sampling.RANDOM_DATA_GENERATOR.randomGenerator.setSeed(12345)
    }

    @Test fun testProbability() {
        val d = CategoricalDistribution.of(values)
        for (value in d.supportLowerBound..d.supportUpperBound) {
            val expected: Double
            if (value in values) {
                expected = 1.0 / 3
            } else {
                expected = 0.0
            }

            assertTrue(Precision.equals(expected, d.probability(value), 1e-10))
        }
    }

    @Test fun testSupport() {
        val d = CategoricalDistribution.of(values)
        assertEquals(1, d.supportLowerBound)
        assertEquals(42, d.supportUpperBound)
    }

    @Test fun testSample() {
        val probabilities = doubleArrayOf(.15, .6, .25).asF64Array()
        val alias = CategoricalDistribution(probabilities)

        val workBuffer = F64Array(probabilities.size)
        for (i in 0..IntMath.pow(2, 14) - 1) {
            val value = alias.sample()
            workBuffer[value] = workBuffer[value] + 1
        }

        workBuffer.rescale()
        assertArrayEquals(probabilities.toDoubleArray(),
                          workBuffer.toDoubleArray(), 1e-2)
    }

    companion object {
        private val values = intArrayOf(42, 3, 7, 3, 7, 42)
    }
}
