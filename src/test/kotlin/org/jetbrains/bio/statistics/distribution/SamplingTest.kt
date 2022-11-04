package org.jetbrains.bio.statistics.distribution

import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import gnu.trove.map.hash.TObjectIntHashMap
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath
import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class SamplingTest {
    @Before
    fun setUp() {
        Sampling.RANDOM_DATA_GENERATOR.reSeed(42)
    }

    @Test
    fun testSampleCombination() {
        val r = Sampling.RANDOM_DATA_GENERATOR.randomGenerator

        val n = r.nextInt(6) + 10     // n \in [10, 15]
        val k = n - r.nextInt(6) - 5  // k \in [5 , 10]
        val numPartitions = CombinatoricsUtils.binomialCoefficient(n - 1, k - 1).toInt()
        val numRuns = numPartitions * numPartitions

        val countsMap = TObjectIntHashMap<TIntList>()
        for (i in 0 until numRuns) {
            countsMap.adjustOrPutValue(TIntArrayList(Sampling.sampleCombination(n - 1, k - 1)), 1, 1)
        }

        val counts = TIntArrayList(countsMap.valueCollection())
        assertEquals(numPartitions, counts.size())
        for (i in 0 until numPartitions) {
            val p = counts[i].toDouble() / numRuns
            assertTrue(FastMath.abs(p - 1.0 / numPartitions) <= 0.1)
        }
    }

    @Test
    fun testSampleGamma() {
        // Edge case: shape = 0.
        assertEquals(0.0, Sampling.sampleGamma(0.0, 10.0))
    }
}