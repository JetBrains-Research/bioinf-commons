package org.jetbrains.bio.statistics.hypothesis

import com.google.common.math.IntMath
import org.apache.commons.math3.distribution.BetaDistribution
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class FdrTest {
    @Test
    fun testQValueEstimation() {
        val numTests = IntMath.pow(2, 16)
        val logMemberships = BetaDistribution(1.0, 1.0).sample(numTests)
            .asF64Array().log()

        val alphas = F64Array.of(0.0001, 0.001, 0.01, 0.1)
        for (i in 0 until alphas.size) {
            val directRejected = Fdr.control(logMemberships, alphas[i]).cardinality()
            val qvalueRejected = Fdr.qvalidate(logMemberships).toDoubleArray()
                .count { it <= alphas[i] }
            assertTrue(directRejected > 0)
            assertEquals(directRejected, qvalueRejected)
        }
    }
}
