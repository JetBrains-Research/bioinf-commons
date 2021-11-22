package org.jetbrains.bio.statistics.hypothesis

import com.google.common.primitives.Ints
import org.jetbrains.bio.viktor.F64Array
import org.junit.Assert.assertArrayEquals
import org.junit.Assert.assertEquals
import org.junit.Before
import org.junit.Test
import java.util.*
import kotlin.math.max

class FisherExactTestTest {
    @Before
    fun setUp() {
        // Keep the seed fixed, because otherwise the precomputed R
        // P-values will change.
        RANDOM.setSeed(42)
    }

    @Test
    fun againstROneSided() {
        val (Ns, Ks, ns, ks) = generate()

        // Computed by R 3.2.2
        // val arg = "t(matrix(c($k, ${K - k}, ${n - k}, ${N - K - (n - k)}), nrow=2))"
        // "fisher.test($arg, alternative = \"less\", simulate.p.value=F)\$p.value"
        val expected = doubleArrayOf(
            0.9103448, 0.9806255,
            1.0, 1.0, 1.0, 0.5757576,
            1.0, 0.3333333, 1.551382e-05,
            1.0
        )

        for (i in Ns.indices) {
            val N = Ns[i]
            val K = Ks[i]
            val n = ns[i]
            val k = ks[i]
            assertEquals(expected[i], FisherExactTest(N, K, n, k)(), 1e-6)
        }
    }

    @Test
    fun againstRTwoSided() {
        val (Ns, Ks, ns, ks) = generate()

        // Computed by R 3.2.2
        val expected = doubleArrayOf(
            0.5862069, 0.1266209,
            1.0, 1.0, 1.0, 1.0,
            1.0, 0.3333333, 2.064322e-05,
            1.0
        )

        for (i in Ns.indices) {
            val N = Ns[i]
            val K = Ks[i]
            val n = ns[i]
            val k = ks[i]
            assertEquals(
                expected[i],
                FisherExactTest(N, K, n, k)(Alternative.TWO_SIDED),
                1e-6
            )
        }
    }

    @Test
    fun twoSidedCutoff() {
        assertEquals(1.0, FisherExactTest(13, 5, 4, 2)(Alternative.TWO_SIDED), 1e-6)
    }

    private fun generate(): List<IntArray> {
        val Ns = RANDOM.ints(8, 36).limit(10).toArray()
        val Ks = (0 until 10).map { RANDOM.nextInt(Ns[it] + 1) }.toIntArray()
        val ns = (0 until 10).map { RANDOM.nextInt(Ns[it] + 1) }.toIntArray()
        val ks = (0 until 10).map {
            val lower = max(0, ns[it] - (Ns[it] - Ks[it]))
            val upper = Ints.min(ns[it], Ks[it])
            lower + RANDOM.nextInt(upper - lower + 1)
        }.toIntArray()
        return listOf(Ns, Ks, ns, ks)
    }

    companion object {
        private val RANDOM = Random()
    }
}

