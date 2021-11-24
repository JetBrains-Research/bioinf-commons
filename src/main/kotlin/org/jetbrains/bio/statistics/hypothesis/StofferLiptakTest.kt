package org.jetbrains.bio.statistics.hypothesis

import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.jetbrains.bio.viktor.KahanSum
import java.util.*
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/**
 * Stoffer-Liptak test for combining and adjustment dependent pvalues.
 *
 * Unlike the adjustment for the Fisher's combined probability test,
 * the Stouffer-Liptak joint p-value can be applied for nonparametric data.
 *
 * References:
 *
 * Pedersen, Brent S., et al. "Comb-p: software for combining, analyzing,
 *      grouping and correcting spatially correlated P-values." Bioinformatics, 2012.
 *
 * Kechris, Katerina J., Brian Biehs, and Thomas B. Kornberg. "Generalizing moving averages for tiling arrays using
 *      combined p-value statistics." Statistical applications in genetics and molecular biology, 2010.
 */
class StofferLiptakTest(pValues: DoubleArray, maxCorrelationDistance: Int = MAX_CORRELATION_DISTANCE) {
    init {
        require(!pValues.all { it == 0.0 }) { "Zero pvalues array" }
    }

    private val distanceCorrelations = computeCorrelations(pValues, maxCorrelationDistance)

    fun combine(pValues: DoubleArray): Double {
        require(pValues.isNotEmpty()) { "Empty array" }
        val n = pValues.size
        if (n == 1) {
            return pValues.first()
        }
        val correctionSum = KahanSum()
        val zSum = KahanSum()
        for (i in 1 until n) {
            zSum.feed(zscore(pValues[i]))
            for (j in i + 1 until n) {
                val distance = min(distanceCorrelations.size - 1, j - i)
                correctionSum.feed(distanceCorrelations[distance])
            }
        }
        val fullCorrection = n + 2.0 * correctionSum.result()
        check(fullCorrection > 0) {
            "Negative correction value for sqrt $fullCorrection"
        }
        val testStat = zSum.result() / sqrt(fullCorrection)
        check(!testStat.isNaN()) {
            "Nan during combining pvalues ${pValues.joinToString(",") { it.toString() }}"
        }
        // cumulativeProbability may return 1.0 in case of big testStat values
        return max(EPSILON, 1.0 - NORMAL.cumulativeProbability(testStat))
    }

    fun zscore(p: Double): Double {
        // Process the cases where there are many pvalues equal to zero or ones
        var pCorrected = p
        if (p > 1.0 - EPSILON)
            pCorrected = 1.0 - EPSILON
        else if (p < EPSILON)
            pCorrected = EPSILON
        // This shouldn't be too close to 0 or 1, otherwise returns Inf
        return NORMAL.inverseCumulativeProbability(1.0 - pCorrected)
    }

    companion object {
        /**
         * EPSILON - min pvalue threshold, so that [zscore] inverseCumulativeProbability doesn't return +Inf
         */
        internal const val EPSILON: Double = 9e-17
        internal val NORMAL = NormalDistribution()

        /**
         * Maximum distance to compute correlations between p-values for Stoffer-Liptak test
         */
        const val MAX_CORRELATION_DISTANCE = 100

        internal fun computeCorrelations(pValues: DoubleArray, maxCorrelationDistance: Int): DoubleArray {
            val distanceCorrelations = DoubleArray(min(pValues.size / 2, maxCorrelationDistance) + 1)
            val shifted = DoubleArray(pValues.size)
            for (i in 1 until distanceCorrelations.size) {
                System.arraycopy(pValues, i, shifted, 0, pValues.size - i - 1)
                Arrays.fill(shifted, pValues.size - i - 1, shifted.size, 0.0)
                var correlation = 0.0
                if (!pValues.all { it == 0.0 } && !shifted.all { it == 0.0 }) {
                    correlation = PearsonsCorrelation().correlation(pValues, shifted)
                }
                check(!correlation.isNaN() && correlation.isFinite()) {
                    "Wrong correlation between pvalues"
                }
                distanceCorrelations[i] = correlation
            }
            return distanceCorrelations
        }

    }
}
