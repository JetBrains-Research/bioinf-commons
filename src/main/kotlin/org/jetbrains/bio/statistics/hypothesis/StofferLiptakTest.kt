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
        // Test assumes positive correlation between p-values
        val testStat = zSum.result() / sqrt(n + 2.0 * max(0.0, correctionSum.result()))
        check(!testStat.isNaN()) {
            "Nan during combining pvalues sum(z)=${zSum.result()} sum(corr)=${correctionSum.result()} " +
                    "pvalues=${pValues.joinToString(",") { it.toString() }}"
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
            for (i in 1 until distanceCorrelations.size) {
                val original = DoubleArray(pValues.size - i - 1)
                System.arraycopy(pValues, 0, original, 0, pValues.size - i - 1)
                val shifted = DoubleArray(pValues.size - i - 1)
                System.arraycopy(pValues, i, shifted, 0, pValues.size - i - 1)
                Arrays.fill(shifted, pValues.size - i - 1, shifted.size, 0.0)
                var correlation = PearsonsCorrelation().correlation(original, shifted)
                // Safeguard against NaN from Pearson correlation for zero or small vectors,
                // See example at: https://github.com/JetBrains-Research/span/issues/34
                if (correlation.isNaN() || !correlation.isFinite()) {
                   correlation = 0.0
                }
                distanceCorrelations[i] = correlation
            }
            return distanceCorrelations
        }
    }
}
