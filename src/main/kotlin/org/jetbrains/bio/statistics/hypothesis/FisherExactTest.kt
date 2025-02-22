package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.statistics.MoreMath
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.max
import kotlin.math.min

/**
 * Hypergeometric distribution describes the number of successes k
 * in n draws from a population of size N with a total of K successes.
 *
 * The corresponding test (also known as Fisher Exact Test) evaluates
 * a null hypothesis of independence. The one-sided P-value can be
 * calculated as follows:
 *
 *              min(k, K)
 *     P(X <= k) = sum    C(K, i) C(N - K, n - i) / C(N, n)
 *          max(0, n + K - N)
 *
 * The two-sided P-value is computed by summing up the probabilities
 * of all the tables with probability less than or equal to that of
 * the observed table.
 */
class FisherExactTest(
    private val N: Int,
    private val K: Int,
    private val n: Int,
    private val k: Int
) {

    /**
     * Support of the Hypergeometric distribution with parameters
     * [N], [K] and [n].
     */
    private val support: IntRange get() = max(0, n - (N - K))..min(n, K)

    init {
        require(N >= 0) { "N must be >= 0 (N = $N)" }
        require(K >= 0 && k <= N) { "K must be from [0; $N] (K = $K)" }
        require(n in 0..N) { "n must be from [0; $N] (n = $n)" }
        require(k in support) {
            "k must be from [${support.first}; ${support.last}] (k = $k)"
        }
    }

    private fun dnhyperNCP(ncp: Double = 1.0): List<Double> {
        val logSupportProbs = support.map { hypergeometricProbabilityLog(N, K, n, it) }
        var d = support.zip(logSupportProbs).map { (sup, prob) ->
            prob + ln(ncp) * sup
        }
        val dMax = d.maxOrNull()!!
        d = d.map { exp(it - dMax) }
        val dSum = d.sum()
        return d.map {it / dSum}

    }

    operator fun invoke(alternative: Alternative = Alternative.LESS): Double {
        return when (alternative) {
            Alternative.LESS -> {
                var acc = 0.0
                for (x in support.first..min(k, K)) {
                    acc += hypergeometricProbability(N, K, n, x)
                }

                acc
            }

            Alternative.TWO_SIDED -> {
                val relativeError = 1 + 1e-7
                val dnhyper = dnhyperNCP()
                val acc = dnhyper.filter {
                    it <= dnhyper[k - support.first] * relativeError
                }.sum()
                if (acc > 1.0) return 1.0 else return acc
            }

            else -> error(alternative.toString())
        }
    }

    companion object {
        /**
         * Construct FET for a given contingency table.
         *
         *   a | b
         *   -----
         *   c | d
         */
        @JvmStatic
        fun forTable(a: Int, b: Int, c: Int, d: Int): FisherExactTest {
            return FisherExactTest(
                N = a + b + c + d, K = a + b,
                n = a + c, k = a
            )
        }
        /**
         * This is up to 30x faster than calling [HypergeometricDistribution#logProbability]
         * because of tabulated factorial
         */
        private fun hypergeometricProbability(N: Int, K: Int, n: Int, k: Int): Double {
            return exp(
                MoreMath.binomialCoefficientLog(K, k)
                        + MoreMath.binomialCoefficientLog(N - K, n - k)
                        - MoreMath.binomialCoefficientLog(N, n)
            )
        }

        private fun hypergeometricProbabilityLog(N: Int, K: Int, n: Int, k: Int): Double {
            return MoreMath.binomialCoefficientLog(K, k) +
                    MoreMath.binomialCoefficientLog(N - K, n - k) -
                    MoreMath.binomialCoefficientLog(N, n)
        }

    }
}