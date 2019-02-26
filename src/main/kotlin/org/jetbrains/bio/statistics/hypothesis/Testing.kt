package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.argSort

enum class Alternative {
    LESS, GREATER, TWO_SIDED
}

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
class FisherExactTest(private val N: Int,
                      private val K: Int,
                      private val n: Int,
                      private val k: Int) {

    /**
     * Support of the Hypergeometric distribution with parameters
     * [N], [K] and [n].
     */
    private val support: IntRange get() = Math.max(0, n - (N - K))..Math.min(n, K)

    init {
        require(N >= 0) { "N must be >= 0 (N = $N)" }
        require(K >= 0 && k <= N) { "K must be from [0; $N] (K = $K)" }
        require(n in 0..N) { "n must be from [0; $N] (n = $n)" }
        require(k in support) {
            "k must be from [${support.first}; ${support.last}] (k = $k)"
        }
    }

    operator fun invoke(alternative: Alternative = Alternative.LESS): Double {
        return when (alternative) {
            Alternative.LESS -> {
                var acc = 0.0
                for (x in support.start..Math.min(k, K)) {
                    acc += hypergeometricProbability(N, K, n, x)
                }

                acc
            }
            Alternative.TWO_SIDED -> {
                val baseline = hypergeometricProbability(N, K, n, k)
                var acc = 0.0
                for (x in support) {
                    val current = hypergeometricProbability(N, K, n, x)
                    // XXX the cutoff is chosen to be consistent with
                    //     'fisher.test' in R.
                    if (current <= baseline + 1e-6) {
                        acc += current
                    }
                }

                acc
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
        @JvmStatic fun forTable(a: Int, b: Int, c: Int, d: Int): FisherExactTest {
            return FisherExactTest(N = a + b + c + d, K = a + b,
                    n = a + c, k = a)
        }
    }
}

// This is up to 30x faster then calling
// [HypergeometricDistribution#logProbability] because of the
// tabulated factorial.
private fun hypergeometricProbability(N: Int, K: Int, n: Int, k: Int): Double {
    return Math.exp(MoreMath.binomialCoefficientLog(K, k)
                    + MoreMath.binomialCoefficientLog(N - K, n - k)
                    - MoreMath.binomialCoefficientLog(N, n))
}

object Multiple {
    /**
     * Applies Benjamini-Hochberg correction to the P-values.
     *
     * See https://en.wikipedia.org/wiki/False_discovery_rate.
     */
    fun adjust(ps: F64Array): F64Array {
        val m = ps.size
        val sorted = ps.argSort(reverse = true)
        val original = IntArray(m)
        for (k in 0 until m) {
            original[sorted[k]] = k
        }

        val adjusted = F64Array(m)
        for (k in 0 until m) {
            adjusted[k] = Math.min(1.0, ps[sorted[k]] * m / (m - k))
        }

        for (k in 1 until m) {
            adjusted[k] = Math.min(adjusted[k], adjusted[k - 1])
        }

        adjusted.reorder(original)
        return adjusted
    }
}