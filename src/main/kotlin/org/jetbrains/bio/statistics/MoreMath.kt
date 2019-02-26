package org.jetbrains.bio.statistics

import org.apache.commons.math3.exception.NotPositiveException
import org.apache.commons.math3.special.Gamma
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.viktor.asF64Array
import org.jetbrains.bio.viktor.logAddExp

/**
 * Useful mathematical routines absent in [java.util.Math]
 * and [org.apache.commons.math3.util.FastMath].
 *
 * When adding new functionality please consider reading
 * http://blog.juma.me.uk/2011/02/23/performance-of-fastmath-from-commons-math.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @since 26/03/15
 */
object MoreMath {
    /**
     * Evaluates log(exp(a) + exp(b)) using the following trick
     *
     *     log(exp(a) + log(exp(b)) = a + log(1 + exp(b - a))
     *
     * assuming a >= b.
     */
    @JvmStatic fun logAddExp(a: Double, b: Double) = a logAddExp b

    /**
     * Computes log |exp(x) - exp(y)| faster and more precise than
     * direct computation, safe from overflow.
     */
    @JvmStatic fun logAbsDiffExp(x: Double, y: Double): Double {
        return if (x.isInfinite() && y.isInfinite()) {
            /* both corresponding values are zero, as is the difference */
            Double.NEGATIVE_INFINITY
        } else {
            Math.max(x, y) + StrictMath.log1p(-FastMath.exp(-Math.abs(x - y)))
        }
    }

    private val LOG_FACTORIAL_1 = DoubleArray(1024)
    private val LOG_FACTORIAL_16 = DoubleArray(1024)
    init {
        for (i in 0..1023) {
            LOG_FACTORIAL_1[i] = CombinatoricsUtils.factorialLog(i);
            LOG_FACTORIAL_16[i] = CombinatoricsUtils.factorialLog(1024 + 16 * i);
        }
    }

    /**
     * Computes log(i!).
     *
     * Uses 16KB of cached values to drastically speed up calculations
     * for most adequate arguments; delegates to [Gamma.logGamma] for
     * others.
     */
    @JvmStatic fun factorialLog(i: Int): Double {
        if (i < 0) {
            throw NotPositiveException(i)
        }
        if (i < 1024) {
            return LOG_FACTORIAL_1[i]
        }
        if (i < 17408) {
            val reference = Math.min((i - 1024) / 16, 1023)
            var res = LOG_FACTORIAL_16[reference]
            var j = reference * 16 + 1025
            while (j < i) {
                // the unroll here reduces the number of logarithms by half,
                // while the floating-point unit guards against overflow.
                res += Math.log(j.toDouble() * (j.toDouble() + 1.0))
                j += 2
            }
            if (j == i) {
                res += Math.log(i.toDouble())
            }
            return res
        }

        return Gamma.logGamma((i + 1).toDouble())
    }

    /**
     * Computes C(n, k) using tabulated [factorialLog].
     */
    @JvmStatic fun binomialCoefficientLog(n: Int, k: Int): Double {
        return factorialLog(n) - factorialLog(k) - factorialLog(n - k)
    }
}


/**
 * A simple static data structure for range sum queries.
 *
 * See http://en.wikipedia.org/wiki/Range_Queries
 *
 * @author Sergei Lebedev
 * @since 10/04/14
 */
class PrefixSumTable(values: DoubleArray) {
    /**
     * Don't use [Arrays#parallelPrefix] here, because [cumSum] uses [KahanSum]
     */
    private val sums = values.clone().asF64Array().apply { cumSum() }

    /**
     * Computes a prefix sum of values up to *right*, i. e.
     *
     *     values[0] + values[1] + ... + values[right]
     *
     * in O(1) time, where n is the number of elements in the initializer
     * array.
     *
     * @param right rightmost index.
     * @return computed prefix sum.
     */
    operator fun get(right: Int) = sums[right]

    /**
     * Computes an interval sum of *[left, right]*, i. e.
     *
     *     values[left] + values[left + 1] + ... + values[right]
     *
     * in O(1) time, where n is the number of elements in the initializer
     * array.
     *
     * @param left leftmost index
     * @param right rightmost index
     * @return computed interval sum.
     */
    operator fun get(left: Int, right: Int): Double {
        return sums[right] - (if (left == 0) 0.0 else sums[left - 1])
    }
}

fun IntArray.standardDeviation(): Double {
    var sum = 0L
    var sumSq = 0L
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return Math.sqrt((sumSq - sum.toDouble() * sum.toDouble() / size) / (size - 1))
}

fun ShortArray.standardDeviation(): Double {
    var sum = 0L
    var sumSq = 0L
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return Math.sqrt((sumSq - sum.toDouble() * sum.toDouble() / size) / (size - 1))
}

