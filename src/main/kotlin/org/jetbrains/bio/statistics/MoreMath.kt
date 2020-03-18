package org.jetbrains.bio.statistics

import org.apache.commons.math3.exception.NotPositiveException
import org.apache.commons.math3.special.Gamma
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.viktor.logAddExp
import kotlin.math.*

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
            max(x, y) + StrictMath.log1p(-FastMath.exp(-abs(x - y)))
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
            val reference = min((i - 1024) / 16, 1023)
            var res = LOG_FACTORIAL_16[reference]
            var j = reference * 16 + 1025
            while (j < i) {
                // the unroll here reduces the number of logarithms by half,
                // while the floating-point unit guards against overflow.
                res += ln(j.toDouble() * (j.toDouble() + 1.0))
                j += 2
            }
            if (j == i) {
                res += ln(i.toDouble())
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
 * Corrected standard deviation from sample of ints
 */
fun IntArray.standardDeviation(): Double {
    var sum = 0L
    var sumSq = 0L
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return sqrt((sumSq - sum.toDouble() * sum.toDouble() / size) / (size - 1))
}

/**
 * Corrected standard deviation from sample of shorts
 */
fun ShortArray.standardDeviation(): Double {
    var sum = 0L
    var sumSq = 0L
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return sqrt((sumSq - sum.toDouble() * sum.toDouble() / size) / (size - 1))
}

/**
 * Corrected standard deviation from sample of longs
 */
fun LongArray.standardDeviation(): Double {
    var sum = 0L
    var sumSq = 0L
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return sqrt((sumSq - sum.toDouble() * sum.toDouble() / size) / (size - 1))
}

