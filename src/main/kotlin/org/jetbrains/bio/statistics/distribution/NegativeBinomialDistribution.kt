package org.jetbrains.bio.statistics.distribution

import gnu.trove.map.hash.TIntObjectHashMap
import org.apache.commons.math3.distribution.AbstractIntegerDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.distribution.PoissonDistribution
import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.math3.random.Well19937c
import org.apache.commons.math3.special.Beta
import org.apache.commons.math3.special.Gamma
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.Precision
import org.apache.log4j.Logger
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import org.jetbrains.bio.viktor.asF64Array
import java.util.*
import java.util.stream.DoubleStream
import kotlin.math.abs
import kotlin.math.ln

/**
 * Implementation of a Negative Binomial distribution.
 * Common distributions are available in Apache Commons Math library [org.apache.commons.math3.distribution]
 *
 * @author lebedev
 * @since 17/10/13
 * @see [Negative binomial distribution](http://en.wikipedia.org/wiki/Negative_binomial_distribution)
 */
class NegativeBinomialDistribution(rng: RandomGenerator,
                                   private val mean: Double,
                                   internal val failures: Double) : AbstractIntegerDistribution(rng) {

    /**
     * Access the probability of success for this distribution.
     *
     * @return the probability of success.
     */
    val probabilityOfSuccess: Double

    /* Precalculated constants to speed up probability calculation */
    /**
     * failures * log(1-p) - log Gamma(failures)
     */
    private val fLog1MinusPMinusLogGammaF: Double
    /**
     * (1-p)^failures / Gamma(failures)
     */
    private val oneMinusPToFDivideByGammaF: Double
    private val logP: Double
    private val logMean: Double
    private val expMinusMean: Double

    init {
        require(failures > 0) {
            "Number of failures $failures should be > 0"
        }
        require(mean >= 0) {
            "Mean $mean should be >= 0"
        }
        probabilityOfSuccess = if (failures.isInfinite()) 0.0 else mean / (mean + failures)
        fLog1MinusPMinusLogGammaF = if (failures.isInfinite())
            0.0
        else
            failures * Math.log(1.0 - probabilityOfSuccess) - Gamma.logGamma(failures)
        oneMinusPToFDivideByGammaF = if (failures.isInfinite())
            0.0
        else
            Math.pow(1 - probabilityOfSuccess, failures) / Gamma.gamma(failures)
        logP = Math.log(probabilityOfSuccess)
        logMean = if (mean == 0.0) Double.NEGATIVE_INFINITY else Math.log(mean)
        expMinusMean = Math.exp(-mean)
    }

    override fun probability(k: Int): Double {
        if (k < 0 || k == Integer.MAX_VALUE) {
            return 0.0
        }
        if (mean == 0.0) {
            return if (k == 0) 1.0 else 0.0
        }
        return if (failures.isInfinite()) {
            Math.pow(mean, k.toDouble()) * expMinusMean / CombinatoricsUtils.factorialDouble(k)
        } else Gamma.gamma(k + failures) / CombinatoricsUtils.factorialDouble(k) *
                oneMinusPToFDivideByGammaF * Math.pow(probabilityOfSuccess, k.toDouble())
    }

    override fun logProbability(k: Int): Double {
        if (k < 0 || k == Integer.MAX_VALUE) {
            return Double.NEGATIVE_INFINITY
        }
        if (mean == 0.0) {
            return if (k == 0) 0.0 else Double.NEGATIVE_INFINITY
        }
        return if (failures.isInfinite()) {
            k * logMean - mean - CombinatoricsUtils.factorialLog(k)
        } else Gamma.logGamma(k + failures) - CombinatoricsUtils.factorialLog(k) +
                fLog1MinusPMinusLogGammaF + k * logP
    }

    override fun cumulativeProbability(k: Int): Double {
        if (k < 0) {
            return 0.0
        } else if (k == Integer.MAX_VALUE) {
            return 1.0
        }

        return 1 - Beta.regularizedBeta(probabilityOfSuccess, (k + 1).toDouble(), failures)
    }

    override fun getNumericalMean(): Double {
        return mean
    }

    override fun getNumericalVariance(): Double {
        return mean / (1.0 - probabilityOfSuccess)
    }

    override fun getSupportLowerBound(): Int {
        return 0
    }

    override fun getSupportUpperBound(): Int {
        return Integer.MAX_VALUE
    }

    override fun isSupportConnected(): Boolean {
        return true
    }

    override fun sample(): Int {
        // Note(lebedev): we treat negative binomial distribution as a
        // Gamma-Poisson mixture.
        val lambda = if (failures.isInfinite())
            mean
        else
            GammaDistribution(random, failures, probabilityOfSuccess / (1 - probabilityOfSuccess)).sample()
        return if (Precision.equals(lambda, 0.0)) {
            0
        } else PoissonDistribution(random, lambda,
                PoissonDistribution.DEFAULT_EPSILON,
                PoissonDistribution.DEFAULT_MAX_ITERATIONS).sample()
    }

    override fun toString(): String {
        return "NegativeBinomial(failures r=$failures, probability of success p=$probabilityOfSuccess)"
    }

    companion object {
        private val LOG = Logger.getLogger(NegativeBinomialDistribution::class.java)

        fun usingMean(mean: Double, failures: Double): NegativeBinomialDistribution {
            return NegativeBinomialDistribution(Well19937c(), mean, failures)
        }

        fun usingSuccessProbability(p: Double, failures: Double): NegativeBinomialDistribution {
            return NegativeBinomialDistribution(Well19937c(), p * failures / (1.0 - p), failures)
        }

        /**
         * Uses defaults from *Constants*.
         * @see [of]
         */
        @Throws(IllegalArgumentException::class)
        fun of(values: DoubleArray): NegativeBinomialDistribution {
            return of(DoubleStream.of(*values).mapToInt { x -> x.toInt() }.toArray())
        }

        /**
         * Creates Negative binomial distribution from MLE of parameters.
         *
         * @param values  sample to estimate parameters from.
         * @return Negative binomial distribution.
         */
        fun of(values: IntArray): NegativeBinomialDistribution {
            val weights = DoubleArray(values.size)
            Arrays.fill(weights, 1.0)

            val mean = values.average()
            val sd = values.standardDeviation()

            val failures = fitNumberOfFailures(values,
                    weights.asF64Array(),
                    mean,
                    estimateFailuresUsingMoments(mean, sd * sd))

            return usingMean(mean, failures)
        }

        /** Gamma fit based on
         * http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
         */
        internal fun fitGamma(meanLog: Double, mean: Double, prevA: Double): Double {
            if (mean == 0.0) {
                return Double.NaN
            }
            val logMean = Math.log(mean)
            var aInv = 1 / prevA
            var a = prevA

            for (i in 0..99) {
                aInv += (meanLog - logMean + Math.log(a) - Gamma.digamma(a)) * aInv * aInv / (aInv - Gamma.trigamma(a))

                val aNext = 1 / aInv
                if (Math.abs(a - aNext) < 1E-6 * a) {
                    return aNext
                }
                a = aNext
            }
            return a
        }

        private fun diGammaInPlace(eta: F64Array): F64Array {
            for (i in 0 until eta.size) {
                eta[i] = Gamma.digamma(eta[i])
            }
            return eta
        }

        private fun triGammaInPlace(eta: F64Array): F64Array {
            for (i in 0 until eta.size) {
                eta[i] = Gamma.trigamma(eta[i])
            }
            return eta
        }

        /**
         * NegativeBinomial fit based on
         * http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
         */
        fun fitNumberOfFailures(values: IntArray,
                                weights: F64Array,
                                mean: Double,
                                failures: Double): Double {
            if (failures.isInfinite() || mean == 0.0) {
                /* this is already Poisson or singular distribution, can't be optimized */
                return failures
            }

            val valuesMap = TIntObjectHashMap<KahanSum>()

            for (i in values.indices) {
                val value = values[i]
                val weight = weights[i]
                if (!valuesMap.contains(value)) {
                    valuesMap.put(value, KahanSum())
                }
                valuesMap[value].feed(weight)
            }

            val valuesSmall = valuesMap.keys()
            val weightsSmall = F64Array(valuesMap.size()) { i ->
                valuesMap[valuesSmall[i]].result()
            }

            var a = failures // shape

            val sum = weightsSmall.sum()

            for (i in 0..199) {
                val lambdaExpectationSum = KahanSum()
                val logLambdaExpectationSum = KahanSum()

                for (e in valuesSmall.indices) {

                    val weight = weightsSmall[e]
                    val value = valuesSmall[e]

                    val a2 = a + value     // shape of posterior Gamma
                    val lambdaExpectation = a2
                    val logLambdaExpectation = Gamma.digamma(a2)
                    lambdaExpectationSum.feed(lambdaExpectation * weight)
                    logLambdaExpectationSum.feed(logLambdaExpectation * weight)
                }

                val newMean = lambdaExpectationSum.result() / sum
                val newMeanLog = logLambdaExpectationSum.result() / sum
                val newA = fitGamma(newMeanLog, newMean, a)
                if (newA.isNaN() || abs(a - newA) < a * 1E-6) {
                    return a
                }
                a = newA
            }

            return a
        }

        fun fitNumberOfFailures(values: F64Array,
                                  weights: F64Array,
                                  mean: F64Array,
                                  failures: Double): Double {

            if (failures.isInfinite() || mean.equals(0)) {
                /* this is already Poisson or singular distribution, can't be optimized */
                return failures
            }
            var a = weights.sum()/(weights*(values/mean - 1.0)*(values/mean - 1.0)).sum() // shape
            //var a = 1.0
            for (i in 0..99) {
                a = abs(a)
                val fDeriv = (weights * (diGammaInPlace(values + a) -
                        Gamma.digamma(a) -
                        (mean + a).apply { logInPlace() } -
                        (values + a) / (mean + a) + ln(a) + 1.0))
                        .sum()
                val fSecDeriv = (weights *
                        (triGammaInPlace(values + a) -
                            Gamma.trigamma(a) +
                            ((values + a)/((mean + a).apply { timesAssign(mean + a)})) +
                            (1.0 / a)) -
                        weights * 2.0 / (mean + a))
                        .sum()

                val aNext = a - fDeriv/fSecDeriv
                if (Math.abs(a - aNext) < 1E-4) {
                    println(aNext)
                    return aNext
                }
                a = aNext
            }
            return a
        }

        fun estimateFailuresUsingMoments(mean: Double, variance: Double): Double {
            val p1 = if (Precision.equals(variance, 0.0)) 1.0 else mean / variance
            val p0 = 1 - p1
            if (p0 < 0.0001) {
                LOG.debug("Failures = ${Double.POSITIVE_INFINITY} for negative binomial distribution: " +
                        "mean = $mean is greater than variance = $variance")
                return Double.POSITIVE_INFINITY
            }
            return mean * p1 / p0
        }
    }
}
