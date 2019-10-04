package org.jetbrains.bio.statistics.distribution

import com.google.common.primitives.Ints
import gnu.trove.set.hash.TIntHashSet
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.PoissonDistribution
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.random.Well19937c
import org.apache.commons.math3.special.Beta
import org.apache.commons.math3.special.Gamma
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.viktor.F64Array

/**
 * Sampling procedures for common statistical distributions.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 18/03/13
 */
object Sampling {
    val RANDOM_DATA_GENERATOR: RandomDataGenerator = RandomDataGenerator()

    fun sampleBeta(alpha: Double, beta: Double): Double {
        val x = sampleGamma(alpha, 1.0)
        val y = sampleGamma(beta, 1.0)
        return x / (x + y)
    }

    fun sampleGamma(shape: Double, rate: Double): Double {
        if (shape == 0.0) {
            return 0.0  // Special case.
        }

        return RANDOM_DATA_GENERATOR.nextGamma(shape, 1.0 / rate)
    }

    fun sampleShiftedGeometric(shift: Int, p: Double): Int {
        require(p > 0.0 && p <= 1.0) { "p" }
        require(shift >= 0) { "shift must be >=0" }
        var k = 0
        while (!sampleBernoulli(p)) {
            k++
        }

        return k + shift
    }

    fun sampleBinomial(n: Int, p: Double): Int {
        return RANDOM_DATA_GENERATOR.nextBinomial(n, p)
    }


    fun samplePoisson(rate: Double): Int {
        require(rate >= 0) { "rate must be >=0" }
        return if (Precision.equals(rate, 0.0)) {
            0
        } else {
            Ints.checkedCast(RANDOM_DATA_GENERATOR.nextPoisson(rate))
        }
    }

    fun sampleNegBinomial(mean: Double, failures: Double): Int {
        require(mean >= 0 && failures > 0) {"mean and failures must be > 0"}
        if (mean == 0.0) return 0
        val rate = sampleGamma(failures, failures / mean)
        return samplePoisson(rate)
    }

    @JvmStatic
    fun sampleUniform(a: Int = 0, b: Int = Int.MAX_VALUE): Int {
        return RANDOM_DATA_GENERATOR.nextInt(a, b)
    }

    @JvmStatic
    fun sampleUniform(a: Double = 0.0, b: Double = 1.0): Double {
        return RANDOM_DATA_GENERATOR.nextUniform(a, b, true)
    }

    /**
     * Samples a Dirichlet distributed random variate.
     *
     * @param concentration concentration parameters (an array).
     * @return an array of weights which sum to one.
     */
    fun sampleDirichlet(concentration: F64Array): F64Array {
        val k = concentration.size
        val weights = F64Array(k)
        for (component in 0 until k) {
            weights[component] = sampleGamma(concentration[component], 1.0)
        }

        weights.rescale()
        return weights
    }

    /**
     * Samples a uniformly random combination of size k from n objects.
     * The implementation follow Floyd algorithm, described in: Bentley, Floyd.,
     * "Programming pearls: a sample of brilliance.",
     * Communications of the ACM 30.9 (1987): 754-757.
     *
     * @param n number of elements available.
     * @param k number of elements to sample.
     * @return sampled combination as k-sized array of indices [0..n-1].
     */
    fun sampleCombination(n: Int, k: Int): IntArray {
        val indicesSet = TIntHashSet()

        for (i in n - k until n) {
            val pos = RANDOM_DATA_GENERATOR.nextInt(0, i)
            if (pos in indicesSet) {
                indicesSet.add(i)
            } else {
                indicesSet.add(pos)
            }
        }

        val indices = indicesSet.toArray()
        indices.sort()
        return indices
    }

    /**
     * Samples a Bernoulli trial variable; returns `true` with probability `probability`
     * and `false` otherwise.
     *
     * @param probability The trial success probability.
     * @return `true` with probability `probability` and `false` otherwise.
     */
    @JvmStatic
    fun sampleBernoulli(probability: Double): Boolean = when {
        Precision.equals(probability, 0.0) -> false
        Precision.equals(probability, 1.0) -> true
        else -> sampleUniform(0.0, 1.0) < probability
    }

}

/**
 * Methods for computing density and p.m.f. values of common distributions.
 *
 * @author Sergei Lebedev
 * @since 25/09/13
 */
object Densities {
    fun logDirichletDensity(weights: F64Array,
                            concentration: F64Array): Double {
        val k = weights.size
        val totalConcentration = concentration.sum()
        var acc = Gamma.logGamma(totalConcentration)
        for (component in 0 until k) {
            val alpha = concentration[component]
            acc += (alpha - 1) * Math.log(weights[component]) - Gamma.logGamma(alpha)
        }

        return acc
    }

    fun logBetaDensity(x: Double, alpha: Double, beta: Double): Double {
        return BetaDistribution(null, alpha, beta).logDensity(x)
    }


    fun logPoissonDensity(n: Int, rate: Double): Double {
        if (rate == 0.0) {
            return if (n == 0) 0.0 else Double.NEGATIVE_INFINITY
        }

        return PoissonDistribution(null, rate,
                PoissonDistribution.DEFAULT_EPSILON,
                PoissonDistribution.DEFAULT_MAX_ITERATIONS)
                .logProbability(n)
    }

    fun logBinomialDensity(k: Int, n: Int, p: Double): Double {
        return BinomialDistribution(null, n, p).logProbability(k)
    }

    fun logNegativeBinomialDensity(n: Int, mean: Double, failures: Double): Double {
        return if (failures.isInfinite()) {
            logPoissonDensity(n, mean)
        } else {
            NegativeBinomialDistribution(Well19937c(), mean, failures).logProbability(n)
        }
    }
}

/**
 * Functions for computing `E[X]` and `E[log X]` for different random variables.
 *
 * @see http://en.wikipedia.org/wiki/Expected_value
 * @author Sergei Lebedev
 * @since 27/02/14
 */
object Expectations {

    fun logDirichlet(concentration: F64Array, dst: F64Array) {
        val totalConcentration = concentration.sum()
        requireNotNaN("concentration", totalConcentration)
        val totalConcentrationDigamma = Gamma.digamma(totalConcentration)
        for (i in 0 until concentration.size) {
            dst[i] = Gamma.digamma(concentration[i]) - totalConcentrationDigamma
        }
    }

    fun beta(alpha: Double, beta: Double): Double = alpha / (alpha + beta)

    fun logBeta(alpha: Double, beta: Double): Double {
        requireNotNaN("alpha", alpha)
        requireNotNaN("beta", beta)
        return Gamma.digamma(alpha) - Gamma.digamma(alpha + beta)
    }
}

private fun requireNotNaN(role: String, x: Double) {
    require(!x.isNaN()) { "$role is not a number" }
}

/**
 * Functions for computing Kullback-Leibler divergence for common
 * statistical distributions.
 *
 * @author Sergei Lebedev
 * @since 20/08/14
 */
object KullbackLeibler {
    fun dirichlet(concentration1: F64Array, concentration2: F64Array): Double {
        val n = concentration1.size
        require(concentration2.size == n)
        val meanLogWeights = F64Array(n)
        Expectations.logDirichlet(concentration2, meanLogWeights)

        var acc = -Gamma.logGamma(concentration1.sum()) + Gamma.logGamma(concentration2.sum())
        for (component in 0 until n) {
            acc -= (concentration1[component] - 1) * meanLogWeights[component] -
                    Gamma.logGamma(concentration1[component])
            acc += (concentration2[component] - 1) * meanLogWeights[component] -
                    Gamma.logGamma(concentration2[component])
        }

        return acc
    }

    fun beta(alpha1: Double, beta1: Double, alpha2: Double, beta2: Double): Double {
        val meanLogSuccessProbability = Expectations.logBeta(alpha2, beta2)
        val meanLogFailureProbability = Expectations.logBeta(beta2, alpha2)
        val logNum = (alpha1 - 1) * meanLogSuccessProbability + (beta1 - 1) *
                meanLogFailureProbability - Beta.logBeta(alpha1, beta1)
        val logDen = (alpha2 - 1) * meanLogSuccessProbability + (beta2 - 1) *
                meanLogFailureProbability - Beta.logBeta(alpha2, beta2)
        return -logNum + logDen
    }

}
