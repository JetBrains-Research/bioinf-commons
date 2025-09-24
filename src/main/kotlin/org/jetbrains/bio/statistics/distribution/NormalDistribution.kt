package org.jetbrains.bio.statistics.distribution

import org.apache.commons.math3.distribution.AbstractIntegerDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.math3.random.Well19937c
import org.jetbrains.bio.statistics.distribution.NormalIntDistribution.Companion.of
import org.jetbrains.bio.statistics.standardDeviation
import java.util.*
import java.util.stream.DoubleStream
import kotlin.math.ln
import kotlin.math.sqrt

class NormalIntDistribution(
    val mean: Double,
    val variance: Double,
    rng: RandomGenerator = Well19937c(),
) : AbstractIntegerDistribution(rng) {

    init {
        require(variance > 0) {
            "Variance $variance should be > 0"
        }
        require(mean >= 0) {
            "Mean $mean should be >= 0"
        }
    }

    private val backNormalDistribution = NormalDistribution(mean, sqrt(variance))

    override fun probability(k: Int): Double {
        if (k < 0 || k == Integer.MAX_VALUE) {
            return 0.0
        }
        return backNormalDistribution.probability(k.toDouble(), k.toDouble() + 1)
    }

    override fun logProbability(k: Int): Double {
        if (k < 0 || k == Integer.MAX_VALUE) {
            return Double.NEGATIVE_INFINITY
        }
        return ln(backNormalDistribution.probability(k.toDouble(), k.toDouble() + 1))
    }

    override fun cumulativeProbability(k: Int): Double {
        if (k == Integer.MAX_VALUE) {
            return 1.0
        }
        return backNormalDistribution.cumulativeProbability(k.toDouble() + 1)
    }

    override fun getNumericalMean(): Double {
        return mean
    }

    override fun getNumericalVariance(): Double {
        return variance * variance
    }

    override fun getSupportLowerBound(): Int {
        return Integer.MIN_VALUE
    }

    override fun getSupportUpperBound(): Int {
        return Integer.MAX_VALUE
    }

    override fun isSupportConnected(): Boolean {
        return true
    }

    override fun sample(): Int {
        return backNormalDistribution.sample().toInt()
    }

    override fun toString(): String {
        return "NormalDistribution(mean=$mean, sd=$variance)"
    }

    companion object {

        /**
         * Uses defaults from *Constants*.
         * @see [of]
         */
        @Throws(IllegalArgumentException::class)
        fun of(values: DoubleArray): NormalIntDistribution {
            return of(DoubleStream.of(*values).mapToInt { x -> x.toInt() }.toArray())
        }

        /**
         * Creates Negative binomial distribution from MLE of parameters.
         *
         * @param values  sample to estimate parameters from.
         * @return Negative binomial distribution.
         */
        fun of(values: IntArray): NormalIntDistribution {
            val weights = DoubleArray(values.size)
            Arrays.fill(weights, 1.0)

            val mean = values.average()
            val sd = values.standardDeviation()

            return NormalIntDistribution(mean, sd * sd)
        }
    }
}

