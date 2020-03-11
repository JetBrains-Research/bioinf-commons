package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.special.Gamma
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.gson.NotDirectlyDeserializable
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import java.util.function.IntSupplier

class NegBinEmissionScheme(mean: Double, failures: Double) :
        IntegerEmissionScheme, NotDirectlyDeserializable {

    override val degreesOfFreedom: Int = 2

    var mean: Double = mean
        private set
    var failures: Double = failures
        private set

    @Transient
    private var logMean: Double = 0.0
    @Transient
    private var fLog1MinusPMinusLogGammaF: Double = 0.0
    @Transient
    private var logP: Double = 0.0

    init {
        updateTransients()
    }

    val successProbability: Double get() = mean / (mean + failures)

    override fun sampler(): IntSupplier {
        val distribution =
                NegativeBinomialDistribution(Sampling.RANDOM_DATA_GENERATOR.randomGenerator, mean, failures)
        return IntSupplier { distribution.sample() }
    }

    override fun logProbability(value: Int): Double = when {
        mean.isNaN() || failures.isNaN() || value < 0 || value == Integer.MAX_VALUE -> Double.NEGATIVE_INFINITY
        mean == 0.0 -> if (value == 0) 0.0 else Double.NEGATIVE_INFINITY
        failures.isInfinite() -> value * logMean - mean - MoreMath.factorialLog(value)
        else -> {
            Gamma.logGamma(value + failures) - MoreMath.factorialLog(value) +
                    fLog1MinusPMinusLogGammaF + value * logP
        }
    }

    override fun update(sample: IntArray, weights: F64Array) {
        mean = weights.dot(sample) / weights.sum()
        failures = NegativeBinomialDistribution.fitNumberOfFailures(sample, weights, mean, failures)
        LOG.debug("Mean $mean\tFailures $failures")
        updateTransients()
    }

    override fun updateTransients() {
        logMean = Math.log(mean)
        val p = successProbability
        logP = Math.log(p)
        fLog1MinusPMinusLogGammaF = if (failures.isInfinite()) {
            0.0
        } else {
            failures * Math.log(1.0 - p) - Gamma.logGamma(failures)
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(NegBinEmissionScheme::class.java)
    }
}
