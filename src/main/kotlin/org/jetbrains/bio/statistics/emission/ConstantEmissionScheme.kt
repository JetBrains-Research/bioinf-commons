package org.jetbrains.bio.statistics.emission

import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import java.util.function.IntSupplier

class ConstantIntegerEmissionScheme(private val emission: Int) : IntegerEmissionScheme {
    override val degreesOfFreedom: Int = 0

    override fun sampler(): IntSupplier = IntSupplier { emission }

    override fun logProbability(value: Int): Double {
        return if (value == emission) 0.0 else Double.NEGATIVE_INFINITY
    }

    override fun update(sample: IntArray, weights: F64Array) {
        /* do nothing */
        LOG.debug("skip")
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(ConstantIntegerEmissionScheme::class.java)
    }
}
