package org.jetbrains.bio.statistics.model

import com.google.common.primitives.Shorts
import org.apache.commons.math3.distribution.AbstractIntegerDistribution
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.distribution.Sampling
import java.util.function.IntToDoubleFunction

/**
 * A chain of data frame to data frame transformations used when building
 * a sampler.
 *
 * @author Sergei Lebedev
 * @since 12/09/14
 */
class SamplingChain private constructor(
        val states: IntArray,
        private val df: DataFrame) {

    private fun wrap(newDf: DataFrame) = SamplingChain(states, newDf)

    fun result() = df.with("state", states)

    fun map(f: (DataFrame) -> DataFrame): SamplingChain = wrap(f(df))

    // Note(lebedev): this is too common not to be factored in a method.
    fun binomial(kField: String, nField: String,
                 nDistribution: AbstractIntegerDistribution,
                 successProbability: IntToDoubleFunction): SamplingChain {
        val numObservations = states.size
        val ns = ShortArray(numObservations)
        val ks = ShortArray(numObservations)
        for (t in 0 until numObservations) {
            // We sample a zero-truncated version of the binomial distribution,
            // because n = 0 doesn't make sense for our models.
            ns[t] = Shorts.checkedCast((nDistribution.sample() + 1).toLong())
            ks[t] = Shorts.checkedCast(Sampling.sampleBinomial(
                    ns[t].toInt(),
                    successProbability.applyAsDouble(states[t])).toLong())
        }

        return wrap(df.with(nField, ns).with(kField, ks))
    }

    companion object {
        @JvmStatic
        fun start(states: IntArray): SamplingChain {
            require(states.isNotEmpty()) { "no data" }
            return SamplingChain(states, DataFrame())
        }
    }
}