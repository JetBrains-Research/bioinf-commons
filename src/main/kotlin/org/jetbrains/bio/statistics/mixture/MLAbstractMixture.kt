package org.jetbrains.bio.statistics.mixture

import com.google.common.base.MoreObjects
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.CategoricalDistribution
import org.jetbrains.bio.statistics.forking
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.MLMonitor
import org.jetbrains.bio.statistics.model.SamplingChain
import org.jetbrains.bio.util.MultitaskProgress
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.ln

/**
 * A generic mixture model with parameters estimated via ML.
 *
 * @author Sergei Lebedev
 * @since 06/09/13
 */
abstract class MLAbstractMixture(protected val numComponents: Int, weights: F64Array) : ClassificationModel {
    val logWeights: F64Array = weights.log()

    val weights: F64Array get() = logWeights.exp()

    override fun degreesOfFreedom() = numComponents - 1  // component weights.

    open fun context(df: DataFrame): MixtureIterationContext {
        return MLMixtureIterationContext(df)
    }

    override fun fit(
        preprocessed: Preprocessed<DataFrame>, title: String,
        threshold: Double, maxIter: Int
    ) {
        val df = preprocessed.get()
        val context = context(df)
        val monitor = MLMonitor(title, threshold, maxIter)
        while (true) {
            context.iterate()
            val logLikelihood = MixtureInternals.logLikelihood(df, context.logJointProbabilities)
            MultitaskProgress.reportTask(title)
            if (monitor.monitor(logLikelihood)) {
                monitor.finish(this)
                MultitaskProgress.finishTask(title)
                break
            }

            context.logGammas.expInPlace()
            updateParameters(df, context.logGammas)
        }
    }

    override fun predict(preprocessed: Preprocessed<DataFrame>): IntArray {
        val df = preprocessed.get()
        val context = context(df)
        context.iterate()
        return MixtureInternals.predict(context.logGammas)
    }

    override fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array {
        val df = preprocessed.get()
        val context = context(df)
        context.iterate()
        return context.logGammas
    }

    override fun logLikelihood(preprocessed: Preprocessed<DataFrame>): Double {
        val df = preprocessed.get()
        val context = context(df)
        context.iterate()
        return MixtureInternals.logLikelihood(df, context.logJointProbabilities)
    }

    protected fun samplingChain(numObservations: Int): SamplingChain {
        return SamplingChain.start(sampleStates(numObservations))
    }

    protected fun sampleStates(numObservations: Int): IntArray {
        return CategoricalDistribution(weights.toDoubleArray()).sample(numObservations)
    }

    protected open fun updateParameters(df: DataFrame, gammas: F64Array) {
        for (i in 0 until numComponents) {
            val rowView = gammas.V[i]
            logWeights[i] = ln(rowView.sum()) - ln(rowView.size.toDouble())
        }
    }

    protected abstract fun logProbability(i: Int, df: DataFrame, t: Int): Double

    protected fun toStringHelper(): MoreObjects.ToStringHelper {
        return MoreObjects.toStringHelper(this).add("weights", weights)
    }

    protected open inner class MLMixtureIterationContext(df: DataFrame) :
        MixtureIterationContext(numComponents, df) {

        override fun refill() {
            val numObservations = df.rowsNumber
            (0 until numComponents).forking { i ->
                for (observation in 0 until numObservations) {
                    logJointProbabilities[observation, i] =
                        logProbability(i, df, observation) + logWeights[i]
                }
            }
        }

        override fun expect() = MixtureInternals.evaluate(logJointProbabilities, logGammas)
    }
}
