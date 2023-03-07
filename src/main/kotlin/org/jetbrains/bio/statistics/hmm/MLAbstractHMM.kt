package org.jetbrains.bio.statistics.hmm

import com.google.common.base.MoreObjects
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.CategoricalDistribution
import org.jetbrains.bio.statistics.forking
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.MLMonitor
import org.jetbrains.bio.statistics.model.SamplingChain
import org.jetbrains.bio.util.MultitaskProgress
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor._I

/**
 * A generic HMM with parameters estimated via ML.
 *
 * @author Sergei Lebedev
 * @since 15/08/13
 */
abstract class MLAbstractHMM(
    protected val numStates: Int,
    priorProbabilities: F64Array,
    transitionProbabilities: F64Array
) : ClassificationModel {

    init {
        require(numStates > 1) { "expected at least two states" }
    }

    val logPriorProbabilities: F64Array = priorProbabilities.log()
    val logTransitionProbabilities: F64Array = transitionProbabilities.log()

    val priorProbabilities: F64Array
        get() {
            return logPriorProbabilities.exp()
        }

    val transitionProbabilities: F64Array
        get() {
            return logTransitionProbabilities.exp()
        }

    override fun degreesOfFreedom(): Int {
        return numStates - 1 +              // prior,
                (numStates - 1) * numStates  // transitions.
    }

    fun context(df: DataFrame): HMMIterationContext {
        return MLAbstractHMMIterationContext(df)
    }

    /**
     * Fits a HMM to a given sample using Baum-Welch algorithm.
     *
     * @param preprocessed a sample to fit a model to.
     * @param title human-readable title for the sample, e.g. `"chr1"`.
     * @param threshold convergence threshold.
     * @param maxIterations an upper bound on EM iterations.
     */
    override fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String,
        threshold: Double,
        maxIterations: Int
    ) {
        val df = preprocessed.get()
        val context = context(df)
        val monitor = MLMonitor(title, threshold, maxIterations)
        while (true) {
            context.iterate()
            val logLikelihood = HMMInternals.logLikelihood(df, context.logForwardProbabilities)
            MultitaskProgress.reportTask(title)
            if (monitor.monitor(logLikelihood)) {
                monitor.finish(this)
                MultitaskProgress.finishTask(title)
                break
            }

            updateParameters(df, context.logXiSums, context.logGammas)
        }
    }

    override fun fit(
        preprocessed: List<Preprocessed<DataFrame>>, title: String,
        threshold: Double, maxIterations: Int
    ) {
        val dfs = preprocessed.map { it.get() }
        val contexts = dfs.map { context(it) }
        val monitor = MLMonitor(title, threshold, maxIterations)
        while (true) {
            contexts.parallelStream().forEach(HMMIterationContext::iterate)

            var logLikelihood = 0.0
            for ((df, context) in dfs.zip(contexts)) {
                val partialLL = HMMInternals.logLikelihood(df, context.logForwardProbabilities)
                logLikelihood += partialLL
            }
            MultitaskProgress.reportTask(title)
            if (monitor.monitor(logLikelihood)) {
                monitor.finish(this)
                MultitaskProgress.finishTask(title)
                break
            }

            updateParameters(dfs, contexts)
        }
    }

    override fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array {
        val df = preprocessed.get()
        val context = context(df)
        context.iterate()
        return context.logGammas
    }

    override fun predict(preprocessed: Preprocessed<DataFrame>): IntArray {
        val df = preprocessed.get()
        val context = context(df)
        context.iterate()
        return HMMInternals.viterbi(
            df, numStates, context.logPriorProbabilities,
            context.logTransitionProbabilities,
            context.logObservationProbabilities
        )
    }

    override fun logLikelihood(preprocessed: Preprocessed<DataFrame>): Double {
        return logLikelihoodInPlace(preprocessed)
    }

    // An less allocating version of log-likelihood suitable for use in
    // Bayesian likelihood computations. See issue #191 on GitHub.
    private fun logLikelihoodInPlace(preprocessed: Preprocessed<DataFrame>): Double {
        val df = preprocessed.get()
        val logForwardProbabilities = F64Array(2, numStates)
        for (state in 0 until numStates) {
            logForwardProbabilities[0, state] = logPriorProbabilities[state] + logProbability(state, df, 0)
        }

        val numObservations = df.rowsNumber
        val workBuffer = F64Array(numStates)
        for (observation in 1 until numObservations) {
            for (nextState in 0 until numStates) {
                for (priorState in 0 until numStates) {
                    workBuffer[priorState] =
                        logForwardProbabilities[(observation - 1) % 2, priorState] +
                                logTransitionProbabilities[priorState, nextState]
                }

                logForwardProbabilities[observation % 2, nextState] =
                    workBuffer.logSumExp() + logProbability(nextState, df, observation)
            }
        }

        return logForwardProbabilities.V[(numObservations - 1) % 2].logSumExp()
    }

    protected fun samplingChain(numObservations: Int): SamplingChain {
        return SamplingChain.start(sampleStates(numObservations))
    }

    protected fun sampleStates(numObservations: Int): IntArray {
        val distributions = Array(numStates) { state ->
            CategoricalDistribution(transitionProbabilities.V[state])
        }

        var state = CategoricalDistribution(priorProbabilities).sample()
        val states = IntArray(numObservations)
        states[0] = state
        for (t in 1 until numObservations) {
            state = distributions[state].sample()
            states[t] = state
        }
        return states
    }

    protected abstract fun logProbability(state: Int, df: DataFrame, observation: Int): Double

    private fun updateParameters(df: DataFrame, logXiSums: F64Array, logGammas: F64Array) {
        for (i in 0 until numStates) {
            logTransitionProbabilities.V[i] = logXiSums.V[i]
            logTransitionProbabilities.V[i].logRescale()
        }

        logPriorProbabilities.V[_I] = logGammas.V[_I, 0]
        logPriorProbabilities.logRescale()

        logGammas.expInPlace()
        updateParameters(df, logGammas)
    }

    protected abstract fun updateParameters(df: DataFrame, gammas: F64Array)


    private fun updateParameters(dfs: List<DataFrame>, contexts: List<HMMIterationContext>) {
        val logXiSums = F64Array.full(numStates, numStates, init = Double.NEGATIVE_INFINITY)

        for (context in contexts) {
            logXiSums.logAddExpAssign(context.logXiSums)
        }

        for (state in 0 until numStates) {
            logPriorProbabilities[state] = java.lang.Double.NEGATIVE_INFINITY
        }

        for (context in contexts) {
            for (state in 0 until numStates) {
                logPriorProbabilities[state] = MoreMath.logAddExp(
                    logPriorProbabilities[state], context.logGammas[state, 0]
                )
            }
        }

        logPriorProbabilities.logRescale()

        for (state in 0 until numStates) {
            logTransitionProbabilities.V[state] = logXiSums.V[state]
            logTransitionProbabilities.V[state].logRescale()
        }

        val df = DataFrame.rowBind(dfs.toTypedArray())
        val numObservations = contexts.sumOf { it.logGammas.shape[1] }
        val logGammas = F64Array(numStates, numObservations)

        var observation = 0
        for (context in contexts) {
            val localLogGammas = context.logGammas
            for (localObservation in 0 until localLogGammas.shape[1]) {
                for (state in 0 until numStates) {
                    logGammas[state, observation] = localLogGammas[state, localObservation]
                }
                observation++
            }
        }
        assert(observation == numObservations) { "Incomplete iterations contexts" }

        logGammas.expInPlace()
        updateParameters(df, logGammas)
    }

    fun sampleStates(preprocessed: Preprocessed<DataFrame>): IntArray {
        val df = preprocessed.get()
        val context = context(df)
        val logForwardProbabilities = context.calculateLogForwardProbabilities()
        val rowsNumber = logForwardProbabilities.shape[0]
        val workBuffer = logForwardProbabilities.V[rowsNumber - 1]
        workBuffer.logRescale()
        var state = CategoricalDistribution(workBuffer.exp()).sample()
        val states = IntArray(rowsNumber)
        states[rowsNumber - 1] = state

        for (t in rowsNumber - 2 downTo 0) {
            workBuffer.V[_I] = logForwardProbabilities.V[t] +
                    logTransitionProbabilities.V[_I, state]
            workBuffer.logRescale()
            state = CategoricalDistribution(workBuffer.exp()).sample()
            states[t] = state
        }

        return states
    }

    protected fun toStringHelper(): MoreObjects.ToStringHelper {
        return MoreObjects.toStringHelper(this)
            .add("priorProbabilities", priorProbabilities)
            .add("transitionProbabilities", transitionProbabilities)
    }

    inner class MLAbstractHMMIterationContext(df: DataFrame) :
        HMMIterationContext(
            numStates, logPriorProbabilities,
            logTransitionProbabilities, df
        ) {

        override fun refill() {
            (0 until numStates).forking { state ->
                val numObservations = df.rowsNumber
                for (observation in 0 until numObservations) {
                    logObservationProbabilities[observation, state] = logProbability(state, df, observation)
                }
            }
        }
    }
}
