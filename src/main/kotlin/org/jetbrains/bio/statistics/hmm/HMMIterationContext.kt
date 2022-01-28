package org.jetbrains.bio.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.chunked
import org.jetbrains.bio.statistics.model.IterationContext
import org.jetbrains.bio.viktor.F64Array
import java.util.concurrent.ForkJoinTask

/**
 * Iteration context for the both frequentist and Bayesian hidden
 * Markov models.
 *
 * @author Sergei Lebedev
 * @since 09/10/14
 */
abstract class HMMIterationContext(
    numStates: Int,
    val logPriorProbabilities: F64Array,
    val logTransitionProbabilities: F64Array,
    df: DataFrame
) :
    IterationContext(numStates, df) {

    private val numObservations = df.rowsNumber

    // Note(lebedev): for Bayesian models all quantities below are expected
    // values under the variational posterior distribution.
    val logXiSums: F64Array = F64Array(numStates, numStates)
    val logGammas: F64Array = F64Array(numStates, numObservations)
    val logForwardProbabilities: F64Array = F64Array(numObservations, numStates)
    val logBackwardProbabilities: F64Array = F64Array(numObservations, numStates)
    val logObservationProbabilities: F64Array = F64Array(numObservations, numStates)  // filled.

    /**
     * Computes expectations of the latent indicator variables.
     *
     *     \gamma_i(t) = E[I{s_t = i}]
     *     \xi_{ij}(t) = E[I{s_{t - 1} = i} I{s_t = j}]
     */
    override fun expect() {
        ForkJoinTask.invokeAll(
            ForkJoinTask.adapt {
                HMMInternals.logForward(
                    df, numStates, logPriorProbabilities,
                    logTransitionProbabilities,
                    logObservationProbabilities,
                    logForwardProbabilities
                )
            },
            ForkJoinTask.adapt {
                HMMInternals.logBackward(
                    df, numStates, logTransitionProbabilities,
                    logObservationProbabilities,
                    logBackwardProbabilities
                )
            })

        val negativeInfinity = F64Array.full(
            numStates, numStates,
            init = Double.NEGATIVE_INFINITY
        )
        val s = (1 until df.rowsNumber).chunked().map { chunk ->
            val localLogXiSums = negativeInfinity.copy()
            val logXit = negativeInfinity.copy()

            for (s in chunk.lo until chunk.hi) {
                for (priorState in 0 until numStates) {
                    for (nextState in 0 until numStates) {
                        logXit[priorState, nextState] = logForwardProbabilities[s - 1, priorState] +
                                logTransitionProbabilities[priorState, nextState] +
                                logObservationProbabilities[s, nextState] +
                                logBackwardProbabilities[s, nextState]
                    }
                }

                logXit.logRescale()
                localLogXiSums.logAddExpAssign(logXit)
            }

            localLogXiSums
        }

        val initial = negativeInfinity.copy()
        s.reduce(initial) { left, right ->
            val res = left.copy()
            res.logAddExpAssign(right)
            res
        }.copyTo(logXiSums)

        // TODO(lebedev): we can evaluate 'logGammas' simultaneously with
        // 'logXiSums' in the above loop.
        HMMInternals.evaluate(df, logForwardProbabilities, logBackwardProbabilities, logGammas)
    }

    // This is ad-hoc. Please inline.
    fun calculateLogForwardProbabilities(): F64Array {
        refill()
        HMMInternals.logForward(
            df, numStates, logPriorProbabilities, logTransitionProbabilities,
            logObservationProbabilities, logForwardProbabilities
        )
        return logForwardProbabilities
    }
}

