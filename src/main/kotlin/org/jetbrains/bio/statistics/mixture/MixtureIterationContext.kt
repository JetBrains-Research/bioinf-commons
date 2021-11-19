package org.jetbrains.bio.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.model.IterationContext
import org.jetbrains.bio.viktor.F64Array

/**
 * Iteration context for both frequentist and Bayesian mixture models.
 *
 * @author Sergei Lebedev
 * @since 04/12/14
 */
abstract class MixtureIterationContext(numStates: Int, df: DataFrame) :
    IterationContext(numStates, df) {

    // Note(lebedev): for Bayesian models all quantities bellow are expected
    // values under the variational posterior distribution.
    val logGammas: F64Array = F64Array(numStates, df.rowsNumber)

    /**
     * Joint log-probability of the observation x_t being generated
     * by component i.
     *
     * \log P(x_t, i) = P(i) P(x_t|i)
     *
     * This field is filled by [refill].
     */
    val logJointProbabilities: F64Array = F64Array(df.rowsNumber, numStates)

    override fun expect() = MixtureInternals.evaluate(logJointProbabilities, logGammas)

    /**
     * Prepares "refilled" fields for the next iteration.
     */
    abstract override fun refill()
}