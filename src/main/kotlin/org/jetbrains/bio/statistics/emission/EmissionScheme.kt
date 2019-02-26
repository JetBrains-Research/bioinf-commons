package org.jetbrains.bio.statistics.emission

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate

/**
 * @author Alexey Dievsky
 * @date 4/20/15
 */
interface EmissionScheme {
    val degreesOfFreedom: Int

    /**
     * Sample multiple observations into a sample column into rows
     * specified by a predicate.
     *
     * @param df sample
     * @param d dimension number (column)
     * @param fill a predicate which returns `true` for row numbers that
     *             are to be sampled into.
     */
    fun sample(df: DataFrame, d: Int, fill: IntPredicate)

    /**
     * Calculate log P(X=x), where X is distributed according to the scheme,
     * and x is an observation in row `t` and column `d` of the sample.
     *
     * @param df sample
     * @param t observation number (row)
     * @param d dimension number (column)
     * @return log P(X=x)
     */
    fun logProbability(df: DataFrame, t: Int, d: Int): Double

    /**
     * Replace the distribution with MLE from the same family based on weighted
     * observations.
     *
     * @param df sample
     * @param d dimension number (column)
     * @param weights weights vector
     */
    fun update(df: DataFrame, d: Int, weights: F64Array)
}