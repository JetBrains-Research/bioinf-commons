package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.distribution.Sampling
import java.util.function.IntPredicate
import kotlin.math.exp

/**
 *
 * Poisson regression.
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
class PoissonRegressionEmissionScheme(
        covariateLabels: List<String>,
        regressionCoefficients: DoubleArray,
        override val degreesOfFreedom: Int) : IntegerRegressionEmissionScheme(covariateLabels, regressionCoefficients) {

    override val link: (Double) -> Double = { Math.exp(it) }
    override val linkDerivative: (Double) -> Double = { Math.exp(it) }
    override val linkVariance: (Double) -> Double = { it }
    var sample = emptyArray<Double>()
        private set


    /**
     * @param d - number of column in which we want to sample
     * @param fill - predicate which marks rows we need to use for sampling
     */
    override fun sample(df: DataFrame, d: Int, fill: IntPredicate) {

        val observations = df.sliceAsInt(df.labels[d])
        val observation = DoubleArray(covariateLabels.size + 1) { 1.0 }
        (0 until df.rowsNumber).forEach { row ->
            if (fill.test(row)) {
                this.covariateLabels.forEachIndexed { index, label ->
                    observation[index + 1] = df.getAsDouble(row, label)
                }
                observations[row] = Sampling
                        .samplePoisson(exp(regressionCoefficients.zip(observation) { a, b -> a * b }.sum()))
            }
        }
    }

    /**
     * @param t - number of row
     * @param d - number of column, should be a column with observations.
     */
    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        val logLambda = regressionCoefficients
                .zip(
                        doubleArrayOf(1.0)
                                .plus(
                                        df.rowAsDouble(t)
                                                .filterIndexed
                                                { index, _ -> covariateLabels.contains(df.labels[index]) }))
                { a, b -> a * b }
                .sum()
        return (
                df.getAsInt(t, df.labels[d])
                        * logLambda
                        - MoreMath.factorialLog(df.getAsInt(t, df.labels[d]))
                        - FastMath.exp(logLambda))

    }
}