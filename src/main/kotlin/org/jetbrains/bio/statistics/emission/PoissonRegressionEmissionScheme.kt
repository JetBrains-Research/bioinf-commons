package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.distribution.Sampling

/**
 *
 * Poisson regression.
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
class PoissonRegressionEmissionScheme(
        covariateLabels: List<String>,
        regressionCoefficients: DoubleArray
) : IntegerRegressionEmissionScheme(covariateLabels, regressionCoefficients) {

    override val link: (Double) -> Double = { Math.exp(it) }
    override val linkDerivative: (Double) -> Double = { Math.exp(it) }
    override val linkVariance: (Double) -> Double = { it }
    override val sampler: (Double) -> Int = { Sampling.samplePoisson(it) }

    /**
     * @param t - number of row
     * @param d - number of column, should be a column with observations.
     */
    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        val logLambda = getLogObservation(df, t)
        return (
                df.getAsInt(t, df.labels[d])
                        * logLambda
                        - MoreMath.factorialLog(df.getAsInt(t, df.labels[d]))
                        - FastMath.exp(logLambda))

    }
}