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

    override fun link(x: Double): Double { return Math.exp(x) }
    override fun linkDerivative(x: Double): Double { return Math.exp(x) }
    override fun linkVariance(x: Double): Double { return x }
    override fun sampler(x: Double): Int { return Sampling.samplePoisson(x) }

    /**
     * @param t - number of row
     * @param d - number of column, should be a column with observations.
     */
    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        val logLambda = getPredictor(df, t)
        return (
                df.getAsInt(t, df.labels[d])
                        * logLambda
                        - MoreMath.factorialLog(df.getAsInt(t, df.labels[d]))
                        - FastMath.exp(logLambda))

    }
}