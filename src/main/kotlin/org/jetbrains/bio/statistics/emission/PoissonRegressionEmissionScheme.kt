package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.viktor.F64Array
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
        regressionCoefficients: DoubleArray
) : IntegerRegressionEmissionScheme(covariateLabels, regressionCoefficients) {

    override fun mean(eta: Double) = exp(eta)
    override fun meanDerivative(eta: Double) = exp(eta)
    override fun meanVariance(mean: Double) = mean

    override fun sampler(mean: Double) = Sampling.samplePoisson(mean)

    override fun meanInPlace(eta: F64Array) = eta.apply { expInPlace() }
    override fun meanDerivativeInPlace(eta: F64Array) = eta.apply { expInPlace() }
    override fun meanVarianceInPlace(mean: F64Array) = mean

    override fun zW(y: F64Array, eta: F64Array): Pair<F64Array, F64Array> {
        // Since h(η) = h'(η) = var(h(η)), we can skip h'(η) and var(h(η)) calculations and simplify W:
        // W = diag(h'(η)^2 / var(h(η))) = h(η)
        val countedLink = meanInPlace(eta.copy())
        eta += (y - countedLink).apply { divAssign(countedLink) }
        return eta to countedLink
    }

    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        // We don't use the existing Poisson log probability because that saves us one logarithm.
        // We would have to provide lambda = exp(logLambda), and the Poisson implementation would then have to
        // calculate log(lambda) again.
        val logLambda = getPredictor(df, t)
        val y = df.getAsInt(t, df.labels[d])
        return y * logLambda - MoreMath.factorialLog(y) - FastMath.exp(logLambda)
    }
}