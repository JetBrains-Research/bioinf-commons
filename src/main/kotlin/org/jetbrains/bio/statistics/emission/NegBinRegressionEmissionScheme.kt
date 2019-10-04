package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.special.Gamma
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.ln

/**
 *
 * Negative Binomial regression.
 *
 * @author Elena Kartysheva
 * @date 9/13/19
 */
class NegBinRegressionEmissionScheme(
        covariateLabels: List<String>,
        regressionCoefficients: DoubleArray,
        failures: Double
) : IntegerRegressionEmissionScheme(covariateLabels, regressionCoefficients) {

    var failures = failures
        private set

    private var fLogf = failures * ln(failures)

    override fun mean(eta: Double) = exp(eta)
    override fun meanDerivative(eta: Double) = exp(eta)
    override fun meanVariance(mean: Double) = mean + mean * mean / failures

    override fun sampler(mean: Double) = Sampling.sampleNegBinomial(mean, failures)

    override fun meanInPlace(eta: F64Array) = eta.apply { expInPlace() }
    override fun meanDerivativeInPlace(eta: F64Array) = eta.apply { expInPlace() }
    override fun meanVarianceInPlace(mean: F64Array) = mean + mean * mean / failures

    override fun zW(y: F64Array, eta: F64Array): Pair<F64Array, F64Array> {
        val countedLink = meanInPlace(eta.copy())
        eta += (y - countedLink).apply { divAssign(countedLink) }
        return eta to countedLink
    }

    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        val mean = getPredictor(df, t)
        val y = df.getAsInt(t, df.labels[d])
        val logMeanPlusFailure = ln(mean + failures)

        return when {
            failures.isNaN() || y < 0 || y == Integer.MAX_VALUE -> Double.NEGATIVE_INFINITY
            mean == 0.0 -> if (y == 0) 0.0 else Double.NEGATIVE_INFINITY
            failures.isInfinite() -> y * ln(y.toFloat()) - mean - MoreMath.factorialLog(y)
            else -> Gamma.logGamma(y + failures) - MoreMath.factorialLog(y) + fLogf - failures * logMeanPlusFailure +
                    y * (ln(mean) - logMeanPlusFailure)
        }
    }

    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        val X = generateDesignMatrix(df)
        val yInt = df.sliceAsInt(df.labels[d])
        val y = DoubleArray(yInt.size) { yInt[it].toDouble() }.asF64Array()
        val iterMax = 100
        val tol = 1e-8
        var beta0 = regressionCoefficients
        var beta1 = regressionCoefficients
        for (i in 0 until iterMax) {
            val eta = WLSRegression.calculateEta(X, beta0)
            val prevFailures = failures
            failures = NegativeBinomialDistribution.fitNumberOfFailures(y, weights, meanInPlace(eta.copy()), failures)
            fLogf = failures * ln(failures)
            val (z, W) = zW(y, eta)
            W *= weights
            beta1 = WLSRegression.calculateBeta(X, z, W)
            if ((beta1.zip(beta0) { a, b -> abs(a - b) }).sum() + abs(failures - prevFailures) < tol) {
                break
            }
            beta0 = beta1
        }
        regressionCoefficients = beta1
    }
}