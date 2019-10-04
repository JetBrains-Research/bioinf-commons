package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.*
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.util.function.IntPredicate
import kotlin.math.abs

/**
 *
 * Regression with integer support.
 * Examples: Poisson, binomial.
 *
 * @param covariateLabels the labels that will be used to procure covariates from the supplied [DataFrame]
 * @param regressionCoefficients the regression coefficients. Note that
 *     regressionCoefficients.size == covariateLabels.size + 1
 * because the 0th coefficient corresponds to the intercept, the 1st coefficient corresponds to 0th label etc.
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
abstract class IntegerRegressionEmissionScheme(
        covariateLabels: List<String>, regressionCoefficients: DoubleArray
) : EmissionScheme {
    val covariateLabels = covariateLabels.map { it.intern() }
    var regressionCoefficients: DoubleArray = regressionCoefficients
        protected set
    @Transient var W = DoubleArray(0)

    /**
     * Provides the mean function (also known as the inverted link function), h = g^{-1}.
     *
     * The generalized linear model can be defined through the link function.
     *
     * η = g(μ),
     * where μ is the distribution mean, η = Xβ is the linear predictor, and g is the link function.
     *
     * Conversely,
     * μ = h(η),
     * where h = g^{-1} is called the inverse link function, or sometimes the mean function.
     *
     * @param eta η = Xβ, the linear predictor
     * @return μ = h(η), the distribution mean
     *
     * Note that it's often more efficient to provide an implementation for the vector-based [meanInPlace].
     * It's also possible to override [zW] and skip other link-related implementations altogether.
     */
    abstract fun mean(eta: Double): Double

    /**
     * The mean derivative function, h'(η).
     *
     * η = Xβ is the linear predictor. See [mean] for more information.
     *
     * @param eta η = Xβ, the linear predictor
     * @return h'(η)
     *
     * Note that it's often more efficient to provide an implementation for the vector-based [meanDerivativeInPlace].
     * It's also possible to override [zW] and skip other link-related implementations altogether.
     */
    abstract fun meanDerivative(eta: Double): Double

    /**
     * The mean variance function, var h(η)
     *
     * h(η) is the mean function, η = Xβ is the linear predictor. See [mean] for more information.
     *
     * @param mean μ = h(η)
     * @return var(h(η))
     *
     * Note that it's often more efficient to provide an implementation for the vector-based [meanVarianceInPlace].
     * It's also possible to override [zW] and skip other mean-related implementations altogether.
     */
    abstract fun meanVariance(mean: Double): Double

    /**
     * Samples from the regression given the distribution mean.
     *
     * @param mean the distribution mean μ = h(η)
     */
    abstract fun sampler(mean: Double): Int

    override val degreesOfFreedom: Int = regressionCoefficients.size

    /**
     * Apply [mean] function h(η) in place.
     *
     * Can be overridden to use more efficient vector operations.
     * It's also possible to override [zW] and skip other mean-related implementations altogether.
     *
     * @param eta η = Xβ, the linear predictor vector. Will be modified to contain μ = h(η).
     * @return [eta] for easier call chaining.
     */
    open fun meanInPlace(eta: F64Array): F64Array {
        for(i in 0 until eta.size) {
            eta[i] = mean(eta[i])
        }
        return eta
    }

    /**
     * Apply [meanDerivative] function h'(η) in place.
     *
     * Can be overridden to use more efficient vector operations.
     * It's also possible to override [zW] and skip other mean-related implementations altogether.
     *
     * @param eta η = Xβ, the linear predictor vector. Will be modified to contain h'(η).
     * @return [eta] for easier call chaining.
     */
    open fun meanDerivativeInPlace(eta: F64Array): F64Array {
        for(i in 0 until eta.size) {
            eta[i] = meanDerivative(eta[i])
        }
        return eta
    }

    /**
     * Apply [meanVariance] function var(h(η)) in place.
     *
     * Can be overridden to use more efficient vector operations.
     * It's also possible to override [zW] and skip other mean-related implementations altogether.
     *
     * @param mean μ = h(η), the distribution mean vector. Will be modified to contain var(h(η)).
     * @return [mean] for easier call chaining.
     */
    open fun meanVarianceInPlace(mean: F64Array): F64Array {
        for(i in 0 until mean.size) {
            mean[i] = meanVariance(mean[i])
        }
        return mean
    }

    /**
     * Calculates z and W for IRLS as defined in https://bwlewis.github.io/GLM/ .
     *
     * Note that the document linked above uses definitions and denotations that differ from
     * the widely accepted ones.
     * In particular, "GLMs, abridged" uses the term "link function" and the symbol "g" to describe
     * what is normally called the mean function (or the inverse link function) and denoted by "h".
     * It also uses "b" instead of "y" for the response vector and "A" instead of "X" for the model matrix.
     *
     * z = η + (y - h(η)) / h'(η)
     * W = diag(h'(η)^2 / var(h(η)))
     *
     * This is the only function necessary for IRLS and [update] to function. You can override it directly
     * instead of providing implementations for [mean] and friends if there's a more efficient way to compute
     * z and W. For example, with a logarithmic link, h = h' = var h, which allows to simplify W.
     *
     * @param y the response vector
     * @param eta η = Xβ, the linear predictor. WARNING: [eta] can be modified to reduce copying!
     * @return z and W as a [Pair] (the diagonal matrix W is stored as a vector)
     */
    open fun zW(y: F64Array, eta: F64Array): Pair<F64Array, F64Array> {
        val countedLink = meanInPlace(eta.copy())
        val countedLinkDerivative = meanDerivativeInPlace(eta.copy())

        eta += (y - countedLink).apply { divAssign(countedLinkDerivative)}

        meanVarianceInPlace(countedLink)

        countedLinkDerivative *= countedLinkDerivative
        countedLinkDerivative /= countedLink
        return eta to countedLinkDerivative
    }

    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        val X = generateDesignMatrix(df)
        val yInt = df.sliceAsInt(df.labels[d])
        val y = DoubleArray (yInt.size) {yInt[it].toDouble()}.asF64Array()
        val iterMax = 5
        val tol = 1e-8
        var beta0 = regressionCoefficients
        var beta1 = regressionCoefficients
        for (i in 0 until iterMax) {
            val eta = WLSRegression.calculateEta(X, beta0)
            val (z, W) = zW(y, eta)
            W *= weights
            beta1 = WLSRegression.calculateBeta(X, z, W)
            if ((beta1.zip(beta0) { a, b -> abs(a - b) }).sum() < tol) {
                break
            }
            beta0 = beta1
        }
        regressionCoefficients = beta1
    }

    /**
     * Generates the model matrix from the dataframe.
     *
     * Extracts the labelled covariates from the supplied dataframe and prepends the intercept column.
     */
    protected fun generateDesignMatrix(df: DataFrame) =
            WLSRegression.designMatrix(Array(covariateLabels.size) { df.sliceAsDouble(covariateLabels[it]) })

    /**
     * Calculates η = Xβ, the linear predictor.
     *
     * Uses a row from the provided dataframe to get the covariate values.
     *
     * @param df dataframe
     * @param t row number
     */
    fun getPredictor(df: DataFrame, t: Int): Double {
        var res = regressionCoefficients[0]
        covariateLabels.forEachIndexed { index, label ->
            res += df.getAsDouble(t, label) * regressionCoefficients[index + 1]
        }
        return res
    }

    override fun sample(df: DataFrame, d: Int, fill: IntPredicate) {
        val observations = df.sliceAsInt(df.labels[d])
        (0 until df.rowsNumber).forEach { row ->
            if (fill.test(row)) {
                observations[row] = sampler(mean(getPredictor(df, row)))
            }
        }
    }

    /**
     * Calculates the f-test p-value.
     *
     * Performs the f-test for the null hypothesis that Rβ = r and returns the p-value.
     */
    fun Ftest(df: DataFrame, d: Int, R: RealMatrix, r: RealVector): Double {
        val x = Array(covariateLabels.size) {df.sliceAsDouble(covariateLabels[it])}
        val yInt = df.sliceAsInt(df.labels[d])
        val y = DoubleArray (yInt.size) {yInt[it].toDouble()}
        val X = WLSRegression.designMatrix(x)
        val residuals = WLSRegression.calculateEta(X, regressionCoefficients).apply { minusAssign(y.asF64Array()) }
        val sigma2 = residuals.dot(W.asF64Array() * residuals) / (X[0].size - X.size)
        val XTWXI = WLSRegression.calculateBetaVariance(X, W.asF64Array())
        val RBeta = R.transpose().operate(ArrayRealVector(regressionCoefficients))
        val RBetaMinusr = RBeta.subtract(r)
        val inverse = LUDecomposition(R.transpose().multiply(Array2DRowRealMatrix(XTWXI)).multiply(R)).solver.inverse
        val Fstat = inverse.operate(RBetaMinusr).dotProduct(RBetaMinusr) / (r.dimension * sigma2)
        val cumulativeProb = FDistribution(
            r.dimension.toDouble(), (X[0].size - X.size).toDouble()
        ).cumulativeProbability(Fstat)
        return 1 - cumulativeProb
    }
}
