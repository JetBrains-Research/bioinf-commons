package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.RealVector
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.util.function.IntPredicate
import kotlin.math.abs

/**
 *
 *Regression with integer support.
 * Examples: Poisson, binomial.
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
    abstract fun link(x: Double): Double
    abstract fun linkDerivative(x: Double): Double
    abstract fun linkVariance(x: Double): Double
    abstract fun sampler(x: Double): Int
    override val degreesOfFreedom: Int = regressionCoefficients.size

    open fun linkInPlace(x: F64Array) {
        for(i in 0 until x.size) {
            x[i] = link(x[i])
        }
    }

    open fun linkDerivativeInPlace(x: F64Array) {
        for(i in 0 until x.size) {
            x[i] = linkDerivative(x[i])
        }
    }

    open fun linkVarianceInPlace(x: F64Array) {
        for(i in 0 until x.size) {
            x[i] = linkVariance(x[i])
        }
    }

    open fun zW(y: F64Array, eta: F64Array): Pair<F64Array, F64Array> {
        val countedLink = eta.copy().apply { linkInPlace(this) }
        val countedLinkDerivative = eta.copy().apply { linkDerivativeInPlace(this) }
        val z =
                eta
                        .apply { plusAssign(
                                y.copy()
                                        .apply { minusAssign(countedLink)}
                                        .apply { divAssign(countedLinkDerivative)})}
        val countedLinkVar = countedLink.copy().apply { linkVarianceInPlace(this) }
        val W = countedLink.apply {timesAssign(countedLink)}.apply {divAssign(countedLinkVar)}
        return z to W
    }
    /**
     * IRLS algorithm is used for coefficients prediction.
     * For more details: https://bwlewis.github.io/GLM/
     *
     * @param d - number of column which contains observations
     */
    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        val X = WLSRegression.designMatrix(Array(covariateLabels.size) { df.sliceAsDouble(covariateLabels[it]) })
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
     * @param t - number of row
     */
    fun getPredictor(df: DataFrame, t: Int): Double {
        var res = regressionCoefficients[0]
        covariateLabels.forEachIndexed { index, label ->
            res += df.getAsDouble(t, label) * regressionCoefficients[index + 1]
        }
        return res
    }

    /**
     * @param d - number of column in which we want to sample
     * @param fill - predicate which marks rows we need to use for sampling
     */
    override fun sample(df: DataFrame, d: Int, fill: IntPredicate) {
        val observations = df.sliceAsInt(df.labels[d])
        (0 until df.rowsNumber).forEach { row ->
            if (fill.test(row)) {
                observations[row] = sampler(link(getPredictor(df, row)))
            }
        }
    }

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
        val inverse = LUDecomposition(R.transpose().multiply(XTWXI).multiply(R)).solver.inverse
        val Fstat = inverse.operate(RBetaMinusr).dotProduct(RBetaMinusr) / (r.dimension*sigma2)
        val pVal = 1 - FDistribution(r.dimension.toDouble(), (X[0].size - X.size).toDouble()).cumulativeProbability(Fstat)
        return pVal
    }
}
