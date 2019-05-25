package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate

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
    abstract val link: (Double) -> Double
    abstract val linkDerivative: (Double) -> Double
    abstract val linkVariance: (Double) -> Double
    abstract val sampler: (Double) -> Int
    override val degreesOfFreedom: Int = regressionCoefficients.size

    /**
     * IRLS algorithm is used for coefficients prediction.
     * For more details: https://bwlewis.github.io/GLM/
     *
     * @param d - number of column which contains observations
     */
    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        val x = Array2DRowRealMatrix(covariateLabels.map { df.sliceAsDouble(it) }.toTypedArray())
                .transpose()
                .data

        // needed here to add intercept to X
        val wlr = WLSMultipleLinearRegression()
        wlr.newSampleData(DoubleArray(x.size), x, DoubleArray(x.size))
        val X = wlr.getx()
        val y = ArrayRealVector(df.sliceAsInt(df.labels[d]).map { it.toDouble() }.toDoubleArray())

        val iterMax = 5
        val tol = 1e-8


        var beta0: RealVector = ArrayRealVector(regressionCoefficients)
        var beta1: RealVector = ArrayRealVector(regressionCoefficients)

        for (i in 0 until iterMax) {
            val eta = X.operate(beta0)
            val countedLink = eta.map { link(it) }
            val countedLinkDerivative = eta.map { linkDerivative(it) }
            val z: RealVector = eta.add(y.subtract(countedLink).ebeDivide(countedLinkDerivative))
            val countedLinkVar = countedLink.map { linkVariance(it) }
            val W = countedLinkDerivative
                    .ebeMultiply(countedLinkDerivative)
                    .ebeDivide(countedLinkVar)
                    .ebeMultiply(ArrayRealVector(weights.toDoubleArray()))
                    .toArray()

            wlr.newSampleData(z.toArray(), x, W)

            beta1 = ArrayRealVector(wlr.calculateBeta())
            if (beta1.getL1Distance(beta0) < tol) {
                break
            }
            beta0 = beta1
        }
        regressionCoefficients = beta1.toArray()
    }

    /**
     * @param t - number of row
     */
    fun getPredictor(df: DataFrame, t: Int): Double {
        val observation = DoubleArray(covariateLabels.size + 1) { 1.0 }
        covariateLabels.forEachIndexed { index, label ->
            observation[index + 1] = df.getAsDouble(t, label)
        }
        return regressionCoefficients.zip(observation) { a, b -> a * b }.sum()
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
}