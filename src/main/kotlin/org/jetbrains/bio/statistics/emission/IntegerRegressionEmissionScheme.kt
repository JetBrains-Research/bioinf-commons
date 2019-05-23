package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import java.util.*

/**
 *
 *Regression with integer support.
 * Examples: Poisson, binomial.
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
abstract class IntegerRegressionEmissionScheme(val covariateLabels: List<String>, regressionCoefficients: DoubleArray)
    : EmissionScheme {
    var regressionCoefficients: DoubleArray = regressionCoefficients
        protected set
    abstract val link: (Double) -> Double
    abstract val linkDerivative: (Double) -> Double
    abstract val linkVariance: (Double) -> Double

    /**
     * IRLS algorithm is used for coefficients prediction.
     * For more details: https://bwlewis.github.io/GLM/
     *
     * @param d - number of column which contains observations
     */
    override fun update(df: DataFrame, d:Int, weights: F64Array){
        val y = Arrays.stream(df.sliceAsInt(df.labels[d].intern())).asDoubleStream().toArray()
        var xFromCovariates = emptyArray<DoubleArray>()
        for (label in covariateLabels){
            xFromCovariates += df.sliceAsDouble(label)
        }

        val x  = Array2DRowRealMatrix(xFromCovariates).transpose().data

        val wlr = WLSMultipleLinearRegression()
        val weights_0 = DoubleArray(x.size)
        Arrays.fill(weights_0, 0.0)
        wlr.newSampleData(y, x, weights_0)

        val iterMax = 5
        val tol = 1e-8


        var X0:RealVector = ArrayRealVector(regressionCoefficients)
        var X1:RealVector = ArrayRealVector(regressionCoefficients)

        //
        val X = wlr.getx()
        val Y = ArrayRealVector(y)

        for (i in 0 until iterMax) {
            val eta = X.operate(X0)
            val countedLink = eta.map { link(it) }
            val countedLinkDeriv = eta.map { linkDerivative(it) }
            val z:RealVector = eta.add(Y.subtract(countedLink).ebeDivide(countedLinkDeriv))
            val countedLinkVar = countedLink.map { linkVariance(it) }
            val W = countedLinkDeriv
                    .ebeMultiply(countedLinkDeriv)
                    .ebeDivide(countedLinkVar)
                    .ebeMultiply(ArrayRealVector(weights.toDoubleArray()))
                    .toArray()

            wlr.newSampleData(z.toArray(), x, W)

            X1 = ArrayRealVector(wlr.calculateBeta())
            if (X1.subtract(X0).l1Norm < tol) {
                break
            }
            X0 = X1
        }
        this.regressionCoefficients = X1.toArray()
        println(regressionCoefficients.toList())
    }
}