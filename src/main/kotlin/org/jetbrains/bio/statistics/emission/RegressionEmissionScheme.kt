package org.jetbrains.bio.statistics.emission

import htsjdk.samtools.util.SequenceUtil.n
import org.apache.commons.math3.distribution.PoissonDistribution
import org.apache.commons.math3.linear.*
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.DataFrameBuilder
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.mixture.MLFreeMixture
import java.util.*
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.jetbrains.bio.viktor.toF64Array
import org.tukaani.xz.check.None
import java.util.function.IntPredicate
import java.util.function.IntSupplier
import javax.lang.model.type.ArrayType
import kotlin.math.exp
import kotlin.random.Random

abstract class RegressionEmissionScheme(covariateLabels: List<String>, regressionCoefficients: DoubleArray) : EmissionScheme {
    val covariateLabels: List<String> = covariateLabels
    var regressionCoefficients: DoubleArray = regressionCoefficients
        protected set
    abstract val link: (Double) -> Double
    abstract val linkDerivative: (Double) -> Double
    abstract val linkVariance: (Double) -> Double

    override fun update(df: DataFrame, d:Int, weights: F64Array){
        var y = df.sliceAsDouble(df.labels[d].intern())
        var xFromCovariates = emptyArray<DoubleArray>()
        for (label in covariateLabels){
            xFromCovariates += df.sliceAsDouble(label)
        }

        var x  = Array2DRowRealMatrix(xFromCovariates).transpose().data

        var wlr = WLSMultipleLinearRegression()
        val weights_0 = DoubleArray(x.size)
        Arrays.fill(weights_0, 0.0)
        wlr.newSampleData(y, x, weights_0)

        val iterMax = 250
        val tol = 1e-8

        //fill x0 with zeros
        val x0 = DoubleArray(x[0].size+1)
        Arrays.fill(x0, 0.0)

        var X0:RealVector = ArrayRealVector(x0)
        var X1:RealVector = ArrayRealVector(x0)

        //
        val X = wlr.getx()
        val Y = ArrayRealVector(y)

        for (i in 0 until iterMax) {
            val eta = X.operate(X0)
            val countedLink = eta.map { link(it) }
            val countedLinkDeriv = eta.map { linkDerivative(it) }
            val z:RealVector = eta.add(Y.subtract(countedLink).ebeDivide(countedLinkDeriv))
            val countedLinkVar = countedLink.map { linkVariance(it) }
            val W = countedLinkDeriv.ebeMultiply(countedLinkDeriv).ebeDivide(countedLinkVar).ebeMultiply(ArrayRealVector(weights.data)).toArray()
            
            wlr.newSampleData(z.toArray(), x, W)

            X1 = ArrayRealVector(wlr.calculateBeta())
            if (X1.subtract(X0).l1Norm < tol) {
                break
            }
            X0 = X1
        }
        this.regressionCoefficients = X1.toArray()
    }
}

class PoissonRegressionEmissionScheme (
        covariateLabels: List<String>,
        regressionCoefficients: DoubleArray,
        override val degreesOfFreedom: Int): RegressionEmissionScheme(covariateLabels, regressionCoefficients) {

    override val link: (Double) -> Double = {Math.exp(it)}
    override val linkDerivative: (Double) -> Double = {Math.exp(it)}
    override val linkVariance: (Double) -> Double = {it}
    var sample =  emptyArray<Double>()
        private set


    override fun sample(df: DataFrame, d: Int, fill: IntPredicate) {
        val observations = df.sliceAsDouble(df.labels[d])
        var observation = DoubleArray(covariateLabels.size)
        (0 until df.rowsNumber).forEach { row ->
            if(fill.test(row)){
                this.covariateLabels.forEachIndexed { index, label ->
                    observation[index] = df.getAsDouble(row, label)
                }
            }
            observations[row] = Sampling.samplePoisson(exp(regressionCoefficients.zip(observation) { a, b -> a*b }.sum())).toDouble()
        }
    }

    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        TODO("not implemented")
    }
}


// 0 - zero-inflated component
// 1 - LOW
// 2 - HIGH
class ZeroPoissonMixture (weights: F64Array, covariateLabels: List<String>, regressionCoefficients: Array<DoubleArray>): MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights){
    val covariateLabels: List<String> = covariateLabels
    var regressionCoefficients: Array<DoubleArray> = regressionCoefficients
        protected set

    val components: MutableList<EmissionScheme> = mutableListOf<EmissionScheme>(
            ConstantIntegerEmissionScheme(0),
            PoissonRegressionEmissionScheme(
                    covariateLabels = covariateLabels,
                    regressionCoefficients = regressionCoefficients[0],
                    degreesOfFreedom = covariateLabels.size),
            PoissonRegressionEmissionScheme(
                    covariateLabels = covariateLabels,
                    regressionCoefficients = regressionCoefficients[1],
                    degreesOfFreedom = covariateLabels.size))


    override fun getEmissionScheme(i: Int, d: Int): EmissionScheme {
        return components[i]
    }
}


fun main(args: Array<String>) {
    // 1
    var covar = DataFrame()
            .with("x1", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })
            .with("x2", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })
            .with("y", DoubleArray(1000000))
    //проверка, что с внешними весами все еще работает
    var regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(4.0, -2.0), 2)
    val pred = IntPredicate {true}


    regrES.sample(covar, 2, pred)
    regrES.update(covar, 2, DoubleArray(1000000, {1.0}).asF64Array())

    print("Beta: ${regrES.regressionCoefficients.asList()}")

    // MLFreeMixture
}
