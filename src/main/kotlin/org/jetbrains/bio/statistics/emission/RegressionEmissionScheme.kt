package org.jetbrains.bio.statistics.emission

import htsjdk.samtools.util.SequenceUtil.n
import org.apache.commons.math3.distribution.PoissonDistribution
import org.apache.commons.math3.linear.*
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.DataFrameBuilder
import org.jetbrains.bio.statistics.distribution.Sampling
import java.util.*
import org.jetbrains.bio.viktor.F64Array
import org.tukaani.xz.check.None
import java.util.function.IntPredicate
import kotlin.math.exp
import kotlin.random.Random

abstract class RegressionEmissionScheme(covariateLabels: List<String>, regressionCoefficients: DoubleArray) : EmissionScheme {
    val covariateLabels: List<String> = covariateLabels
    var regressionCoefficients: DoubleArray = regressionCoefficients
        protected set
    abstract val link: (Double) -> Double
    abstract val linkDerivative: (Double) -> Double
    abstract val linkVariance: (Double) -> Double

    fun update(df: DataFrame, t:Int){
        var y = df.sliceAsDouble(df.labels[t].intern())
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

        var X0: RealVector = ArrayRealVector(x0)
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
            val W = countedLinkDeriv.ebeMultiply(countedLinkDeriv).ebeDivide(countedLinkVar).toArray()

            wlr.newSampleData(z.toArray(), x, W)

            X1 = wlr.calculateBeta()
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
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun update(df: DataFrame, d: Int, weights: F64Array) {
    }





}


fun main(args: Array<String>) {
    // 1
    var covar = DataFrame()
            .with("x1", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })
            .with("x2", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })
            .with("y", DoubleArray(1000000))
    //2
    var regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(4.0, -2.0), 2)
    val pred = IntPredicate {true}

    val startTime = System.nanoTime()
    var maxDiff = doubleArrayOf(0.0, 0.0, 0.0)
    val trueValues = doubleArrayOf(0.0, 4.0, -2.0)
    for(j in 1..100){

        regrES.sample(covar, 2, pred)

        regrES.update(covar, 2)

        regrES.regressionCoefficients.forEachIndexed { index, d ->
            if (maxDiff[index] <= Math.abs(d - trueValues[index])) maxDiff[index] = Math.abs(d - trueValues[index])
        }


        regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(4.0, -2.0), 2)

    }
    print(maxDiff.toList())
    println("Mean working time: ${(System.nanoTime() - startTime)/100}")
}
