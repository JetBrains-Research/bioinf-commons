package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.linear.*
import org.apache.commons.math3.util.FastMath
import org.apache.log4j.Level
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.MoreMath
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.mixture.MLFreeMixture
import org.jetbrains.bio.util.Logs
import java.util.*
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.util.function.IntPredicate
import kotlin.math.*
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
            val W = countedLinkDeriv.ebeMultiply(countedLinkDeriv).ebeDivide(countedLinkVar).ebeDivide(ArrayRealVector(weights.toDoubleArray())).toArray()
            
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

        val observations = df.sliceAsInt(df.labels[d])
        var observation = DoubleArray(covariateLabels.size + 1, {1.0})
        (0 until df.rowsNumber).forEach { row ->
            if(fill.test(row)){
                this.covariateLabels.forEachIndexed { index, label ->
                    observation[index + 1] = df.getAsDouble(row, label)
                }
                observations[row] = Sampling.samplePoisson(exp(regressionCoefficients.zip(observation) { a, b -> a*b }.sum())).toInt()
            }
        }
    }

    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        val logLambda = regressionCoefficients
                .zip(
                        doubleArrayOf(1.0)
                                .plus(
                                        df.rowAsDouble(t)
                                                .filterIndexed
                                                { index, d -> covariateLabels.contains(df.labels[index]) }))
                { a, b -> a * b }
                .sum()
        return(df.getAsInt(t, df.labels[d])* logLambda - MoreMath.factorialLog(df.getAsInt(t, df.labels[d]))) - FastMath.exp(logLambda)

    }

    override fun update(df: DataFrame, d:Int, weights: F64Array){
        var y = Arrays.stream(df.sliceAsInt(df.labels[d].intern())).asDoubleStream().toArray()
        var xFromCovariates = emptyArray<DoubleArray>()
        for (label in covariateLabels){
            xFromCovariates += df.sliceAsDouble(label)
        }

        var x  = Array2DRowRealMatrix(xFromCovariates).transpose().data

        var wlr = WLSMultipleLinearRegression()
        val weights_0 = DoubleArray(x.size)
        Arrays.fill(weights_0, 0.0)
        wlr.newSampleData(y, x, weights_0)

        val iterMax = 5
        val tol = 1e-8

        //fill x0 with zeros
        //val x0 = DoubleArray(x[0].size+1)
        //Arrays.fill(x0, 0.0)

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
            val W = countedLinkDeriv.ebeMultiply(countedLinkDeriv).ebeDivide(countedLinkVar).ebeMultiply(ArrayRealVector(weights.toDoubleArray())).toArray()

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



// 0 - zero-inflated component
// 1 - LOW
// 2 - HIGH
class ZeroPoissonMixture (weights: F64Array, covariateLabels: List<String>, regressionCoefficients: Array<DoubleArray>): MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights){
    val covariateLabels: List<String> = covariateLabels
    var regressionCoefficients: Array<DoubleArray> = regressionCoefficients
        protected set

    val components: List<EmissionScheme> = listOf<EmissionScheme>(
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

    fun sample(df:DataFrame, d: IntArray) {
        val states = sampleStates(df.rowsNumber)
        for (t in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, t).sample(df, d[t], IntPredicate { states[it] == i })
            }
        }
    }

    //For this moment suppose that numDimension = 1
    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        super.updateParameters(df, gammas)
        for (d in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, d).update(df, d, gammas.V[i])
            }
        }
    }
}


fun main(args: Array<String>) {
    //1

    var covar = DataFrame()
            .with("y", IntArray(1000000))
            .with("x1", DoubleArray(1000000) { Random.nextDouble(0.5, 1.0) })
            /*.with("x2", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })*/
    /*
    //проверка, что с внешними весами все еще работает
    var regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(4.0, -2.0), 2)
    val pred = IntPredicate {true}


    regrES.sample(covar, 0, pred)
    println("Update")
    regrES.update(covar, 0, DoubleArray(1000000, {1.0}).asF64Array())

    print("Beta: ${regrES.regressionCoefficients.asList()}")
*/
    // MLFreeMixture
    Logs.addConsoleAppender(Level.DEBUG)
    var mix = ZeroPoissonMixture(
            doubleArrayOf(0.4, 0.5, 0.1).asF64Array(),
            listOf("x1"/*, "x2"*/),
            arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(2.0, 2.0)))

    mix.sample(covar, intArrayOf(0))

    var yaMix = ZeroPoissonMixture(
            doubleArrayOf(1/3.0, 1/3.0, 1/3.0).asF64Array(),
            listOf("x1"/*, "x2"*/),
            arrayOf(doubleArrayOf(0.0, 1.5), doubleArrayOf(1.6, 0.0)))
    yaMix.fit(Preprocessed.of(covar), 10e-3, 12)
}
