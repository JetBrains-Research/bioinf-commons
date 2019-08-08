package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.*
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
    @Transient var W = DoubleArray(0)
    abstract fun link(x: Double): Double
    abstract fun linkDerivative(x: Double): Double
    abstract fun linkVariance(x: Double): Double
    abstract fun sampler(x: Double): Int
    override val degreesOfFreedom: Int = regressionCoefficients.size

    open fun linkInPlace(x: DoubleArray) {
        for(i in 0 until x.size) {
            x[i] = link(x[i])
        }
    }

    open fun linkDerivativeInPlace(x:DoubleArray) {
        for(i in 0 until x.size) {
            x[i] = linkDerivative(x[i])
        }
    }

    open fun linkVarianceInPlace(x: DoubleArray) {
        for(i in 0 until x.size) {
            x[i] = linkVariance(x[i])
        }
    }
    /**
     * IRLS algorithm is used for coefficients prediction.
     * For more details: https://bwlewis.github.io/GLM/
     *
     * @param d - number of column which contains observations
     */
    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        val x = Array(covariateLabels.size) {df.sliceAsDouble(covariateLabels[it])}
        val wlr = WLSMultipleLinearRegression()
        wlr.newSampleData(DoubleArray(x[0].size), x, DoubleArray(x[0].size))
        val X = wlr.x
        val yInt = df.sliceAsInt(df.labels[d])
        val y = DoubleArray (yInt.size) {yInt[it].toDouble()}
        val iterMax = 5
        val tol = 1e-8
        var beta0: RealVector = ArrayRealVector(regressionCoefficients, false)
        var beta1: RealVector = ArrayRealVector(regressionCoefficients, false)
        for (i in 0 until iterMax) {
            val eta = X.operate(beta0)
            val countedLink = eta.toArray().apply { linkInPlace(this) }
            val countedLinkDerivative = eta.toArray().apply { linkDerivativeInPlace(this) }
            val z = DoubleArray(eta.dimension)
            {eta.getEntry(it) + (y[it] - countedLink[it]) / countedLinkDerivative[it]}
            val countedLinkVar = countedLink.apply { linkVarianceInPlace(this) }
            W = DoubleArray(countedLink.size)
            {countedLink[it] * countedLink[it] / countedLinkVar[it] * weights[it]}
            wlr.newSampleData(z, W)
            beta1 = wlr.calculateBeta()
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
        var res = regressionCoefficients[0]
        covariateLabels.forEachIndexed { index, label ->
            res += df.getAsDouble(t, label)*regressionCoefficients[index + 1]
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
        val wlr = WLSMultipleLinearRegression()
        wlr.newSampleData(y, x, W)
        val X = wlr.x
        val residuals = X
                .operate(ArrayRealVector(regressionCoefficients, false))
                .subtract(ArrayRealVector(y, false))
        val sigma = residuals
                .dotProduct(DiagonalMatrix(W)
                .operate(residuals)) / (X.rowDimension - X.columnDimension)
        val inverseXTX = wlr.calculateBetaVariance()
        val RBeta = R.transpose().operate(ArrayRealVector(regressionCoefficients))
        val RBetaMinusr = RBeta.subtract(r)
        val inverse = LUDecomposition(R.transpose().multiply(inverseXTX).multiply(R)).solver.inverse
        val Fstat = inverse.operate(RBetaMinusr).dotProduct(RBetaMinusr) / (r.dimension*sigma)
        val pVal = 1 - FDistribution(r.dimension.toDouble(), (X.rowDimension - X.columnDimension).toDouble()).cumulativeProbability(Fstat)

        return pVal
    }
}

/*
fun getLocalBGEstimate(chr1: Chromosome, path_mappability: Path, coverage: Coverage): DoubleArray {
    val shiftLength = 2500
    val slidingWindowSize = 100000
    val numOfSlidingWindows = (chr1.length - slidingWindowSize)/shiftLength + 1
    val mapSource = BigWigFile.read(path_mappability)
    (0 until numOfSlidingWindows).forEach {
        coverage.getBothStrandsCoverage(
                ChromosomeRange(
                        it*shiftLength,
                        it*shiftLength + slidingWindowSize,
                        chr1))* slidingWindowSize / mapSource.summarize(
                        chr1.name,
                        it*shiftLength,
                        it*shiftLength + slidingWindowSize)[0].sum
    }
    val len = (chr1.length - 1) / 200 + 1
    val result = DoubleArray (len).forEach { it =  }
} */


