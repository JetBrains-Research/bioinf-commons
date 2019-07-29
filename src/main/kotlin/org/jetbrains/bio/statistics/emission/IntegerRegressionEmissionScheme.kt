package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.*
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.sequence.TwoBitSequence
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path
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
            val countedLink = eta.map { link(it) }
            val countedLinkDerivative = eta.map { linkDerivative(it) }
            val z = DoubleArray(eta.dimension)
            {eta.getEntry(it) + (y[it] - countedLink.getEntry(it)) / countedLinkDerivative.getEntry(it)}
            val countedLinkVar = countedLink.map { linkVariance(it) }
            W = DoubleArray(countedLink.dimension)
            {countedLink.getEntry(it) * countedLink.getEntry(it) / countedLinkVar.getEntry(it) * weights[it]}
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

fun getIntCover(chr1: Chromosome, coverage: Coverage, bin: Int): IntArray {
    val len = (chr1.length - 1) / bin + 1
    val cover = IntArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * bin, (i + 1) * bin, chr1))
    }
    cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * bin, chr1.length, chr1))
    return cover
}
fun getDoubleCover(chr1: Chromosome, coverage: Coverage, bin: Int): DoubleArray {
    val len = (chr1.length - 1) / bin + 1
    val cover = DoubleArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage
                .getBothStrandsCoverage(ChromosomeRange(i * bin, (i + 1) * bin, chr1))
                .toDouble()
    }
    cover[len - 1] = coverage
            .getBothStrandsCoverage(ChromosomeRange((len-1) * bin, chr1.length, chr1))
            .toDouble()
    return cover
}
fun getGC(chr1: Chromosome, bin: Int): DoubleArray {
    val len = (chr1.length - 1) / bin + 1
    val seq: TwoBitSequence = chr1.sequence
    val GCcontent = DoubleArray(len)
    for (i in 0 until len - 1) {
        GCcontent[i] = seq.substring(i*bin, (i + 1)*bin).count { it == 'c' || it == 'g' }.toDouble()/bin
    }
    GCcontent[len - 1] = seq
            .substring((len-1)*bin, seq.length)
            .count { it == 'c'|| it == 'g' }
            .toDouble()/( seq.length - (len-1)*bin)
    return GCcontent
}

fun getMappability(chr1: Chromosome, path_mappability: Path, bin: Int): DoubleArray {
    val mapSummary = BigWigFile
            .read(path_mappability)
            .summarize(chr1.name, 0, chr1.length - chr1.length%bin, numBins = (chr1.length - 1)/bin)
    val result = DoubleArray (mapSummary.size + 1) {
        if (it < mapSummary.size) mapSummary[it].sum/bin else 1.0}
    result[mapSummary.size] = BigWigFile
            .read(path_mappability)
            .summarize(chr1.name, chr1.length - chr1.length%bin, 0)[0].sum / chr1.length%bin
    return result}
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


