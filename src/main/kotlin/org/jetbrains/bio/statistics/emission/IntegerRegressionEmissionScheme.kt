package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.*
import org.apache.log4j.Level
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.sequence.TwoBitSequence
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.mixture.ZeroPoissonMixture
import org.jetbrains.bio.util.Logs
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.lang.System.arraycopy
import java.nio.file.Paths
import java.util.function.IntPredicate
import org.jetbrains.bio.big.BigWigFile
import java.nio.file.Path
import kotlin.math.ln

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
    var omega = emptyArray<Double>()
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
        val X = wlr.getX()
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
            val W = DoubleArray(countedLink.dimension)
            {countedLink.getEntry(it) * countedLink.getEntry(it) / countedLinkVar.getEntry(it) * weights[it]}
            wlr.newSampleData(z, W)
            beta1 = wlr.calculateBeta()
            omega = W.toTypedArray()
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
        val lnY = DoubleArray (y.size) {ln(y[it])}
        wlr.newSampleData(y, x, omega.toDoubleArray())
        val X = wlr.getX()
        val residuals = X
                .operate(ArrayRealVector(regressionCoefficients, false))
                .subtract(ArrayRealVector(y, false))
        val sigma = residuals
                .dotProduct(DiagonalMatrix(omega.toDoubleArray())
                .operate(residuals)) / (X.getRowDimension() - X.getColumnDimension())
        val inverseXTX = wlr.calculateBetaVariance()
        val RBeta = R.transpose().operate(ArrayRealVector(regressionCoefficients))
        val RBetaMinusr = RBeta.subtract(r)
        val inverse = LUDecomposition(R.transpose().multiply(inverseXTX).multiply(R)).solver.inverse
        val Fstat = inverse.operate(RBetaMinusr).dotProduct(RBetaMinusr) / (r.dimension*sigma)
        val pVal = 1 - FDistribution(r.dimension.toDouble(), (X.rowDimension - X.columnDimension).toDouble()).cumulativeProbability(Fstat)

        return pVal
    }
}

fun getIntCover(chr1: Chromosome, coverage: Coverage): IntArray {
    val len = (chr1.length - 1) / 200 + 1
    val cover = IntArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * 200, (i + 1) * 200, chr1))
    }
    cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * 200, chr1.length, chr1))
    return cover
}
fun getDoubleCover(chr1: Chromosome, coverage: Coverage): DoubleArray {
    val len = (chr1.length - 1) / 200 + 1
    val cover = DoubleArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage
                .getBothStrandsCoverage(ChromosomeRange(i * 200, (i + 1) * 200, chr1))
                .toDouble()
    }
    cover[len - 1] = coverage
            .getBothStrandsCoverage(ChromosomeRange((len-1) * 200, chr1.length, chr1))
            .toDouble()
    return cover
}
fun getGC(chr1: Chromosome): DoubleArray {
    val len = (chr1.length - 1) / 200 + 1
    val seq: TwoBitSequence = chr1.sequence
    val GCcontent = DoubleArray(len)
    for (i in 0 until len - 1) {
        GCcontent[i] = seq.substring(i*200, (i + 1)*200).count { it == 'c' || it == 'g' }.toDouble()/200
    }
    GCcontent[len - 1] = seq
            .substring((len-1)*200, seq.length)
            .count { it == 'c'|| it == 'g' }
            .toDouble()/( seq.length - (len-1)*200)
    return GCcontent
}

fun getMappability(chr1: Chromosome, path_mappability: Path): DoubleArray {
    val mapSummary = BigWigFile
            .read(path_mappability)
            .summarize(chr1.name, 0, chr1.length - chr1.length%200, numBins = (chr1.length - 1)/200)
    val result = DoubleArray (mapSummary.size + 1) {
        if (it < mapSummary.size) mapSummary[it].sum/200 else 1.0}
    result[mapSummary.size] = BigWigFile
            .read(path_mappability)
            .summarize(chr1.name, chr1.length - chr1.length%200, 0)[0].sum / chr1.length%200
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

fun fitK4Me3Bam(
        dirIn: String,
        dirOut: String,
        fileMe: String,
        fileInput: String,
        bic: DoubleArray,
        aic: DoubleArray,
        pval: Array<DoubleArray>,
        index: Int) {
    println("Start $fileMe")
    val path_me = Paths.get("$dirIn$fileMe")
    val path_input = Paths.get("$dirIn$fileInput")
    val path_mappability = Paths.get("/home/elena.kartysheva/Documents/test_data/1/wgEncodeCrgMapabilityAlign36mer.bigWig")
    val genomeQuery = GenomeQuery(Genome["hg19"])
    val readsQueryMe = ReadsQuery(genomeQuery, path_me, true)
    val readsQueryInput = ReadsQuery(genomeQuery, path_input, true)
    val coverageMe = readsQueryMe.get()
    val coverageInput = readsQueryInput.get()

    val chrList: List<Chromosome> = (1..23).map { if (it < 23) genomeQuery["chr$it"]!! else genomeQuery["chrx"]!!}

    val coverLength = chrList.map { it.length/200 + 1 }.sum()
    val coverMe = IntArray (coverLength)
    val coverInput = DoubleArray (coverLength)
    val GCcontent = DoubleArray (coverLength)
    val mappability = DoubleArray (coverLength)
    var prevIdx = 0
    chrList.forEach {
        arraycopy(getIntCover(it, coverageMe), 0, coverMe, prevIdx, it.length/200 + 1)
        arraycopy(getDoubleCover(it, coverageInput), 0, coverInput, prevIdx, it.length/200 + 1)
        arraycopy(getGC(it), 0, GCcontent, prevIdx, it.length/200 + 1)
        arraycopy(getMappability(it, path_mappability), 0, mappability, prevIdx, it.length/200 + 1)
        prevIdx += (it.length/200 + 1)
    }

    val covar = DataFrame()
            .with("y", coverMe)
            .with("x1", coverInput)
            .with("x2", GCcontent)
            .with("x3", mappability)
            .with("x4", DoubleArray(GCcontent.size) {GCcontent[it]*GCcontent[it]})
    Logs.addConsoleAppender(Level.DEBUG)
    val yaMix = ZeroPoissonMixture(
            doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
            listOf("x1", "x2", "x3", "x4"),
            arrayOf(doubleArrayOf(0.0, 0.0, 0.0, 0.0, 0.0), doubleArrayOf(1.0, 0.0, 0.0, 0.0, 0.0)))
    yaMix.fit(Preprocessed.of(covar), 1e-3, 2)
    bic[index] = yaMix.BIC(covar)
    aic[index] = yaMix.AIC(covar)
    yaMix.save(Paths.get("$dirOut$fileMe.json"))
    pval[index] = yaMix.Ftest(covar,
            0,
            Array2DRowRealMatrix(doubleArrayOf(0.0, 0.0, 1.0, 0.0, 1.0)),
            ArrayRealVector(doubleArrayOf(0.0)))
    println("P-value: ${pval[index].toList()}")
    println("Done $path_me")
}
