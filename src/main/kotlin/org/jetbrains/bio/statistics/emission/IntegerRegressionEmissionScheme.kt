package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.distribution.FDistribution
import org.apache.commons.math3.linear.*
import org.apache.log4j.Level
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.coverage.SingleEndCoverage
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
import java.io.*
import java.lang.System.arraycopy
import java.nio.file.Paths
import java.util.function.IntPredicate
import kotlin.random.Random
import org.jetbrains.bio.big.BigWigFile
import java.nio.file.Path
import kotlin.math.ln
import kotlin.math.pow

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
            println("Iter $i")
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
        println(regressionCoefficients.toList())
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

fun main(args: Array<String>) {
    //1
    /*
    var covar = DataFrame()
           .with("y", IntArray(1000000))
           .with("x1", DoubleArray(1000000) { Random.nextDouble(0.5, 1.0) })
    .with("x2", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })

    //проверка, что с внешними весами все еще работает
    var regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(0.0, 4.0, -2.0))
    val pred = IntPredicate {true}
    regrES.sample(covar, 0, pred)
    println("Update")
    regrES.update(covar, 0, DoubleArray(1000000) {1.0}.asF64Array())
    print("Beta: ${regrES.regressionCoefficients.asList()}")
    */
    // MLFreeMixture
    /*
    Logs.addConsoleAppender(Level.DEBUG)
    var mix = ZeroPoissonMixture(
            doubleArrayOf(0.4, 0.5, 0.1).asF64Array(),
            listOf("x1"/*, "x2"*/),
            arrayOf(doubleArrayOf(1.0, 1.0), doubleArrayOf(2.0, 2.0)))
    var start = System.currentTimeMillis()
    mix.sample(covar, intArrayOf(0))
    println("Sample time: ${(System.currentTimeMillis() - start)}")
    var yaMix = ZeroPoissonMixture(
            doubleArrayOf(1/3.0, 1/3.0, 1/3.0).asF64Array(),
            listOf("x1"/*, "x2"*/),
            arrayOf(doubleArrayOf(0.0, 1.5), doubleArrayOf(1.6, 0.0)))
    start = System.currentTimeMillis()
    yaMix.fit(Preprocessed.of(covar), 1e-3, 12)
    println("Fit time: ${(System.currentTimeMillis() - start)}")
    */

    //verify F-test
    /*
    val model = PoissonRegressionEmissionScheme(listOf("x2", "x3", "x4", "x5"), doubleArrayOf(1.0, 0.0, 0.0, 0.0, 0.0))
    var dat = DataFrame()
            .with("y", intArrayOf(2, 6, 2, 8, 5, 6))
            .with("x2", doubleArrayOf(1.0, 5.0, 3.0, 8.0, 4.0, 6.0))
            .with("x3", doubleArrayOf(0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
            .with("x4", doubleArrayOf(55.2, 51.4, 47.2, 50.2, 49.0, 49.5))
            .with("x5", doubleArrayOf(3047.04, 2641.96, 2227.84, 2520.04, 2401.00, 2450.25))

    model.update(dat, 0, DoubleArray(6){1.0}.asF64Array())
    val pv = model.Ftest(dat, 0, Array2DRowRealMatrix(doubleArrayOf(0.0, 0.0, 1.0, 0.0, 1.0)),
            ArrayRealVector(doubleArrayOf(0.0)))
    println(pv)
    */


    val directIn = "/mnt/stripe/bio/experiments/aging/chipseq/k4me3/k4me3_20vs20_bams/"
    val bic = DoubleArray (40)
    val aic = DoubleArray (40)
    val pval = Array(40) { doubleArrayOf(1.0, 1.0)}
    val fileList = File(directIn)
            .list()
            .filter { it.endsWith(".bam") && !it.contains("unique") && it != "input.bam" }

    fileList.forEachIndexed {index, item -> fitK4Me3Bam(
            directIn,
            "/home/elena.kartysheva/Documents/aging_json_40/",
            item,
            "input.bam",
            bic,
            aic,
            pval,
            index) }

    var fileOut = "/home/elena.kartysheva/Documents/test_data/bic_with_GC^2.txt"
    var myfile = File(fileOut)

    myfile.printWriter().use { out ->
        bic.forEachIndexed { index, d ->
            out.println(fileList[index] + ": " + d.toString())
        }
    }

    fileOut = "/home/elena.kartysheva/Documents/test_data/aic__with_GC^2.txt"
    myfile = File(fileOut)

    myfile.printWriter().use { out ->
        aic.forEachIndexed { index, d ->
            out.println(fileList[index] + ": " + d.toString())
        }
    }

    fileOut = "/home/elena.kartysheva/Documents/test_data/pval.txt"
    myfile = File(fileOut)

    myfile.printWriter().use { out ->
        pval.forEachIndexed { index, d ->
            out.println(fileList[index] + ": " + d.toList())
        }
    }

    println("Wrote to file")
}
