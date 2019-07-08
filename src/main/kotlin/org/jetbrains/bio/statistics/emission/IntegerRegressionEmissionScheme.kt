package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
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
import java.io.*
import java.nio.file.Paths
import java.util.function.IntPredicate
import kotlin.random.Random

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

        val x = Array2DRowRealMatrix(covariateLabels.map { df.sliceAsDouble(it) }.toTypedArray(), false)
                .transpose()
                .data
/*
        val x = Array<DoubleArray> (covariateLabels.size + 1)
        { if (it == 0)
            DoubleArray(df.rowsNumber) {1.0}
        else
            df.sliceAsDouble(covariateLabels[it-1])}
        val X = Array2DRowRealMatrix(x, false).transpose() */

        /*
        val x = Array<DoubleArray> (df.rowsNumber) { DoubleArray(covariateLabels.size) {1.0} }

         covariateLabels.forEachIndexed { index, label ->
             (0 until df.rowsNumber).forEach { i ->
                 x[i][index] = df.getAsDouble(i, label)
             }
         } */



        //val xx = covariateLabels.map { df.sliceAsDouble(it) }
        // needed here to add intercept to X
        val wlr = WLSMultipleLinearRegression()
        wlr.newSampleData(DoubleArray(x.size), x, DoubleArray(x.size))
        val X = wlr.getx()
        val yInt = df.sliceAsInt(df.labels[d])
        val y = DoubleArray (yInt.size) {yInt[it].toDouble()}
        df.sliceAsInt(df.labels[d])
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
            wlr.newSampleData(z, x, W)
            beta1 = wlr.calculateBeta()
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
fun getIntCover(chr1: Chromosome, coverage: Coverage): IntArray {
    val len = chr1.length / 200 + 1
    val cover = IntArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * 200, (i + 1) * 200, chr1))
    }
    cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * 200, chr1.length, chr1))
    return cover
}
fun getDoubleCover(chr1: Chromosome, coverage: Coverage): DoubleArray {
    val len = chr1.length / 200 + 1
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
    val len = chr1.length / 200 + 1
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

fun fitK4Me3Bam(dirIn: String, dirOut: String, fileMe: String, fileInput: String) {
    println("Start $fileMe")
    val path_me = Paths.get("$dirIn$fileMe")
    val path_input = Paths.get("$dirIn$fileInput")
    val genomeQuery = GenomeQuery(Genome["hg19"])
    val readsQueryMe = ReadsQuery(genomeQuery, path_me, true)
    val readsQueryInput = ReadsQuery(genomeQuery, path_input, true)
    val coverageMe = readsQueryMe.get()
    val coverageInput = readsQueryInput.get()

    val chrList: List<Chromosome> = (1..23).map { if (it < 23) genomeQuery["chr$it"]!! else genomeQuery["chrx"]!!}

    val coverMe = chrList.flatMap { getIntCover(it, coverageMe).toList() }.toIntArray()
    val coverInput = chrList.flatMap { getDoubleCover(it, coverageInput).toList() }.toDoubleArray()
    val GCcontent = chrList.flatMap { getGC(it).toList() }.toDoubleArray()

    val covar = DataFrame()
            .with("y", coverMe)
            .with("x1", coverInput)
            .with("x2", GCcontent)
    Logs.addConsoleAppender(Level.DEBUG)
    val yaMix = ZeroPoissonMixture(
            doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
            listOf("x1", "x2"),
            arrayOf(doubleArrayOf(0.0, 0.0, 0.0), doubleArrayOf(1.0, 0.0, 0.0)))
    yaMix.fit(Preprocessed.of(covar), 1e-3, 30)
    yaMix.save(Paths.get("$dirOut$fileMe.json"))
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

    //обучение на 2х параметрах

    val directIn = "/mnt/stripe/bio/experiments/aging/chipseq/k4me3/k4me3_20vs20_bams/"

    File(directIn)
            .list()
            .filter { it.endsWith(".bam") && !it.contains("unique") && it != "input.bam" }
            .subList(20, 40)
            .forEach { fitK4Me3Bam(
                    directIn,
                    "/home/elena.kartysheva/Documents/aging_json_40/",
                    it,
                    "input.bam") }

}
