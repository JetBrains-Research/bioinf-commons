package org.jetbrains.bio.statistics.emission
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
import org.apache.log4j.Level
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
import java.nio.file.Paths
import java.util.function.IntPredicate
import java.nio.file.Path
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
        // needed here to add intercept to X
        val wlr = WLSMultipleLinearRegression()
        wlr.newSampleData(DoubleArray(x.size), x, DoubleArray(x.size))
        val X = wlr.getx()
        val y = ArrayRealVector(df.sliceAsInt(df.labels[d]).map { it.toDouble() }.toDoubleArray(), false)
        val iterMax = 5
        val tol = 1e-8
        var beta0: RealVector = ArrayRealVector(regressionCoefficients, false)
        var beta1: RealVector = ArrayRealVector(regressionCoefficients, false)
        for (i in 0 until iterMax) {
            val eta = X.operate(beta0)
            val countedLink = eta.map { link(it) }
            val countedLinkDerivative = eta.map { linkDerivative(it) }
            val z = DoubleArray(eta.dimension)
            {i -> eta.getEntry(i) + (y.getEntry(i) - countedLink.getEntry(i)) / countedLinkDerivative.getEntry(i)}
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
fun getIntCover(chr1: Chromosome, path: Path, genomeQuery: GenomeQuery): IntArray {
    val readsQuery = ReadsQuery(genomeQuery, path, false)
    val coverage = readsQuery.get()
    val len = chr1.length / 200 + 1
    val cover = IntArray(len)
    for (i in 0 until len - 1) {
        cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * 200, (i + 1) * 200, chr1))
    }
    cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * 200, chr1.length, chr1))
    return cover
}
fun getDoubleCover(chr1: Chromosome, path: Path, genomeQuery: GenomeQuery): DoubleArray {
    val readsQuery = ReadsQuery(genomeQuery, path, false)
    val coverage = readsQuery.get()
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
fun main(args: Array<String>) {
    //1
    //var covar = DataFrame()
    //        .with("y", IntArray(1000000))
    //        .with("x1", DoubleArray(1000000) { Random.nextDouble(0.5, 1.0) })
    //.with("x2", DoubleArray(1000000) { Random.nextDouble(0.0, 1.0) })
    /*
    //проверка, что с внешними весами все еще работает
    var regrES = PoissonRegressionEmissionScheme(listOf("x1", "x2"), doubleArrayOf(4.0, -2.0))
    val pred = IntPredicate {true}
    regrES.sample(covar, 0, pred)
    println("Update")
    regrES.update(covar, 0, DoubleArray(1000000, {1.0}).asF64Array())
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
    //H3K4me3
    val path_me = Paths.get("/home/elena.kartysheva/Documents/test_data/1/GSM409308_UCSD.H3K4me3.bam")
    val path_input = Paths.get("/home/elena.kartysheva/Documents/test_data/1/GSM605333_UCSD.H1.Input.LL-H1-I1.bed.gz")
    val genomeQuery = GenomeQuery(Genome["hg19"])
    var chr1 = genomeQuery["chr1"]!!

    var coverMe = getIntCover(chr1, path_me, genomeQuery)
    var coverInput = getDoubleCover(chr1, path_input, genomeQuery)
    var GCcontent = getGC(chr1)

    for(i in 2..22){
        println("chr$i")
        chr1 = genomeQuery["chr$i"]!!
        coverMe += getIntCover(chr1, path_me, genomeQuery)
        coverInput += getDoubleCover(chr1, path_input, genomeQuery)
        GCcontent += getGC(chr1)
    }

    chr1 = genomeQuery["chrX"]!!
    coverMe += getIntCover(chr1, path_me, genomeQuery)
    coverInput += getDoubleCover(chr1, path_input, genomeQuery)
    GCcontent += getGC(chr1)


    val covar = DataFrame()
            .with("y", coverMe)
            .with("x1", coverInput)
            .with("x2", GCcontent)
    Logs.addConsoleAppender(Level.DEBUG)
    val yaMix = ZeroPoissonMixture(
            doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
            listOf("x1", "x2"),
            arrayOf(doubleArrayOf(0.0, 1.5, 1.0), doubleArrayOf(1.0, 0.5, 0.0)))
    val start = System.currentTimeMillis()
    yaMix.fit(Preprocessed.of(covar), 1e-3, 50)
    println("Fit time: ${(System.currentTimeMillis() - start)}")
    println("weights: ${yaMix.weights}")
}
