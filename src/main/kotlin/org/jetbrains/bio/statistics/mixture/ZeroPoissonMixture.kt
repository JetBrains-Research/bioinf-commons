package org.jetbrains.bio.statistics.mixture
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.RealVector
import org.apache.log4j.Level
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.*
import org.jetbrains.bio.util.Logs
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Paths
import java.util.function.IntPredicate
import kotlin.math.ln

/**
 * Mixture of 3 components:
 * 0 - zero-inflated component
 * 1 - LOW (poisson with small parameters)
 * 2 - HIGH (poisson with high parameters)
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
class ZeroPoissonMixture(
        weights: F64Array, covariateLabels: List<String>, regressionCoefficients: Array<DoubleArray>
) : MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights) {

    val covariatesNum = covariateLabels.size

    private val components: List<EmissionScheme> = listOf(
            ConstantIntegerEmissionScheme(0),
            PoissonRegressionEmissionScheme(
                    covariateLabels = covariateLabels,
                    regressionCoefficients = regressionCoefficients[0]
            ),
            PoissonRegressionEmissionScheme(
                    covariateLabels = covariateLabels,
                    regressionCoefficients = regressionCoefficients[1]
            )
    )


    override fun getEmissionScheme(i: Int, d: Int): EmissionScheme = components[i]

    override fun sample(df: DataFrame, ds: IntArray) {
        val states = sampleStates(df.rowsNumber)
        for (d in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, d).sample(df, ds[d], IntPredicate { states[it] == i })
            }
        }
    }

    fun BIC(df: DataFrame): Double {
        return ln(df.rowsNumber.toDouble()) *(2*covariatesNum) - 2*this.logLikelihood(Preprocessed.of(df))
    }

    fun AIC(df: DataFrame): Double {
        return 2 *(2*covariatesNum) - 2*this.logLikelihood(Preprocessed.of(df))
    }

    fun Ftest(df: DataFrame, d: Int, R: RealMatrix, r: RealVector): DoubleArray {
        return doubleArrayOf(
                (components[1] as PoissonRegressionEmissionScheme).Ftest(df, d, R, r),
                (components[2] as PoissonRegressionEmissionScheme).Ftest(df, d, R, r))
    }

    companion object {
        @Transient @JvmField var VERSION = 1

        fun fitter() = object : Fitter<ZeroPoissonMixture> {
            override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                               threshold: Double, maxIter: Int, attempt: Int): ZeroPoissonMixture =
                    guess(listOf(preprocessed), title, threshold, maxIter, attempt)
            override fun guess(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                               threshold: Double, maxIter: Int, attempt: Int): ZeroPoissonMixture {
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                return ZeroPoissonMixture(
                        doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
                        listOf("x1", "x2", "x3", "x4"),
                        arrayOf(doubleArrayOf(0.0, 0.0, 0.0, 0.0, 0.0), doubleArrayOf(1.0, 0.0, 0.0, 0.0, 0.0)))
            }
        }

        fun fit(
                dirIn: String,
                dirOut: String,
                fileMe: String,
                fileInput: String,
                fileMappability: String,
                bic: DoubleArray,
                aic: DoubleArray,
                pval: Array<DoubleArray>,
                index: Int,
                bin: Int) {
            println("Start $fileMe")
            val path_me = Paths.get("$dirIn$fileMe")
            val path_input = Paths.get("$dirIn$fileInput")
            val path_mappability = Paths.get(fileMappability)
            val genomeQuery = GenomeQuery(Genome["hg19"])
            val readsQueryMe = ReadsQuery(genomeQuery, path_me, true)
            val readsQueryInput = ReadsQuery(genomeQuery, path_input, true)
            val coverageMe = readsQueryMe.get()
            val coverageInput = readsQueryInput.get()

            val chrList: List<Chromosome> = (1..23).map { if (it < 23) genomeQuery["chr$it"]!! else genomeQuery["chrx"]!!}

            val coverLength = chrList.map { it.length/bin + 1 }.sum()
            val coverMe = IntArray (coverLength)
            val coverInput = DoubleArray (coverLength)
            val GCcontent = DoubleArray (coverLength)
            val mappability = DoubleArray (coverLength)
            var prevIdx = 0
            chrList.forEach {
                System.arraycopy(getIntCover(it, coverageMe, bin), 0, coverMe, prevIdx, it.length / bin + 1)
                System.arraycopy(getDoubleCover(it, coverageInput, bin), 0, coverInput, prevIdx, it.length / bin + 1)
                System.arraycopy(getGC(it, bin), 0, GCcontent, prevIdx, it.length / bin + 1)
                System.arraycopy(getMappability(it, path_mappability, bin), 0, mappability, prevIdx, it.length / bin + 1)
                prevIdx += (it.length/bin + 1)
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
            println("Done $path_me")
        }
    }
}
