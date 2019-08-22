package org.jetbrains.bio.statistics.mixture
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.emission.PoissonRegressionEmissionScheme
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
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

    private val components: Array<EmissionScheme> = arrayOf(
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

    operator fun get(i: Int): EmissionScheme {
        return components[i]
    }

    operator fun set(i: Int, e: PoissonRegressionEmissionScheme) {
        if (i == 0) {
            throw IllegalArgumentException()
        } else {
            components[i] = e
        }
    }

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

    /*
    fun Ftest(df: DataFrame, d: Int, R: RealMatrix, r: RealVector): DoubleArray {
        return doubleArrayOf(
                (components[1] as PoissonRegressionEmissionScheme).Ftest(df, d, R, r),
                (components[2] as PoissonRegressionEmissionScheme).Ftest(df, d, R, r))
    }
    */

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
                        preprocessed.get(0).get().labels.drop(1),
                        arrayOf(
                                DoubleArray(preprocessed.get(0).get().columnsNumber) { 0.0 },
                                DoubleArray(preprocessed.get(0).get().columnsNumber){if (it == 0) 1.0 else 0.0}))
            }
        }
    }
}
