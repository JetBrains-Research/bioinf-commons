package org.jetbrains.bio.statistics.mixture
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.RealVector
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.emission.PoissonRegressionEmissionScheme
import org.jetbrains.bio.viktor.F64Array
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
    }
}
