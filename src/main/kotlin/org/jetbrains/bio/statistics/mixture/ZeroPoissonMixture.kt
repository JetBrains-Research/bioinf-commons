package org.jetbrains.bio.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.emission.PoissonRegressionEmissionScheme
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate

/**
 * Mixture of 3 components:
 * 0 - zero-inflated component
 * 1 - LOW (poisson with small parameters)
 * 2 - HIGH (poisson with high parameters)
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
class ZeroPoissonMixture(weights: F64Array, covariateLabels: List<String>, regressionCoefficients: Array<DoubleArray>) :
        MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights) {

    private val components: List<EmissionScheme> = listOf(
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

    override fun sample(df: DataFrame, d: IntArray) {
        val states = sampleStates(df.rowsNumber)
        for (t in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, t).sample(df, d[t], IntPredicate { states[it] == i })
            }
        }
    }
}
