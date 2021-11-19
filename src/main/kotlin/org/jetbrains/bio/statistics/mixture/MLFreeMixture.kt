package org.jetbrains.bio.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate

/**
 * @author Evgeny Kurbatsky
 * @since 04/07/15
 */
abstract class MLFreeMixture(
    numComponents: Int,
    protected val numDimensions: Int,
    weights: F64Array = F64Array.stochastic(numComponents)
) : MLAbstractMixture(numComponents, weights) {

    protected abstract fun getEmissionScheme(i: Int, d: Int): EmissionScheme

    override fun logProbability(i: Int, df: DataFrame, t: Int): Double {
        var res = 0.0
        for (d in 0 until numDimensions) {
            res += getEmissionScheme(i, d).logProbability(df, t, d)
        }
        return res
    }

    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        super.updateParameters(df, gammas)
        for (d in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, d).update(df, d, gammas.V[i])
            }
        }
    }

    open fun sample(numObservations: Int): DataFrame {
        val states = sampleStates(numObservations)
        var res = DataFrame()
        for (d in 0 until numDimensions) {
            val column = IntArray(numObservations)
            res = res.with("d$d", column)
            for (i in 0 until numComponents) {
                getEmissionScheme(i, d).sample(res, d, IntPredicate { states[it] == i })
            }
        }
        return res.with("state", states)
    }

    open fun sample(df: DataFrame, d: IntArray) {
        val states = sampleStates(df.rowsNumber)
        for (t in 0 until numDimensions) {
            for (i in 0 until numComponents) {
                getEmissionScheme(i, t).sample(df, d[t], IntPredicate { states[it] == i })
            }

        }
    }

    override fun degreesOfFreedom(): Int {
        var res = super.degreesOfFreedom()
        for (i in 0 until numComponents) {
            for (d in 0 until numDimensions) {
                res += getEmissionScheme(i, d).degreesOfFreedom
            }
        }
        return res
    }
}
