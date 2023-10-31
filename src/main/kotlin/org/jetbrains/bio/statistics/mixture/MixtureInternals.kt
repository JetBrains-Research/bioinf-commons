package org.jetbrains.bio.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor._I

/**
 * Generic algorithms for fitting mixtures.
 *
 * @author Sergei Lebedev
 * @since 10/10/14
 */
object MixtureInternals {
    @JvmStatic
    fun predict(logGammas: F64Array): IntArray {
        return IntArray(logGammas.shape[1]) { logGammas.V[_I, it].argMax() }
    }

    @JvmStatic
    fun evaluate(
        logJointProbabilities: F64Array,
        logGammas: F64Array
    ) {
        logJointProbabilities.matrixTransposeTo(logGammas)
        logGammas.along(1).forEach(F64Array::logRescale)
    }

    @JvmStatic
    fun evaluateTransposed(
        logJointProbabilities: F64Array,
        logTransposedGammas: F64Array
    ) {
        logJointProbabilities.copyTo(logTransposedGammas)
        logTransposedGammas.along(0).forEach(F64Array::logRescale)
    }

    @JvmStatic
    fun logLikelihood(
        df: DataFrame,
        logJointProbabilities: F64Array
    ): Double {
        val v = F64Array.Viewer(logJointProbabilities)
        return (0 until df.rowsNumber).sumOf {
            v[it].logSumExp()
        }
    }
}

internal fun F64Array.matrixTransposeTo(other: F64Array) {
    check(nDim == 2) { "not a matrix" }
    check(other.nDim == 2) { "not a matrix" }
    check(shape[0] == other.shape[1] && shape[1] == other.shape[0]) { "matrices not conforming" }
    for (i in 0 until shape[1]) {
        other.V[i] = V[_I, i]
    }
}