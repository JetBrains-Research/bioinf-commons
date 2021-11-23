package org.jetbrains.bio.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate

/**
 * A hidden Markov model with multidimensional integer-valued
 * emissions.
 *
 * A multidimensional HMM works as follows: each combination of state
 * and dimension is associated with an emission scheme (i.e. an integer-valued
 * distribution family such as Poisson, negative binomial or singular) which
 * describes P(d_{k,t} | s_t = s) for any state s, dimension k and
 * observation t. Each emission scheme corresponds to a single pair (contrast
 * with constrained model).
 *
 * For the model to work one must implement exactly one method: getEmissionScheme,
 * which will return the emission scheme for the given state and dimension.
 *
 * Note that instead of an abstract getter we could potentially introduce an
 * array of abstract emission schemes. This, however, would make easy
 * serialization impossible.
 *
 * @author Alexey Dievsky
 * @author Oleg Shpynov
 * @date 01/04/15
 */
abstract class MLFreeHMM(
    numStates: Int,
    private val numDimensions: Int,
    priorProbabilities: F64Array = F64Array.stochastic(numStates),
    transitionProbabilities: F64Array = F64Array.stochastic(numStates, numStates)
) : MLAbstractHMM(numStates, priorProbabilities, transitionProbabilities) {

    protected abstract fun getEmissionScheme(i: Int, d: Int): EmissionScheme

    override fun degreesOfFreedom(): Int {
        var res = super.degreesOfFreedom()
        for (i in 0 until numStates) {
            for (d in 0 until numDimensions) {
                res += getEmissionScheme(i, d).degreesOfFreedom
            }
        }
        return res
    }

    fun sample(numObservations: Int): DataFrame {
        val states = sampleStates(numObservations)
        var res = DataFrame()
        for (d in 0 until numDimensions) {
            val column = IntArray(numObservations)
            res = res.with("d$d", column)
            for (i in 0 until numStates) {
                getEmissionScheme(i, d).sample(res, d) { states[it] == i }
            }
        }
        return res.with("state", states)
    }

    override fun logProbability(i: Int, df: DataFrame, t: Int): Double {
        var res = 0.0
        for (d in 0 until numDimensions) {
            res += getEmissionScheme(i, d).logProbability(df, t, d)
        }
        return res
    }

    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        for (d in 0 until numDimensions) {
            for (i in 0 until numStates) {
                getEmissionScheme(i, d).update(df, d, gammas.V[i])
            }
        }
    }
}