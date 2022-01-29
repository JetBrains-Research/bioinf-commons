package org.jetbrains.bio.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array

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
        for (state in 0 until numStates) {
            for (dimension in 0 until numDimensions) {
                res += getEmissionScheme(state, dimension).degreesOfFreedom
            }
        }
        return res
    }

    fun sample(numObservations: Int): DataFrame {
        val states = sampleStates(numObservations)
        var res = DataFrame()
        for (dimension in 0 until numDimensions) {
            val column = IntArray(numObservations)
            res = res.with("d$dimension", column)
            for (state in 0 until numStates) {
                getEmissionScheme(state, dimension).sample(res, dimension) { states[it] == state }
            }
        }
        return res.with("state", states)
    }

    override fun logProbability(state: Int, df: DataFrame, observation: Int): Double {
        var res = 0.0
        for (d in 0 until numDimensions) {
            res += getEmissionScheme(state, d).logProbability(df, observation, d)
        }
        return res
    }

    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        for (d in 0 until numDimensions) {
            for (state in 0 until numStates) {
                getEmissionScheme(state, d).update(df, d, gammas.V[state])
            }
        }
    }
}