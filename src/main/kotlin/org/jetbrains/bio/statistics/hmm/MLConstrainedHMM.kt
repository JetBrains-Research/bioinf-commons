package org.jetbrains.bio.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.gson.NotDirectlyDeserializable
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor._I
import java.util.function.IntPredicate

/**
 * A hidden Markov model with multidimensional integer-valued emissions
 * and arbitrary state-dimension constraints.
 *
 * A constrained multidimensional HMM works as follows: a number of
 * emission schemes E is fixed, each emission scheme being an integer-valued
 * distribution family (e.g., Poisson, negative binomial or singular),
 * then each combination of state and dimension is assigned an emission
 * scheme which thus describes P(d_{k,t} | s_t = s) for any state s,
 * dimension k and observation t. In most cases it is implied that one
 * emission scheme doesn't cover more than one dimension; this makes
 * its learning rather straightforward. An exception can be made for
 * the schemes that are not supposed to be fitted, such as the singular
 * one.
 *
 * For the model to work one must implement exactly one method: [getEmissionScheme],
 * which will return the emission scheme with the given number. The
 * constructor must also provide the emission scheme assignment.
 *
 * Note that instead of an abstract getter we could potentially introduce
 * an array of abstract emission schemes. This, however, would make easy
 * serialization impossible.
 *
 * @author Alexey Dievsky
 * @date 26/03/15
 */
abstract class MLConstrainedHMM(
        /**
         * A map which assigns an emission scheme to each pair of state
         * and dimension:
         *
         * stateDimensionEmissionMap[state number][dimension number] = emission scheme number
         */
        protected val stateDimensionEmissionMap: Array<IntArray>,
        priorProbabilities: F64Array,
        transitionProbabilities: F64Array)
    : MLAbstractHMM(stateDimensionEmissionMap.size, priorProbabilities, transitionProbabilities),
        NotDirectlyDeserializable {

    protected val numDimensions: Int = stateDimensionEmissionMap[0].size
    protected val numEmissionSchemes: Int = generateEmissionCount(stateDimensionEmissionMap)

    @Transient
    private var emissionDimensionMap: IntArray =
            generateEmissionDimensionMap(stateDimensionEmissionMap)

    init {
        @Suppress("LeakingThis")
        updateTransients()
    }

    protected abstract fun getEmissionScheme(e: Int): EmissionScheme

    override fun logProbability(i: Int, df: DataFrame, t: Int): Double {
        var res = 0.0
        val dimensionEmissionMap = stateDimensionEmissionMap[i]
        for (d in 0 until numDimensions) {
            res += getEmissionScheme(dimensionEmissionMap[d]).logProbability(df, t, d)
        }

        return res
    }


    fun sample(numObservations: Int): DataFrame {
        val states = sampleStates(numObservations)
        var res = DataFrame()
        for (d in 0 until numDimensions) {
            val column = IntArray(numObservations)
            res = res.with(("d" + d).intern(), column)
            for (e in 0 until numEmissionSchemes) {
                if (emissionDimensionMap[e] != d) {
                    continue
                }

                getEmissionScheme(e).sample(res, d, IntPredicate {
                    stateDimensionEmissionMap[states[it]][d] == e
                })
            }
        }
        return res.with("state", states)
    }

    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        val weights = F64Array(df.rowsNumber)
        for (e in 0 until numEmissionSchemes) {
            val emissionScheme = getEmissionScheme(e)
            val d = emissionDimensionMap[e]
            for (i in 0 until numStates) {
                if (stateDimensionEmissionMap[i][d] != e) {
                    continue
                }

                weights += gammas.V[i]
            }

            emissionScheme.update(df, d, weights)
            weights.V[_I] = 0.0
        }
    }

    override fun updateTransients() {
        emissionDimensionMap = generateEmissionDimensionMap(stateDimensionEmissionMap)
    }

    override fun degreesOfFreedom(): Int {
        var res = super.degreesOfFreedom()
        for (e in 0 until numEmissionSchemes) {
            res += getEmissionScheme(e).degreesOfFreedom
        }
        return res
    }

    companion object {
        fun generateEmissionCount(stateDimensionEmissionMap: Array<IntArray>): Int {
            return stateDimensionEmissionMap.map { it.max()!! }.max()!! + 1
        }

        fun generateEmissionDimensionMap(stateDimensionEmissionMap: Array<IntArray>): IntArray {
            val res = IntArray(generateEmissionCount(stateDimensionEmissionMap))
            for (dimensionEmissionMap in stateDimensionEmissionMap) {
                for (d in dimensionEmissionMap.indices) {
                    res[dimensionEmissionMap[d]] = d
                }
            }
            return res
        }
    }
}
