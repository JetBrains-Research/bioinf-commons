package org.jetbrains.bio.statistics

import com.google.common.primitives.Ints
import java.util.*

/**
 * A container for internal parameters used by pairwise comparison models.
 *
 * @author Dmitry Groshev
 * @author Sergei Lebedev
 * @since 16/08/14
 */
class Constraints (
        /** Mapping from (dimension, state) pairs to equivalence classes.  */
        val constraintMap: Array<IntArray>,
        /**
         * Mapping from (dimension, 1D state) to the corresponding
         * equivalence classes.
         */
        val constraintMap1D: Array<IntArray>) {

    /**
     * Number of samples being compared.
     */
    val numDimensions = constraintMap.size

    /**
     * Number of equivalence classes. Should be the same as the number of
     * different elements in [constraintMap].
     */
    val numClasses: Int

    /**
     * Number of [DifferenceType] states used by the model. Should be
     * a square of [numStates1D].
     */
    val numStates = constraintMap.first().size

    /**
     * Number of one dimensional states used by the model. E.g. for
     * `MethylationDifference2` the number of 1D states is `2`.
     */
    val numStates1D = constraintMap1D.first().size

    init {
        require(constraintMap.isNotEmpty()) { "zero dimensions" }
        require(constraintMap1D.size == constraintMap.size) {
            "incompatible dimensions in constraint maps"
        }

        // Check that two constraint maps agree on the number of classes?
        numClasses = Arrays.stream(constraintMap).mapToInt { Ints.max(*it) }.max().asInt + 1
    }
}


/**
 * A common interface for pairwise comparison enums.
 *
 * @author Sergei Lebedev
 * @since 27/08/14
 */
interface DifferenceType<out State> {
    val isDifferent: Boolean

    fun component(index: Int) = when (index) {
        0 -> first()
        1 -> second()
        else -> throw IndexOutOfBoundsException()
    }

    fun first(): State
    fun second(): State
}

