package org.jetbrains.bio.statistics.model

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
