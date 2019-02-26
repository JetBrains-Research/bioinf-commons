package org.jetbrains.bio.statistics.hypothesis

import com.google.common.collect.ImmutableSet
import org.jetbrains.bio.viktor.F64Array

/**
 * A null hypothesis for statistical testing.
 *
 * @author Sergei Lebedev
 * @since 13/02/15
 */
interface NullHypothesis<State> {
    val nullStates: Set<State>

    fun apply(logMemberships: Map<State, F64Array>): F64Array {
        return nullStates.map { logMemberships[it]!! }.reduceRight(F64Array::logAddExp)
    }

    companion object {
        fun <State> of(vararg nullStates: State): NullHypothesis<State> = of(nullStates.toSet())

        /** Creates a hypothesis for a set of null states. */
        fun <State> of(nullStates: Set<State>) = object : NullHypothesis<State> {
            override val nullStates = ImmutableSet.copyOf(nullStates)
        }
    }
}