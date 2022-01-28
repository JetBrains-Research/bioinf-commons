package org.jetbrains.bio.statistics.model

import org.jetbrains.bio.dataframe.DataFrame

/**
 * Abstract iteration context for a statistical model.
 *
 * @author Sergei Lebedev
 * @since 13/10/14
 */
abstract class IterationContext(
    protected val numStates: Int,
    protected val df: DataFrame
) {
    open fun iterate() {
        refill()
        expect()
    }

    abstract fun expect()

    /**
     * Prepares "refilled" fields for the next iteration.
     */
    abstract fun refill()
}
