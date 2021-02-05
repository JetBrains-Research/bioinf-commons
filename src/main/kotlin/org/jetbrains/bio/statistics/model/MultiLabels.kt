package org.jetbrains.bio.statistics.model

import org.jetbrains.bio.dataframe.DataFrame

/**
 * A class with constants useful for writing multi-dimensional models.
 *
 * @author Sergei Lebedev
 * @since 16/09/14
 */
object MultiLabels {
    const val MAX_LABELS = 40

    /**
     * Labels should interned, because [DataFrame.getLabelIndex] works with reference equality for speed? reasons.
     */
    fun generate(prefix: String, numLabels: Int = MAX_LABELS): Array<String> {
        require(numLabels > 0) { "number of labels $numLabels must be >0" }
        return if (numLabels == 1) {
            arrayOf(prefix.intern())
        } else {
            Array(numLabels) { "${prefix}_${it + 1}".intern() }
        }
    }
}
