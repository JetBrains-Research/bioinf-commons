package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.argSort
import kotlin.math.min


object BenjaminiHochberg {
    /**
     * Applies Benjamini-Hochberg correction to the P-values.
     *
     * See [Fdr.qvalidate] for inplace version.
     */
    fun adjust(ps: F64Array): F64Array {
        val m = ps.size
        val sorted = ps.argSort(reverse = true)
        val original = IntArray(m)
        for (k in 0 until m) {
            original[sorted[k]] = k
        }

        val adjusted = F64Array(m)
        for (k in 0 until m) {
            adjusted[k] = min(1.0, ps[sorted[k]] * m / (m - k))
        }

        for (k in 1 until m) {
            adjusted[k] = min(adjusted[k], adjusted[k - 1])
        }

        adjusted.reorder(original)
        return adjusted
    }
}