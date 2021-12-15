package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.argSort
import kotlin.math.min


object BenjaminiHochberg {
    /**
     * Applies Benjamini-Hochberg correction to the P-values.
     *
     * See [Fdr.qvalidate] for log version for posterior error probabilities.
     */
    fun adjust(ps: F64Array): F64Array {
        val m = ps.size

        // Sort Ps in ascending order and remember the inverse
        // permutation, so that we can return Q-values in the order
        // corresponding to the original Ps.
        val sorted = ps.argSort(reverse = true)
        val original = IntArray(m)
        for (k in 0 until m) {
            original[sorted[k]] = k
        }

        val qvalues  = F64Array(m)
        for (k in 0 until m) {
            qvalues[k] = min(1.0, ps[sorted[k]] * m / (m - k))
        }

        for (k in 1 until m) {
            qvalues[k] = min(qvalues[k], qvalues[k - 1])
        }

        qvalues.reorder(original)
        return qvalues
    }
}