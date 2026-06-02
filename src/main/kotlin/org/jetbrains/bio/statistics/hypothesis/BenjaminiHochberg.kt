package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.argSort
import kotlin.math.ln
import kotlin.math.min


object BenjaminiHochberg {
    /**
     * Applies Benjamini-Hochberg correction to the P-values.
     *
     * See [adjustLogPValues] for the log-space version, and [Fdr.qvalidatePEPs]
     * for the variant operating on posterior error probabilities.
     */
    fun adjustPValues(ps: F64Array): F64Array {
        val m = ps.length
        // Q-value for ascending rank `rank`: min(1, m * p / rank).
        return adjust(ps) { p, rank -> min(1.0, p * m / rank) }
    }

    /**
     * Applies the Benjamini-Hochberg FDR correction to log P-values and
     * returns the resulting Q-values.
     *
     * The log-space equivalent of [adjustPValues], numerically stable for the
     * extremely small P-values that arise at genome scale.
     *
     * See [Fdr.qvalidatePEPs] for the variant operating on posterior error probabilities.
     *
     * @param logPs log P-values.
     * @param logResults if `true`, returns log Q-values; otherwise Q-values.
     */
    fun adjustLogPValues(logPs: F64Array, logResults: Boolean = false): F64Array {
        val m = logPs.length
        val logM = ln(m.toDouble())
        // Q-value for ascending rank `rank`: min(1, m * p / rank), in log space, capped at log(1) = 0.
        val qvalues = adjust(logPs) { logP, rank -> min(0.0, logP + logM - ln(rank.toDouble())) }
        if (!logResults) {
            qvalues.expInPlace()
        }
        return qvalues
    }

    /**
     * Shared Benjamini-Hochberg machinery: sort the values in descending order,
     * compute the raw Q-value for each rank via [rawQValue], enforce monotonicity
     * from the largest value down to the smallest, and restore the original order.
     *
     * @param rawQValue maps a sorted value and its ascending rank (1..m) to the raw Q-value.
     */
    private inline fun adjust(values: F64Array, rawQValue: (value: Double, rank: Int) -> Double): F64Array {
        val m = values.length
        // Sort in descending order and remember the inverse permutation, so that
        // we can return Q-values in the order corresponding to the original values.
        val sorted = values.argSort(reverse = true)
        val original = IntArray(m)
        for (k in 0 until m) {
            original[sorted[k]] = k
        }
        val qvalues = F64Array(m)
        for (k in 0 until m) {
            // Ascending rank of the k-th largest value is (m - k).
            qvalues[k] = rawQValue(values[sorted[k]], m - k)
        }
        // Enforce monotonicity from the largest value down to the smallest.
        for (k in 1 until m) {
            qvalues[k] = min(qvalues[k], qvalues[k - 1])
        }
        qvalues.reorder(original)
        return qvalues
    }

}