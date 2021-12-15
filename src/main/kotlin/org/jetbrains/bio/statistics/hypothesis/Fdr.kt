package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.argSort
import org.jetbrains.bio.viktor.logAddExp
import kotlin.math.ln
import kotlin.math.min

/**
 * Procedures for estimating and controlling FDR for predictions for [ClassificationModel]
 *
 * @author Sergei Lebedev
 * @since 06/06/14
 */
class Fdr(private val alpha: Double) : Predictor {
    override fun predict(logNullMemberships: F64Array): BitterSet {
        return control(logNullMemberships, alpha)
    }

    companion object {

        /**
         * Constructs null membership log-probabilities and passes them to
         * [control].
         */
        fun <State> control(
            logMemberships: Map<State, F64Array>,
            nullStates: Set<State>, alpha: Double
        ): BitterSet {
            return control(NullHypothesis.of(nullStates).apply(logMemberships), alpha)
        }

        /**
         * Controls FDR at level `alpha`.
         *
         * @param logNullMemberships null membership log-probabilities.
         * @param alpha desired FDR.
         * @return a bit vector, where the i-th bit is 0 if the null was NOT
         *          rejected and 1 otherwise.
         *
         * References:
         *
         * Sun & Cai, "Oracle and Adaptive Compound Decision Rules for
         *             False Discovery Rate Control", ASA, 2007.
         * Sun & Cai, "Large scale multiple testing under dependence", JRSS, 2008.
         */
        fun control(logNullMemberships: F64Array, alpha: Double): BitterSet {
            val indices = logNullMemberships.argSort()
            val res = BitterSet(logNullMemberships.size)
            var logSum = Double.NEGATIVE_INFINITY
            var logFdr: Double
            var count = 0
            do {
                logSum = logSum logAddExp logNullMemberships[indices[count]]
                count++
                logFdr = logSum - ln(count.toDouble())
            } while (count < logNullMemberships.size && logFdr <= ln(alpha))

            // All hypotheses prior to and including 'count' are rejected. The
            // '-1' is to compensate the increment in the while loop above.
            for (t in 0..count - 2) {
                res.set(indices[t])
            }

            return res
        }

        /**
         * Estimates Q-values from given posterior error probabilities.
         *
         * Procedure is similar to [BenjaminiHochberg] but for log PEPs.
         *
         * References:
         *
         * Storey & Tibshirani, "Statistical significance for genomewide studies", PNAS, 2003.
         * Cheng & Zhu, "A classification approach to DNA methylation profiling with
         *               bisulfite next-generation sequencing", Bioinformatics, 2014.
         */
        fun qvalidate(logNullMemberships: F64Array): F64Array {
            val m = logNullMemberships.size
            // Sort PEPs in ascending order and remember the inverse
            // permutation, so that we can return Q-values in the order
            // corresponding to the original PEPs.
            val sorted = logNullMemberships.argSort()
            val original = IntArray(m)
            for (k in 0 until m) {
                original[sorted[k]] = k
            }
            val qvalues = logNullMemberships.copy()
            qvalues.reorder(sorted)
            var acc = Double.NEGATIVE_INFINITY
            for (k in 0 until m) {
                acc = acc logAddExp qvalues[k]
                qvalues[k] = acc - ln((k + 1).toDouble())
            }
            for (k in m - 2 downTo 0) {
                qvalues[k] = min(qvalues[k], qvalues[k + 1])
            }
            qvalues.reorder(original)
            qvalues.expInPlace()
            return qvalues
        }
    }
}