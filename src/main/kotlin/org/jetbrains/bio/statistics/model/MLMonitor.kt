package org.jetbrains.bio.statistics.model

import gnu.trove.list.array.TDoubleArrayList
import org.apache.commons.math3.exception.NotANumberException
import org.apache.commons.math3.util.Precision
import org.slf4j.event.Level
import kotlin.math.abs


/**
 * A monitor for models fitted via frequentist EM-algorithm.
 *
 * The run of EM converged if the difference in log-likelihoods is less
 * than [threshold] or if we've reached [maxIter].
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 03/09/14
 */
class MLMonitor(title: String, threshold: Double, maxIter: Int, level: Level = Level.DEBUG) :
    ConvergenceMonitor(title, threshold, maxIter, level) {

    private val logLikelihoods = TDoubleArrayList(maxIter)

    fun monitor(logLikelihood: Double): Boolean {
        if (logLikelihood.isNaN()) {
            throw NotANumberException()
        }

        logLikelihoods.add(logLikelihood)
        return converged()
    }

    private fun converged(): Boolean {
        val iter = logLikelihoods.size()
        val curr = logLikelihoods[iter - 1]
        if (iter <= 1 || iter == maxIter) {
            log("iteration: %03d  LL: %.6f".format(iter, curr))
            return iter == maxIter
        }

        val prev = logLikelihoods[iter - 2]
        val status = if (Precision.compareTo(curr, prev, threshold) < 0) " *" else ""
        log("iteration: %03d  LL: %.6f%s".format(iter, curr, status))

        return abs(curr - prev) < threshold
    }
}