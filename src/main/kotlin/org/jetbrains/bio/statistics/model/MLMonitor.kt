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
 * than [threshold] or if we've reached [maxIterations].
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 03/09/14
 */
class MLMonitor(title: String, threshold: Double, maxIterations: Int, level: Level = Level.DEBUG) :
    ConvergenceMonitor(title, threshold, maxIterations, level) {

    private val logLikelihoods = TDoubleArrayList(maxIterations)

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
        if (iter <= 1 || iter == maxIterations) {
            log("iteration: %03d  LL: %.2f".format(iter, curr))
            return iter == maxIterations
        }

        val prev = logLikelihoods[iter - 2]
        val result = abs((curr - prev) / (curr + prev)) < threshold
        val status = if (result) " *" else ""
        log("iteration: %03d  LL: %.2f%s".format(iter, curr, status))
        return result
    }
}