package org.jetbrains.bio.statistics

import gnu.trove.list.array.TDoubleArrayList
import org.apache.commons.math3.exception.NotANumberException
import org.apache.commons.math3.util.Precision
import org.apache.log4j.Level
import org.apache.log4j.Logger

/**
 * A common class for model monitors.
 *
 * @author Sergei Lebedev
 * @since 05/09/14
 */
abstract class ConvergenceMonitor(protected val title: String,
                                  protected val threshold: Double,
                                  protected val maxIter: Int,
                                  private val level: Level) {
    private val LOG = Logger.getLogger(ConvergenceMonitor::class.java)

    init {
        require(maxIter > 1) { "maximum number of iterations must be greater than one" }
        require(threshold >= 0) { "threshold must be >= 0" }
    }

    fun finish(model: ClassificationModel) = log(model)

    protected fun log(obj: Any) = LOG.log(level, "{$title} $obj")
}


/**
 * A monitor for models fitted via frequentist EM-algorithm.
 *
 * The run of EM converged if the difference in log-likelihoods is less
 * than [threshold] or if we've reached [maxIter].
 *
 * @author Sergei Lebedev
 * @since 03/09/14
 */
class MLMonitor(title: String, threshold: Double, maxIter: Int,
                level: Level = Level.DEBUG) :
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

        return Math.abs(curr - prev) < threshold
    }
}