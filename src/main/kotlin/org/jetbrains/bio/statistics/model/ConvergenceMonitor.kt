package org.jetbrains.bio.statistics.model

import org.jetbrains.bio.util.log
import org.slf4j.LoggerFactory
import org.slf4j.event.Level

/**
 * A common class for model monitors.
 *
 * @author Sergei Lebedev
 * @since 05/09/14
 */
abstract class ConvergenceMonitor(
    protected val title: String,
    protected val threshold: Double,
    protected val maxIter: Int,
    private val level: Level
) {
    init {
        require(maxIter > 1) { "maximum number of iterations must be greater than one" }
        require(threshold >= 0) { "threshold must be >= 0" }
    }

    fun finish(model: ClassificationModel) = log(model)

    protected fun log(obj: Any) = LOG.log(level, "{$title} $obj")

    companion object {
        private val LOG = LoggerFactory.getLogger(ConvergenceMonitor::class.java)
    }
}