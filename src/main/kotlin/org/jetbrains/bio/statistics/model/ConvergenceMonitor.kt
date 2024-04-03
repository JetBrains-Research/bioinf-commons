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
    protected val maxIterations: Int,
    private val level: Level
) {
    init {
        require(maxIterations > 1) { "maximum number of iterations must be greater than one, got $maxIterations" }
        require(0 < threshold && threshold < 1.0) { "threshold must be in 0..1, got $threshold" }
    }

    fun finish(model: ClassificationModel) = log(model)

    protected fun log(obj: Any) = LOG.log(level, "{$title} $obj")

    companion object {
        private val LOG = LoggerFactory.getLogger(ConvergenceMonitor::class.java)
    }
}