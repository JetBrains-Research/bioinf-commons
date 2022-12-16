package org.jetbrains.bio.statistics.model

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.util.MultitaskProgress
import org.slf4j.LoggerFactory

/**
 * A fitter encapsulates the initialization logic for the classification model.
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 17/06/14
 */
interface Fitter<out Model : ClassificationModel> {
    /**
     * Constructs an initial guess for a classification model.
     *
     * @param preprocessed sample.
     * @param threshold convergence threshold (if applicable).
     * @param maxIterations an upper bound on fitting iterations (if applicable).
     * @return guessed classification model.
     */
    fun guess(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIterations: Int
    ): Model

    fun guess(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIterations: Int
    ): Model = guess(preprocessed.first(), title, threshold, maxIterations)

    fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIterations: Int
    ): Model = fit(listOf(preprocessed), title, threshold, maxIterations)

    fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIterations: Int
    ): Model {
        require(threshold > 0) { "Threshold $threshold must be >0" }
        require(maxIterations > 0) { "Maximum number of iterations $maxIterations must be >0" }

        val model = guess(preprocessed, title, threshold, maxIterations)
        MultitaskProgress.addTask(title, maxIterations.toLong())
        model.fit(preprocessed, title, threshold, maxIterations)
        MultitaskProgress.finishTask(title)
        return model
    }

    companion object {
        val LOG = LoggerFactory.getLogger(Fitter::class.java)
        const val TITLE = "unknown"
        const val THRESHOLD = 0.1
        const val MAX_ITERATIONS = 50
    }
}
