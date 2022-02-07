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
     * @param maxIter an upper bound on fitting iterations (if applicable).
     * Should be used in guess to create different initial parameters.
     * @return guessed classification model.
     */
    fun guess(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIter: Int
    ): Model

    fun guess(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIter: Int
    ): Model = guess(preprocessed.first(), title, threshold, maxIter)

    fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIter: Int
    ): Model = fit(listOf(preprocessed), title, threshold, maxIter)

    fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIter: Int
    ): Model {
        require(threshold > 0) { "threshold $threshold must be >0" }
        require(maxIter > 0) { "maximum number of iterations $maxIter must be >0" }

        val model = guess(preprocessed, title, threshold, maxIter)
        MultitaskProgress.addTask(title, maxIter.toLong())
        model.fit(preprocessed, title, threshold, maxIter)
        MultitaskProgress.finishTask(title)
        return model
    }

    companion object {
        const val TITLE = "unknown"
        const val THRESHOLD = 1.0
        const val MAX_ITERATIONS = 50
    }
}
