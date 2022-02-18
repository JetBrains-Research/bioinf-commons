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
     * @param attempt multistart attempt number (if applicable).
     * Should be used in guess to create different initial parameters.
     * @return guessed classification model.
     */
    fun guess(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIterations: Int, attempt: Int
    ): Model

    fun guess(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIterations: Int, attempt: Int
    ): Model = guess(preprocessed.first(), title, threshold, maxIterations, attempt)

    fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIterations: Int,
        attempt: Int = 0
    ): Model = fit(listOf(preprocessed), title, threshold, maxIterations, attempt)

    fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIterations: Int,
        attempt: Int = 0
    ): Model {
        require(threshold > 0) { "Threshold $threshold must be >0" }
        require(maxIterations > 0) { "Maximum number of iterations $maxIterations must be >0" }

        val model = guess(preprocessed, title, threshold, maxIterations, attempt)
        MultitaskProgress.addTask(title, maxIterations.toLong())
        model.fit(preprocessed, title, threshold, maxIterations)
        MultitaskProgress.finishTask(title)
        return model
    }

    /**
     * Returns a new fitter which runs the original one [multiStarts] times
     * and picks the model producing the highest log-likelihood.
     */
    fun multiStarted(
        multiStarts: Int = MULTISTARTS,
        multiStartIter: Int = MULTISTART_ITERATIONS
    ): Fitter<Model> = object : Fitter<Model> by this {
        init {
            require(multiStarts > 1) { "Number of starts $multiStarts must be >1" }
        }

        override fun fit(
            preprocessed: Preprocessed<DataFrame>, title: String,
            threshold: Double, maxIterations: Int, attempt: Int
        ): Model {
            require(attempt == 0) {
                "Cyclic multistart is not allowed"
            }
            val models = (0 until multiStarts).map {
                LOG.info("Multistart ${it + 1}/$multiStarts: $title")
                super.fit(preprocessed, "Multistart ${it + 1}/$multiStarts: $title", threshold, multiStartIter, it)
            }
            val mostLikelyModel = models.indices.maxByOrNull { models[it].logLikelihood(preprocessed) }!!
            LOG.info("Multistart done: $title, best model #${mostLikelyModel + 1}")

            MultitaskProgress.addTask(title, maxIterations.toLong())
            val result = models[mostLikelyModel]
            result.fit(preprocessed, title, threshold, maxIterations)
            MultitaskProgress.finishTask(title)
            return result
        }

        override fun fit(
            preprocessed: List<Preprocessed<DataFrame>>, title: String,
            threshold: Double, maxIterations: Int, attempt: Int
        ): Model {
            require(attempt == 0) {
                "Cyclic multistart is not allowed"
            }
            val models = (0 until multiStarts).map {
                LOG.info("Multistart ${it + 1}/$multiStarts: $title")
                super.fit(preprocessed, "Multistart ${it + 1}/$multiStarts: $title", threshold, multiStartIter, it)
            }
            val mostLikelyModel =
                models.indices.maxByOrNull { index -> preprocessed.sumOf { models[index].logLikelihood(it) } }!!
            LOG.info("Multistart done: $title, best model is #${mostLikelyModel + 1}")

            MultitaskProgress.addTask(title, maxIterations.toLong())
            val result = models[mostLikelyModel]
            result.fit(preprocessed, title, threshold, maxIterations)
            MultitaskProgress.finishTask(title)
            return result
        }
    }

    companion object {
        val LOG = LoggerFactory.getLogger(Fitter::class.java)
        const val TITLE = "unknown"
        const val THRESHOLD = 1.0
        const val MAX_ITERATIONS = 20
        const val MULTISTARTS = 5
        const val MULTISTART_ITERATIONS = 5
    }
}
