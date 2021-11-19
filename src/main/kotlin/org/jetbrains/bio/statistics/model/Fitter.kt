package org.jetbrains.bio.statistics.model

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.util.MultitaskProgress
import org.slf4j.LoggerFactory

/**
 * A fitter encapsulates the initialization logic for the classification model.
 *
 * @author Sergei Lebedev
 * @since 17/06/14
 */
interface Fitter<out Model : ClassificationModel> {
    /**
     * Constructs an initial guess for a classification model.
     *
     * @param preprocessed sample.
     * @param threshold convergence threshold (if applicable).
     * @param maxIter an upper bound on fitting iterations (if applicable).
     * @param attempt multistart attempt number (if applicable).
     * Should be used in guess to create different initial parameters.
     * @return guessed classification model.
     */
    fun guess(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIter: Int, attempt: Int
    ): Model

    fun guess(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIter: Int, attempt: Int
    ): Model = guess(preprocessed.first(), title, threshold, maxIter, attempt)

    fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIter: Int,
        attempt: Int = 0
    ): Model = fit(listOf(preprocessed), title, threshold, maxIter, attempt)

    fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String, threshold: Double, maxIter: Int,
        attempt: Int = 0
    ): Model {
        require(threshold > 0) { "threshold $threshold must be >0" }
        require(maxIter > 0) { "maximum number of iterations $maxIter must be >0" }

        val model = guess(preprocessed, title, threshold, maxIter, attempt)
        MultitaskProgress.addTask(title, maxIter.toLong())
        model.fit(preprocessed, title, threshold, maxIter)
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
            require(multiStarts > 1) { "number of starts $multiStarts must be >1" }
        }

        override fun fit(
            preprocessed: Preprocessed<DataFrame>, title: String,
            threshold: Double, maxIter: Int, attempt: Int
        ): Model {
            require(attempt == 0) {
                "cyclic multistart is not allowed"
            }
            require(maxIter > multiStartIter) {
                "maximum number of iterations $maxIter must be > multistart $multiStartIter"
            }
            val msModel = (0 until multiStarts).map {
                LOG.info("Multistart ${it + 1}/$multiStarts: $title")
                super.fit(preprocessed, "Multistart ${it + 1}/$multiStarts: $title", threshold, multiStartIter, it)
            }.maxByOrNull { it.logLikelihood(preprocessed) }!!
            LOG.info("Multistart done: $title")
            MultitaskProgress.addTask(title, maxIter.toLong())
            msModel.fit(preprocessed, title, threshold, maxIter)
            MultitaskProgress.finishTask(title)
            return msModel
        }

        override fun fit(
            preprocessed: List<Preprocessed<DataFrame>>, title: String,
            threshold: Double, maxIter: Int, attempt: Int
        ): Model {
            require(attempt == 0) {
                "cyclic multistart is not allowed"
            }
            require(maxIter > multiStartIter) {
                "maximum number of iterations $maxIter must be > multistart $multiStartIter"
            }
            val msModel = (0 until multiStarts).map {
                LOG.info("Multistart ${it + 1}/$multiStarts: $title")
                super.fit(preprocessed, "Multistart ${it + 1}/$multiStarts: $title", threshold, multiStartIter, it)
            }.maxByOrNull { m -> preprocessed.map { m.logLikelihood(it) }.sum() }!!
            LOG.info("Multistart done: $title")
            MultitaskProgress.addTask(title, maxIter.toLong())
            msModel.fit(preprocessed, title, threshold, maxIter)
            MultitaskProgress.finishTask(title)
            return msModel
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(ClassificationModel::class.java)
        const val TITLE = "unknown"
        const val THRESHOLD = 0.1
        const val MAX_ITERATIONS = 100
        const val MULTISTARTS = 5
        const val MULTISTART_ITERATIONS = 5
    }
}
