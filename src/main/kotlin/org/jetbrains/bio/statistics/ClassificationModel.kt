package org.jetbrains.bio.statistics

import com.google.common.primitives.Shorts
import com.google.gson.GsonBuilder
import com.google.gson.JsonParseException
import org.apache.commons.math3.distribution.AbstractIntegerDistribution
import org.apache.log4j.Logger
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.gson.F64ArrayTypeAdapter
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.statistics.gson.NotDirectlyDeserializable
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.viktor.F64Array
import java.io.IOException
import java.nio.file.Path
import java.util.function.IntToDoubleFunction


/**
 * An abstract classification model.
 *
 * @author Sergei Lebedev
 * @since 26/08/13
 */
interface ClassificationModel {
    /**
     * Returns the number of free parameters for the model
     */
    fun degreesOfFreedom(): Int

    /**
     * Fits a model to a given sample.
     *
     * @param preprocessed a sample to fit a model to.
     * @param title human-readable title for the sample, e.g. `"chr1"`.
     * @param threshold convergence threshold.
     * @param maxIter an upper bound on EM iterations.
     */
    fun fit(preprocessed: Preprocessed<DataFrame>, title: String,
            threshold: Double, maxIter: Int)

    fun fit(preprocessed: Preprocessed<DataFrame>,
            threshold: Double = Fitter.THRESHOLD,
            maxIter: Int = Fitter.MAX_ITERATIONS) {
        fit(preprocessed, Fitter.TITLE, threshold, maxIter)
    }

    fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String, threshold: Double, maxIter: Int) =
            fit(preprocessed.first(), title, threshold, maxIter)

    /**
     * Assigns state labels to a given sample (assumes the model is
     * already fitted.
     *
     * @param preprocessed a sample to assign state labels to.
     * @return a list, where i-th element is a label for the i-th observation.
     */
    fun predict(preprocessed: Preprocessed<DataFrame>): IntArray

    /**
     * Computes posterior state membership log-probabilities for a given
     * sample under the model.
     *
     * @param preprocessed a sample to evaluate.
     * @return a matrix where (i, t)-th element is \ln P(s_t = i|data, model).
     */
    fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array

    /**
     * Computes natural logarithm of the likelihood of a given sample under
     * the model.
     *
     * @param preprocessed a sample to evaluate.
     * @return log likelihood.
     */
    fun logLikelihood(preprocessed: Preprocessed<DataFrame>): Double

    /**
     * Serializes a model to JSON.
     *
     * CONTRACT: each Model class specifies "VERSION" int static field
     * describing current serialization format version.
     *
     * @param path to the resulting JSON file.
     */
    fun save(path: Path) {
        val gson = createGson()

        path.parent.createDirectories()
        try {
            path.bufferedWriter().use { gson.toJson(this, it) }
        } catch (e: StackOverflowError) {
            // Don't log exception, error message will be lost (removed from output) due to stacktrace size
            Logger.getRootLogger()
                    .error("Serialization StackOverflowError. Model ${javaClass.name}, path = ${path.toAbsolutePath()}")
            throw e
        }
    }

    companion object {

        /**
         * Loads a model from JSON file. While loading deserializer compares
         * version stored in JSON with version provided "VERSION" field and
         * informs whether they are different.
         *
         * CONTRACT: each Model class specifies "VERSION" int static field
         * describing current serialization format version.
         *
         * @param path to serialized model.
         * @throws IOException, JsonParseException if a given path is not readable.
         */
        @Suppress("unchecked_cast")
        @Throws(IOException::class, JsonParseException::class)
        fun <M : ClassificationModel> load(path: Path): M {
            val gson = createGson()
            return path.bufferedReader().use {
                val model = gson.fromJson(it, ClassificationModel::class.java) as M?
                checkNotNull(model) {
                    "failed to load model from $path"
                }
            }
        }

        /**
         * XXX please do not use this instance for serializing arbitrary objects to
         * JSON. It is specialized for serializing models.
         */
        private fun createGson() = GsonBuilder()
                .setPrettyPrinting()
                .setFieldNamingStrategy(GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY)
                .registerTypeAdapter(F64Array::class.java, F64ArrayTypeAdapter)
                .registerTypeAdapterFactory(NotDirectlyDeserializable.ADAPTER_FACTORY)
                .registerTypeAdapterFactory(
                        GSONUtil.classSpecificFactory(ClassificationModel::class.java) { gson, factory ->
                            GSONUtil.classAndVersionAdapter(gson, factory, "model.class.fqn", "model.class.format")
                        })
                .serializeSpecialFloatingPointValues()
                .create()

    }
}

/**
 * A fitter encapsulates the initialization logic for the classification
 * model. By convention the default fitter is accessible via
 * `Model.fitter` static method.
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
     * @return guessed classification model.
     */
    fun guess(preprocessed: Preprocessed<DataFrame>,
              title: String, threshold: Double, maxIter: Int, attempt: Int): Model

    fun guess(preprocessed: List<Preprocessed<DataFrame>>,
              title: String, threshold: Double, maxIter: Int, attempt: Int): Model =
            guess(preprocessed.first(), title, threshold, maxIter, attempt)

    fun fit(preprocessed: Preprocessed<DataFrame>,
            title: String = TITLE, threshold: Double = THRESHOLD,
            maxIter: Int = MAX_ITERATIONS,
            attempt: Int = 0): Model = fit(listOf(preprocessed), title, threshold, maxIter, attempt)

    fun fit(preprocessed: List<Preprocessed<DataFrame>>,
            title: String = TITLE, threshold: Double = THRESHOLD,
            maxIter: Int = MAX_ITERATIONS,
            attempt: Int = 0): Model {
        require(threshold > 0) { "threshold $threshold must be >0" }
        require(maxIter > 0) { "maximum number of iterations $maxIter must be >0" }

        val model = guess(preprocessed, title, threshold, maxIter, attempt)
        model.fit(preprocessed, title, threshold, maxIter)
        return model
    }

    /**
     * Returns a new fitter which runs the original one [multiStarts] times
     * and picks the model producing the highest log-likelihood.
     */
    fun multiStarted(
            multiStarts: Int = MULTISTARTS,
            multiStartIter: Int = MULTISTART_ITERATIONS): Fitter<Model> = object : Fitter<Model> by this {
        init {
            require(multiStarts > 1) { "number of starts $multiStarts must be >1" }
        }

        override fun fit(preprocessed: Preprocessed<DataFrame>, title: String,
                         threshold: Double, maxIter: Int, attempt: Int): Model {
            require(attempt == 0) {
                "cyclic multistart is not allowed"
            }
            require(maxIter > multiStartIter) {
                "maximum number of iterations $maxIter must be > multistart $multiStartIter"
            }
            val msModel = (0 until multiStarts).map {
                super.fit(preprocessed, "multistart $it: $title", threshold, multiStartIter, it)
            }.maxBy { it.logLikelihood(preprocessed) }!!
            msModel.fit(preprocessed, title, threshold, maxIter - multiStartIter)
            return msModel
        }

        override fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                         threshold: Double, maxIter: Int, attempt: Int): Model {
            require(attempt == 0) {
                "cyclic multistart is not allowed"
            }
            require(maxIter > multiStartIter) {
                "maximum number of iterations $maxIter must be > multistart $multiStartIter"
            }
            val msModel = (0 until multiStarts).map {
                super.fit(preprocessed, "multistart $it: $title", threshold, multiStartIter, it)
            }.maxBy { m -> preprocessed.map { m.logLikelihood(it) }.sum() }!!
            msModel.fit(preprocessed, title, threshold, maxIter - multiStartIter)
            return msModel
        }
    }

    companion object {
        const val TITLE = "unknown"
        const val THRESHOLD = 0.1
        const val MAX_ITERATIONS = 100
        const val MULTISTARTS = 5
        const val MULTISTART_ITERATIONS = 5
    }
}

/**
 * A common interface for hyper parameters of Bayesian models.
 *
 * @author Sergei Lebedev
 * @since 02/03/15
 */
interface HyperParameters<T : HyperParameters<T, InstanceModel>,
        InstanceModel : ClassificationModel> {
    /**
     * Calculates Kullback-Leibler divergence between the current parameter
     * distribution and the given one.
     *
     * @see org.jetbrains.bio.statistics.distribution.KullbackLeibler
     */
    fun divergence(other: T): Double

    fun copy(): T

    /**
     * Returns log density of the given frequentist model parameters
     * under the current parameter distribution.
     */
    fun logDensity(model: InstanceModel): Double {
        throw UnsupportedOperationException()
    }

    fun sample(): InstanceModel = throw UnsupportedOperationException()
}

/**
 * Abstract iteration context for a statistical model.
 *
 * @author Sergei Lebedev
 * @since 13/10/14
 */
abstract class IterationContext(protected val numStates: Int,
                                protected val df: DataFrame) {
    open fun iterate() {
        refill()
        expect()
    }

    protected abstract fun expect()

    /**
     * Prepares "refilled" fields for the next iteration.
     */
    abstract fun refill()
}

/**
 * A class with constants useful for writing multi-dimensional models.
 *
 * @author Sergei Lebedev
 * @since 16/09/14
 */
object MultiLabels {
    const val MAX_LABELS = 40

    fun generate(prefix: String, numLabels: Int = MAX_LABELS): Array<String> {
        require(numLabels > 0) { "number of labels $numLabels must be >0" }
        return if (numLabels == 1) {
            arrayOf(prefix)
        } else {
            Array(numLabels) { "${prefix}_${it + 1}" }
        }
    }
}

/**
 * A chain of data frame to data frame transformations used when building
 * a sampler.
 *
 * @author Sergei Lebedev
 * @since 12/09/14
 */
class SamplingChain private constructor(
        val states: IntArray,
        private val df: DataFrame) {

    private fun wrap(newDf: DataFrame) = SamplingChain(states, newDf)

    fun result() = df.with("state", states)

    fun map(f: (DataFrame) -> DataFrame): SamplingChain = wrap(f(df))

    // Note(lebedev): this is too common not to be factored in a method.
    fun binomial(kField: String, nField: String,
                 nDistribution: AbstractIntegerDistribution,
                 successProbability: IntToDoubleFunction): SamplingChain {
        val numObservations = states.size
        val ns = ShortArray(numObservations)
        val ks = ShortArray(numObservations)
        for (t in 0 until numObservations) {
            // We sample a zero-truncated version of the binomial distribution,
            // because n = 0 doesn't make sense for our models.
            ns[t] = Shorts.checkedCast((nDistribution.sample() + 1).toLong())
            ks[t] = Shorts.checkedCast(Sampling.sampleBinomial(
                    ns[t].toInt(),
                    successProbability.applyAsDouble(states[t])).toLong())
        }

        return wrap(df.with(nField, ns).with(kField, ks))
    }

    companion object {
        @JvmStatic
        fun start(states: IntArray): SamplingChain {
            require(states.isNotEmpty()) { "no data" }
            return SamplingChain(states, DataFrame())
        }
    }
}