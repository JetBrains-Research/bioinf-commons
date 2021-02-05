package org.jetbrains.bio.statistics.model

import com.google.gson.GsonBuilder
import com.google.gson.JsonParseException
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.gson.F64ArrayTypeAdapter
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.statistics.gson.NotDirectlyDeserializable
import org.jetbrains.bio.util.Logs
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.viktor.F64Array
import java.io.IOException
import java.nio.file.Path

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
    fun fit(preprocessed: Preprocessed<DataFrame>, title: String, threshold: Double, maxIter: Int)

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
            Logs.getRootLogger().error(
                    "Serialization StackOverflowError. Model ${javaClass.name}, path = ${path.toAbsolutePath()}")
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
