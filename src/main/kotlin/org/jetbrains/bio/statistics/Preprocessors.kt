package org.jetbrains.bio.statistics

import org.apache.commons.math3.util.CombinatoricsUtils
import org.jetbrains.bio.dataframe.DataFrame
import java.util.function.Supplier
import java.util.stream.IntStream

/**
 * A preprocessor augments a sample in a way required by the model. For
 * example it might pre-compute and cache expensive constants like
 * binomial coefficients or factorials.
 *
 * @author Sergei Lebedev
 * @since 17/06/14
 */
interface Preprocessor {
    /**
     * Returns an augmented version of the given `sample`. The original
     * sample is guaranteed to be unmodified.
     */
    fun apply(sample: DataFrame): Preprocessed<DataFrame>

    /**
     * Returns a composed preprocessor that first applies this preprocessor
     * to its input, and then applies the `after` preprocessor to the
     * result.
     */
    fun andThen(after: Preprocessor): Preprocessor {
        return Preprocessor { after.apply(apply(it).get()) }
    }

    companion object {
        inline operator fun invoke(
            crossinline
            block: (DataFrame) -> Preprocessed<DataFrame>
        ): Preprocessor {
            return object : Preprocessor {
                override fun apply(sample: DataFrame): Preprocessed<DataFrame> = block(sample)
            }
        }
    }
}

/**
 * A wrapper for a data frame preprocessed by [Preprocessor].
 *
 * @author Dmitry Groshev
 * @author Sergei Lebedev
 * @since 14/08/14
 */
class Preprocessed<T> private constructor(private val value: T) : Supplier<T> {
    override fun get() = value

    inline fun map(f: (T) -> T) = of(f(get()))

    companion object {
        @JvmStatic
        fun <T> of(value: T): Preprocessed<T> = Preprocessed(value)
    }
}

object Preprocessors {
    /** Returns a no-op preprocessor.  */
    @JvmStatic
    fun identity(): Preprocessor = Preprocessor { Preprocessed.of(it) }

    /**
     * Constructs a preprocessor which computes (n, k) from the given
     * columns of the data frame and stores the result in `toField`.
     */
    @JvmStatic
    fun binomialCoefficientLog(
        nField: String, kField: String,
        toField: String
    ): Preprocessor {
        return Preprocessor { df ->
            val rowsNumber = df.rowsNumber
            val logbcs = DoubleArray(rowsNumber)
            (0 until rowsNumber).forEach { i ->
                val k = df.getAsShort(i, kField).toInt()
                val n = df.getAsShort(i, nField).toInt()
                logbcs[i] = CombinatoricsUtils.binomialCoefficientLog(n, k)
            }
            Preprocessed.of(df.with(toField, logbcs))
        }
    }
}