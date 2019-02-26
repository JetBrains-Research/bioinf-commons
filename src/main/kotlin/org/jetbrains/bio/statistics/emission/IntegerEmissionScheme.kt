package org.jetbrains.bio.statistics.emission

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.F64Array
import java.util.function.IntPredicate
import java.util.function.IntSupplier

/**
 * This interface is used in emission scheme-based HMMs and mixtures.
 *
 * @author Alexey Dievsky
 * @date 26/03/15
 */
interface IntegerEmissionScheme : EmissionScheme {
    fun sampler(): IntSupplier

    override fun sample(df: DataFrame, d: Int, fill: IntPredicate) {
        val supplier = sampler()
        val column = df.sliceAsInt(df.labels[d])
        for (i in 0 until df.rowsNumber) {
            if (fill.test(i)) {
                column[i] = supplier.asInt
            }
        }
    }

    fun logProbability(value: Int): Double

    override fun logProbability(df: DataFrame, t: Int, d: Int): Double {
        return NaNGuard(logProbability(df.getAsInt(t, df.labels[d])))
    }

    private fun NaNGuard(x: Double): Double = if (x.isNaN()) Double.NEGATIVE_INFINITY else x

    fun update(sample: IntArray, weights: F64Array)

    override fun update(df: DataFrame, d: Int, weights: F64Array) {
        update(df.sliceAsInt(df.labels[d]), weights)
    }
}
