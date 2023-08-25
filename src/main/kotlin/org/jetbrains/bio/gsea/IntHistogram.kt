package org.jetbrains.bio.gsea

import gnu.trove.map.hash.TIntIntHashMap
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.viktor.KahanSum
import java.nio.file.Path
import kotlin.math.pow
import kotlin.math.sqrt

class IntHistogram {
    val data: TIntIntHashMap = TIntIntHashMap()

    fun increment(value: Int) {
        data.adjustOrPutValue(value, 1, 1)
    }

    fun mean(valuesNumber: Int = countValues()): Double {
        // sum(x_i) / n
        val metricMeanAcc = KahanSum()
        data.forEachEntry { metric, count ->
            metricMeanAcc += (count * metric).toDouble()
            true
        }
        return metricMeanAcc.result() / valuesNumber
    }

    fun countValues(): Int {
        var totalCount = 0
        data.forEachEntry { _, count ->
            totalCount += count
            true
        }
        return totalCount
    }

    operator fun plusAssign(other: IntHistogram) {
        other.data.forEachEntry { metric, count ->
            data.adjustOrPutValue(metric, count, count)
            true
        }
    }

    /**
     * Returns '-1' if list is empty
     */
    fun median(valuesNumber: Int = countValues()): Int {
        val medianIdx = valuesNumber / 2
        var idx = 0

        val keys = data.keys()
        keys.sort()

        for (metric in keys) {
            val count = data[metric]
            idx += count
            if (idx > medianIdx) {
                return metric
            }
        }
        return -1
    }

    fun stdev(valuesNumber: Int = countValues(), mean: Double = mean(valuesNumber)): Double {
        //listOf(1,3).stream().mapToInt { it }.average()
        //XXX: https://math.stackexchange.com/questions/857566/how-to-get-the-standard-deviation-of-a-given-histogram-image
        //XXX: online variance - https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance  and https://math.stackexchange.com/questions/198336/how-to-calculate-standard-deviation-with-streaming-inputs

        if (valuesNumber == 0) {
            return Double.NaN
        } else if (valuesNumber == 1) {
            return 0.0
        }
        //  sum((x_i - mean)^2) / (n - 1)
        val metricVarAcc = KahanSum()
        data.forEachEntry { metric, count ->
            metricVarAcc += count * (metric - mean).pow(2)
            true
        }
        return sqrt(metricVarAcc.result() / (valuesNumber - 1))
    }

    fun save(path: Path) {
        val keys = data.keys()
        keys.sort()
        DataFrame()
            .with("metric", keys)
            .with("count", keys.map { data[it] }.toIntArray())
            .save(path)
    }

    override fun toString() = data.toString()
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IntHistogram) return false

        if (data != other.data) return false

        return true
    }

    override fun hashCode(): Int = data.hashCode()

    companion object {
        fun create(data: IntArray): IntHistogram {
            val metricHist = IntHistogram()
            data.forEach { metricHist.increment(it) }
            return metricHist
        }
    }
}