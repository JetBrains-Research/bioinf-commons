package org.jetbrains.bio.genome.containers.intersection

import joptsimple.OptionParser
import joptsimple.OptionSet
import joptsimple.util.RegexMatcher
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.RangesMergingList

enum class IntersectionDoubleMetric(
    val column: String,
    val metricFun: (LocationsMergingList, LocationsMergingList) -> Double
) {

    /**
     * We are using merging list, so overlap could be not correct if adjacent bin are merged in longer bins.
     * But for peak calling results should work ok
     */
    JACCARD("jaccard", { a, b ->
        val intersectionSize =
            a.apply(b, RangesMergingList::and).asLocationSequence().sumByDouble { it.length().toDouble() }
        val aSize = a.asLocationSequence().sumByDouble { it.length().toDouble() }
        val bSize = b.asLocationSequence().sumByDouble { it.length().toDouble() }
        intersectionSize / (aSize + bSize - intersectionSize)
    }),

    /**
     * We are using merging list, so overlap could be not correct if adjacent bin are merged in longer bins.
     * But for peak calling results should work ok
     */
    OVERLAP_FRACTION("overlap_fraction", { a, b ->
        IntersectionMetric.OVERLAP.metricFun(a, b).toDouble() / a.size
    });

    fun calcMetric(a: LocationsMergingList, b: LocationsMergingList) = metricFun(a, b)

    companion object {
        fun parse(options: OptionSet): IntersectionDoubleMetric {
            val metricName = options.valueOf("metric") as String
            val metric = values().find { it.column == metricName }
            checkNotNull(metric) { "Unsupported metric: $metricName" }
            return metric
        }

        fun OptionParser.acceptMetricArg() {
            accepts("metric", "Intersection metric type")
                .withRequiredArg()
                .ofType(String::class.java)
                .withValuesConvertedBy(
                    RegexMatcher.regex(
                        values().joinToString(separator = "|") { "(${it.column})" }
                    ))
                .defaultsTo(OVERLAP_FRACTION.column)
        }
    }
}