package org.jetbrains.bio.genome.containers

import joptsimple.OptionParser
import joptsimple.util.RegexMatcher

enum class IntersectionMetric(
        val column: String,
        val metricFun: (LocationsMergingList, LocationsMergingList) -> Long
) {
    OVERLAP("overlap_number", { a, b ->
        a.apply(b, RangesMergingList::overlap).size.toLong()
    }),

    INTERSECTION_NUMBER("intersection_number", { a, b ->
        a.apply(b, RangesMergingList::and).size.toLong()
    });

    fun calcMetric(a: LocationsMergingList, b: LocationsMergingList) = metricFun(a, b)

    companion object {
        fun parse(metricName: String): IntersectionMetric {
            val metric = IntersectionMetric.values().find { it.column == metricName }
            checkNotNull(metric) { "Unsupported metric: $metricName" }
            return metric!!
        }

        fun OptionParser.acceptMetricArg() {
            accepts("metric", "Intersection metric type")
                    .withRequiredArg()
                    .ofType(String::class.java)
                    .withValuesConvertedBy(RegexMatcher.regex(
                            values().joinToString(separator = "|") { "(${it.column})" }
                    ))
                    .defaultsTo(OVERLAP.column)
        }
    }
}