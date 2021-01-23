package org.jetbrains.bio.genome.containers.intersection

import joptsimple.OptionParser
import joptsimple.util.RegexMatcher
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.RangesMergingList

enum class IntersectionMetric(
        val column: String,
        val metricFun: (LocationsMergingList, LocationsMergingList) -> Long
) {
    OVERLAP("overlap_number", { a, b ->
        a.apply(b) { ra, rb -> ra.overlap(rb) }.size.toLong()
    }),

    INTERSECTION_NUMBER("intersection_number", { a, b ->
        a.apply(b, RangesMergingList::and).size.toLong()
    });

    fun calcMetric(a: LocationsMergingList, b: LocationsMergingList) = metricFun(a, b)

    companion object {
        fun parse(metricName: String): IntersectionMetric {
            val metric = values().find { it.column == metricName }
            checkNotNull(metric) { "Unsupported metric: $metricName" }
            return metric
        }

        fun OptionParser.acceptMetricArg() {
            accepts("metric", "Intersection metric type for metric(a,b)")
                    .withRequiredArg()
                    .ofType(String::class.java)
                    .withValuesConvertedBy(RegexMatcher.regex(
                        values().joinToString(separator = "|") { "(${it.column})" }
                    ))
                    .defaultsTo(OVERLAP.column)
        }
    }
}