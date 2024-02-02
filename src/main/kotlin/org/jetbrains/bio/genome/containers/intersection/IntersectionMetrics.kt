package org.jetbrains.bio.genome.containers.intersection

import joptsimple.OptionParser
import joptsimple.OptionSet
import joptsimple.util.RegexMatcher
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

object IntersectionMetrics {
    val OVERLAP = OverlapNumberMetric()

    val INTERSECTION_NUMBER = IntersectionNumberMetric()

    /**
     * XXX We are using merging list, so overlap could be not correct if adjacent bin are merged in longer bins.
     *   But for peak calling results should work ok
     */
    val JACCARD = object : RegionsMetric {
        override val column = "jaccard"

        override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double {
            val intersectionSize = a.calcAdditiveMetricDouble(b) { ra, rb -> calcMetric(ra, rb) }

            val aSize = a.rangeLists.sumOf { ranges -> ranges.sumOf { it.length().toDouble() } }
            val bSize = b.rangeLists.sumOf { ranges -> ranges.sumOf { it.length().toDouble() } }
            return intersectionSize / (aSize + bSize - intersectionSize)
        }

        override fun calcMetric(ra: RangesList, rb: RangesList): Double =
            ra.intersectRanges(rb).sumOf { it.length().toDouble() }

        //TODO implement to calculate proper regions
        override fun calcRegions(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>) =
            a.calcMarkedLocations(b) {ra , rb -> ra.intersectRanges(rb)}
    }

    val OVERLAP_FRACTION = object : RegionsMetric {
        override val column = "overlap_fraction"

        override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double =
            OVERLAP.calcMetric(a, b) / a.size

        override fun calcMetric(ra: RangesList, rb: RangesList): Double = OVERLAP.calcMetric(ra, rb)

        override fun calcRegions(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>) = OVERLAP.calcRegions(a, b)
    }

    fun parse(options: OptionSet) = (options.valueOf("metric") as String).let { text ->
        when (text) {
            OVERLAP.column -> OVERLAP
            INTERSECTION_NUMBER.column -> INTERSECTION_NUMBER
            JACCARD.column -> JACCARD
            OVERLAP_FRACTION.column -> OVERLAP_FRACTION
            else -> throw IllegalArgumentException("Not supported: $text")
        }
    }

    fun OptionParser.acceptMetricArg(default: String) {
        accepts("metric", "Intersection metric type for metric(a,b)")
            .withRequiredArg()
            .ofType(String::class.java)
            .withValuesConvertedBy(
                RegexMatcher.regex(
                    listOf(OVERLAP, INTERSECTION_NUMBER, JACCARD, OVERLAP_FRACTION)
                        .joinToString(separator = "|") { "(${it.column})" }
                ))
            .defaultsTo(default)
    }
}