package org.jetbrains.bio.genome.containers.intersection

import joptsimple.OptionParser
import joptsimple.OptionSet
import joptsimple.util.RegexMatcher
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.RangesMergingList

object IntersectionMetrics {
    val OVERLAP = OverlapMetric(0)

    val INTERSECTION_NUMBER = object : RegionsMetric {
        // TODO: is it different from overlap?
        override val column = "intersection_number"

        override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>) =
            a.calcAdditiveMetric(b) { ra, rb ->
                require(ra is RangesMergingList) {
                    "Only merged range list supported here, but was ${ra::class.java.simpleName}"
                }
                require(rb is RangesMergingList) {
                    "Only merged range list supported here, but was ${rb::class.java.simpleName}"
                }
                ra.and(rb).size.toLong()
            }.toDouble()
    }

    /**
     * XXX We are using merging list, so overlap could be not correct if adjacent bin are merged in longer bins.
     *   But for peak calling results should work ok
     */
    val JACCARD = object : RegionsMetric {
        override val column = "jaccard"

        override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double {
            val intersectionSize = a.calcAdditiveMetricDouble(b) { ra, rb ->
                require(ra is RangesMergingList) {
                    "Only merged range list supported here, but was ${ra::class.java.simpleName}"
                }
                require(rb is RangesMergingList) {
                    "Only merged range list supported here, but was ${rb::class.java.simpleName}"
                }
                ra.and(rb).sumByDouble { it.length().toDouble() }
            }

            val aSize = a.asLocationSequence().sumByDouble { it.length().toDouble() }
            val bSize = b.asLocationSequence().sumByDouble { it.length().toDouble() }
            return intersectionSize / (aSize + bSize - intersectionSize)
        }
    }

    val OVERLAP_FRACTION = object : RegionsMetric {
        override val column = "overlap_fraction"

        override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double =
            OVERLAP.calcMetric(a, b) / a.size
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
                    listOf(OVERLAP, INTERSECTION_NUMBER, JACCARD, OVERLAP_FRACTION).joinToString(separator = "|") { "(${it.column})" }
                ))
            .defaultsTo(default)
    }
}