package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

/**
 * ~ bedtools intersect -wa -a .. -b .. | wc -l
 */
class IntersectionNumberMetric (val aSetFlankedBothSides: Int = 0) : RegionsMetric {
    override val column =
        "intersection_number" + if (aSetFlankedBothSides == 0) "" else "_flnk_$aSetFlankedBothSides"

    override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double =
        a.calcAdditiveMetricDouble(b) { ra, rb -> calcMetric(ra, rb) }

    override fun calcMetric(ra: RangesList, rb: RangesList): Double =
        ra.intersectRangesNumber(rb, flankBothSides = aSetFlankedBothSides).toDouble()
}