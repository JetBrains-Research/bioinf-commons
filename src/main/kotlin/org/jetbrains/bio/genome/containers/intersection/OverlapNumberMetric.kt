package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

/**
 * ~ bedtools intersect -wa -a .. -b .. | uniq | wc -l
 */
class OverlapNumberMetric(val aSetFlankedBothSides: Int = 0) : RegionsMetric {
    override val column =
        "overlap" + if (aSetFlankedBothSides == 0) "" else "_flnk_$aSetFlankedBothSides"

    override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double =
        a.calcAdditiveMetricDouble(b) { ra, rb -> calcMetric(ra, rb) }

    override fun calcMetric(ra: RangesList, rb: RangesList): Double =
        ra.overlapRangesNumber(rb, flankBothSides = aSetFlankedBothSides).toDouble()

    override fun calcRegions(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>) =
        a.calcMarkedLocations(b) { ra, rb -> ra.overlapRanges(rb, aSetFlankedBothSides) }
}