package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

class OverlapMetric(val aSetFlankedBothSides: Int = 0) : RegionsMetric {
    override val column =
        "overlap" + if (aSetFlankedBothSides == 0) "" else "_flnk_$aSetFlankedBothSides"

    override fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>) =
        a.calcAdditiveMetric(b) { ra, rb ->
            ra.overlapRangesNumber(rb, flankBothSides = aSetFlankedBothSides).toLong()
        }.toDouble()
}