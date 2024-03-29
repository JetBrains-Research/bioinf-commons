package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

interface RegionsMetric {
    val column: String

    fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double
    fun calcMetric(ra: RangesList, rb: RangesList): Double
    fun calcRegions(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): List<Location>
}