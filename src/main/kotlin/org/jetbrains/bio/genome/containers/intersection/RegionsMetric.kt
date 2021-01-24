package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList

interface RegionsMetric {
    val column: String

    fun calcMetric(a: LocationsList<out RangesList>, b: LocationsList<out RangesList>): Double
}