package org.jetbrains.bio.gsea

import kotlin.math.min

data class TestedRegionStats(
    var countSetsWithMetricsAboveThr: Int = 0,
    var countSetsWithMetricsBelowThr: Int = 0,
    var simulationsNumber: Int = 0,
    var metricValueForInput: Long = 0L,
    val metricHist: IntHistogram = IntHistogram()
) {
    fun pvalue(hypAlt: PermutationAltHypothesis): Double {
        val pvalAbove = (countSetsWithMetricsAboveThr + 1).toDouble() / (simulationsNumber + 1)
        val pvalBelow = (countSetsWithMetricsBelowThr + 1).toDouble() / (simulationsNumber + 1)
        return when (hypAlt) {
            PermutationAltHypothesis.TWO_SIDED -> 2 * min(pvalAbove, pvalBelow)
            PermutationAltHypothesis.GREATER -> pvalAbove
            PermutationAltHypothesis.LESS -> pvalBelow
        }
    }
}