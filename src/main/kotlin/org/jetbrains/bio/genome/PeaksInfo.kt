package org.jetbrains.bio.genome

import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.size
import java.net.URI
import java.nio.file.Path
import java.util.stream.Collectors
import java.util.stream.Stream

/**
 * @author Oleg Shpynov
 * @since 11/04/2018.
 */

object PeaksInfo {
    val CT_COUNT = TrackAboutLongColumnType("Count") // peaks or regions count
    val CT_TOTAL_LEN = TrackAboutLongColumnType("Total length") // peaks or regions aggregated length
    val CT_GENOME_COVERAGE = TrackAboutPercentageColumnType("Genome coverage")
    val CT_MIN_LEN = TrackAboutLongColumnType("Min length")
    val CT_MAX_LEN = TrackAboutLongColumnType("Max length")
    val CT_MEAN_LEN = TrackAboutLongColumnType("Mean length")
    val CT_MEDIAN_LEN = TrackAboutLongColumnType("Median length")
    val CT_FRIP = TrackAboutDoubleColumnType("FRIP")

    fun compute(
        genomeQuery: GenomeQuery,
        peaksStream: Stream<Location>,
        src: URI?,
        paths: List<Path>,
        fragment: Fragment = AutoFragment
    ): List<TrackAboutMetricValue<*>> {
        val peaks = peaksStream.collect(Collectors.toList())
        val peaksLengths = peaks.map { it.length().toDouble() }.toDoubleArray()
        val peaksCount = peaksLengths.count()
        val peaksLenSum = peaksLengths.sum()
        val coverage = peaksLenSum / genomeQuery.get().sumOf { it.length.toLong() }

        val result = arrayListOf<TrackAboutMetricValue<*>>()
        if (src != null) {
            result.add(TrackAboutColumnTypes.CT_FILE to src)
            result.add(TrackAboutColumnTypes.CT_FILE_SIZE to src.size)
        }
        result.add(CT_COUNT to peaksCount)
        result.add(CT_TOTAL_LEN to peaksLenSum.toLong())
        result.add(CT_GENOME_COVERAGE to coverage)
        result.add(CT_MIN_LEN to if (peaksLengths.isEmpty()) 0L else StatUtils.min(peaksLengths).toLong())
        result.add(CT_MAX_LEN to if (peaksLengths.isEmpty()) 0L else StatUtils.max(peaksLengths).toLong())
        result.add(CT_MEAN_LEN to if (peaksLengths.isEmpty()) 0L else peaksLengths.average().toLong())
        result.add(
            CT_MEDIAN_LEN to if (peaksLengths.isEmpty()) 0L else StatUtils.percentile(peaksLengths, 50.0).toLong()
        )

        // Don't recompute tags coverage if it is not processed locally
        if (paths.isNotEmpty()) {
            val readQueries = paths.map { ReadsQuery(genomeQuery, it, null, true, fragment, showLibraryInfo = false) }
            if (readQueries.all { it.npzPath().isAccessible() }) {
                val coverages = readQueries.map { it.get() }
                val frip = frip(genomeQuery, peaks, coverages)
                result.add(CT_FRIP to frip)
            }
        }
        return result
    }

    private fun frip(genomeQuery: GenomeQuery, peakLocations: List<Location>, coverages: List<Coverage>): Double {
        return coverages.map { coverage ->
            1.0 * peakLocations.sumOf { coverage.getBothStrandsCoverage(it.toChromosomeRange()).toLong() } /
                    genomeQuery.get().sumOf {
                        coverage.getBothStrandsCoverage(ChromosomeRange(0, it.length, it)).toLong()
                    }
        }.average()
    }
}
