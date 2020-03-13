package org.jetbrains.bio.genome

import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.presentablePath
import org.jetbrains.bio.util.size
import org.slf4j.LoggerFactory
import java.net.URI
import java.nio.file.Path
import java.util.*
import java.util.stream.Collectors
import java.util.stream.Stream

/**
 * @author Oleg Shpynov
 * @since 11/04/2018.
 */

object PeaksInfo {

    private val LOG = LoggerFactory.getLogger(PeaksInfo::class.java)

    private fun Long.formatLongNumber() = String.format("%,d", this).replace(',', ' ')

    fun compute(
            genomeQuery: GenomeQuery,
            peaksStream: Stream<Location>,
            src: URI?,
            paths: List<Path>,
            fragment: Fragment = AutoFragment
    ): Map<String, String> {
        val peaks = peaksStream.collect(Collectors.toList())
        val peaksLengths = peaks.map { it.length().toDouble() }.toDoubleArray()
        val peaksCount = peaksLengths.count()
        val peaksLenSum = peaksLengths.sum()
        val coverage = peaksLenSum / genomeQuery.get().map { it.length.toLong() }.sum()
        val result = LinkedHashMap<String, String>()
        if (src != null) {
            result["Track source"] = src.presentablePath()
            result["Source size"] = "${src.size}${if (!src.isAccessible()) " <not accessible>" else ""}"
        }
        result["Peak count"] = peaksCount.toLong().formatLongNumber()
        result["Total peak length"] = peaksLenSum.toLong().formatLongNumber()
        result["Genome coverage"] = String.format("%.2f%%", coverage)
        result["Min length"] = (if (peaksLengths.isEmpty()) 0L else StatUtils.min(peaksLengths).toLong()).formatLongNumber()
        result["Max length"] = (if (peaksLengths.isEmpty()) 0L else StatUtils.max(peaksLengths).toLong()).formatLongNumber()
        result["Mean length"] = (if (peaksLengths.isEmpty()) 0L else peaksLengths.average().toLong()).formatLongNumber()
        result["Median length"] =
                (if (peaksLengths.isEmpty()) 0L else StatUtils.percentile(peaksLengths, 50.0).toLong()).formatLongNumber()

        // Don't recompute tags coverage if it is not processed locally
        if (paths.isNotEmpty()) {
            val readQueries = paths.map { ReadsQuery(genomeQuery, it, true, fragment) }
            if (readQueries.all { it.npzPath().isAccessible() }) {
                val coverages = readQueries.map { it.get() }
                val frip = frip(genomeQuery, peaks, coverages)
                result["FRIP"] = frip.toString()
            }
        }
        return result
    }

    private fun frip(genomeQuery: GenomeQuery, peakLocations: List<Location>, coverages: List<Coverage>): Double {
        val frip = coverages.map { coverage ->
            1.0 * peakLocations.map { coverage.getBothStrandsCoverage(it.toChromosomeRange()).toLong() }.sum() /
                    genomeQuery.get().map {
                        coverage.getBothStrandsCoverage(ChromosomeRange(0, it.length, it)).toLong()
                    }.sum()
        }.average()
        LOG.debug("Frip: $frip")
        return frip
    }
}
