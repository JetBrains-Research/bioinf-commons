package org.jetbrains.bio.genome.coverage

import com.google.common.collect.MinMaxPriorityQueue
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.max
import kotlin.math.sqrt

object Roughness {

    private val LOG = LoggerFactory.getLogger(Roughness::class.java)

    fun coverage(chromosome: Chromosome,
                 start: Int, end: Int,
                 treatmentCoverage: Coverage,
                 treatmentScale: Double,
                 controlCoverage: Coverage?,
                 controlScale: Double?,
                 beta: Double): Double {
        val chromosomeRange = ChromosomeRange(start, end, chromosome)
        val tc = treatmentCoverage.getBothStrandsCoverage(chromosomeRange) * treatmentScale
        return if (controlCoverage != null && controlScale != null) {
            val cc = controlCoverage.getBothStrandsCoverage(chromosomeRange) * controlScale
            max(0.0, tc - beta * cc)
        } else {
            tc
        }
    }

    val REGION_LEN = 10_000
    val TOP_REGIONS = 10_000
    val WORK_REGIONS = 200
    val RESOLUTION = 100

    /**
     * Detects roughness of the coverage.
     */
    fun detectAverageRoughness(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        treatmentScale: Double,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        blackListPath: Path? = null,
        region_len: Int = REGION_LEN,
        top_regions: Int = TOP_REGIONS,
        work_regions: Int = WORK_REGIONS,
        resolution: Int = RESOLUTION
    ): Double {
        LOG.debug("Compute coverage in regions")
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val genome_regions = limitedQuery.get().sumOf { floor(it.length.toDouble() / region_len).toLong() }
        check(top_regions < genome_regions) {
            "Too many top regions $top_regions > $genome_regions"
        }

        val comparator = Comparator<Triple<Chromosome, Int, Double>> { o1, o2 -> o2.third.compareTo(o1.third) }
        val region_coverages: MinMaxPriorityQueue<Triple<Chromosome, Int, Double>> =
            MinMaxPriorityQueue
                .orderedBy(comparator)
                .maximumSize(top_regions)
                .create()
        var regions = 0
        val blackList = if (blackListPath != null) LocationsMergingList.load(limitedQuery, blackListPath) else null
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            for (i in 0 until floor(chr.length.toDouble() / region_len).toInt()) {
                regions += 1
                val start = region_len * i
                val end = region_len * (i + 1)
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                val c = coverage(
                    chr, start, end,
                    treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
                )
                region_coverages.add(Triple(chr, i, c))
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, regions)
        }

        val region_coverages_array = region_coverages.toTypedArray()
        region_coverages_array.sortWith(comparator)

        val step = if (region_coverages_array.size > work_regions) {
            LOG.debug("Pick $work_regions / $top_regions uniform regions for computation speedup")
            ceil(region_coverages_array.size.toDouble() / work_regions).toInt()
        } else
            1

        val std_means = DoubleArray(work_regions)
        for (i in region_coverages_array.indices) {
            if (i % step != 0) {
                continue
            }
            val (chr, start, _) = region_coverages_array[i]
            val stats = DoubleArray(region_len / resolution) {
                coverage(chr, start + it * resolution, start + (it + 1) * resolution,
                    treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
                )
            }
            val mean = stats.average()
            val std = stats.standardDeviation()
            std_means[i / step] = if (mean > 0) std / mean else 0.0
        }
        return std_means.average()
    }

}

fun DoubleArray.standardDeviation(): Double {
    var sum = 0.0
    var sumSq = 0.0
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return sqrt((sumSq - sum * sum / size) / size)
}
