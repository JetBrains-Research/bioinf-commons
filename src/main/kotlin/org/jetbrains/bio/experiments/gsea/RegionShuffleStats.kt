package org.jetbrains.bio.experiments.gsea

import joptsimple.ValueConversionException
import joptsimple.ValueConverter
import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.intersection.IntersectionMetric
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.sampling.shuffleChromosomeRanges
import org.jetbrains.bio.statistics.hypothesis.Multiple
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.div
import org.jetbrains.bio.viktor.asF64Array
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.stream.Collectors
import java.util.stream.IntStream
import kotlin.math.ceil
import kotlin.math.min

/**
 * TODO: keeps all sampling results in memory, although could work and stream and consume less memory
 */

/**
 * @param genome Genome
 * @param srcRegionsPath Source regions, their length will be used for simulations
 * @param backgroundRegions Simulated regions will be taken from background
 * @param simulationsNumber Simulations number
 * @param maxRetries Max retries if cannot shuffle region from background
 */
class RegionShuffleStats(
    private val genome: Genome,
    private val srcRegionsPath: Path,
    private val backgroundRegions: Path?,
    private val simulationsNumber: Int,
    private val chunkSize: Int,
    private val maxRetries: Int
) {

    private fun sampleRegions(lociFilePath: Path, simulationsNumber: Int): List<LocationsMergingList> {
        val progress = Progress {title="Sampling"}.bounded(simulationsNumber.toLong())

        val loci = readRanges(lociFilePath)
        val background = backgroundRegions?.let { readRanges(it) }
        val sampled = IntStream.range(0, simulationsNumber).parallel().mapToObj { _ ->
            val randLoci = shuffleChromosomeRanges(
                genome.toQuery(),
                loci,
                background,
                maxRetries = maxRetries
            ).map { it.on(Strand.PLUS) }

            progress.report()
            LocationsMergingList.create(genome.toQuery(), randLoci)

        }.collect(Collectors.toList())
        progress.done()

        return sampled
    }

    fun readRanges(dmrFile: Path) =
            readLocations(dmrFile, BedFormat.auto(dmrFile), genome.toQuery())
                    .map { it.toChromosomeRange() }.toList()


    /**
     * @param regionLabelAndLociToTest  Test sampled loci vs given regions lists (label and loci) using given metric
     * @param outputFolderPath Results path
     * @param metric Metric will be applied to  `(sampledLoci, lociToTest)`
     */
    fun calcStatistics(
        regionLabelAndLociToTest: List<Pair<String, LocationsMergingList>>,
        outputFolderPath: Path? = null,
        metric: IntersectionMetric = IntersectionMetric.OVERLAP,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER
    ): DataFrame {
        outputFolderPath?.createDirectories()

        val label2Stats = regionLabelAndLociToTest.map { it.first to TestedRegionStats() }.toMap()
        val sourceLoci = readLocations(srcRegionsPath, genome)
        val nChunks = ceil(simulationsNumber.toDouble() / chunkSize).toInt()

        val progress = Progress {title="Simulation progress"}.bounded(nChunks.toLong() * regionLabelAndLociToTest.size.toLong())
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info("Simulations: Chunk [$chunkId of $nChunks], simulations ${start.formatLongNumber()}..${end.formatLongNumber() } of ${simulationsNumber.formatLongNumber()}")
            val sampledRegions = sampleRegions(srcRegionsPath, simulationsInChunk)

            /*
            // XXX: optional save sampledRegions
            val dumpSampledLoci = false
            if (dumpSampledLoci && outputFolderPath != null) {
                val sampledDir = outputFolderPath / "sampled"
                sampledDir.toFile().mkdirs()
                sampledRegions.forEachIndexed { idx, sampledSet ->
                    (sampledDir / "sampled_$idx.tsv").bufferedWriter().use { w ->
                        sampledSet.asLocationSequence().forEach {
                            w.write("${it.chromosome.name}\t${it.startOffset}\t${it.endOffset}")
                        }
                    }
                }
            }
            */

            for ((regionLabel, lociToTest) in regionLabelAndLociToTest) {
                progress.report()

                val metricValueForSrc = metric.calcMetric(sourceLoci, lociToTest)

                val metricValueForSampled = sampledRegions.stream().parallel().mapToLong { sampledLoci ->
                    metric.calcMetric(sampledLoci, lociToTest)
                }.toArray()

                if (outputFolderPath != null) {
                    DataFrame().with("sample", metricValueForSampled)
                        .save(outputFolderPath / "${regionLabel}_${metric.column}.chunk${chunkId+1}.tsv")
                }

                // calc 2 sided: p-value using definition:
                var countGreater = 0 // count when random regions metric value is above given value
                var countBelow = 0 // count when random regions metric value is above given value
                for (v in metricValueForSampled) {
                    if (v >= metricValueForSrc) {
                        countGreater++
                    }
                    if (v <= countBelow) {
                        countBelow++
                    }
                }

                val stats = label2Stats[regionLabel]!!
                stats.countSetsWithMetricsAboveThr = stats.countSetsWithMetricsAboveThr + countGreater
                stats.countSetsWithMetricsBelowThr = stats.countSetsWithMetricsBelowThr + countBelow
                stats.simulationsNumber = stats.simulationsNumber + simulationsInChunk
                stats.metricValueForSrc = metricValueForSrc
                stats.metricValuesForSampledChunkSets.add(metricValueForSampled)
            }
        }

        progress.done()

        // Save results to DataFrame:
        val regionLabels = regionLabelAndLociToTest.map { it.first }
        val testLociNumber = regionLabelAndLociToTest.map { it.second.size }.toIntArray()
        val sourceLociNumber = IntArray(regionLabels.size) { sourceLoci.size}
        val pValuesList = ArrayList<Double>()
        val srcMetricValues = ArrayList<Long>()
        val sampledSetsMetricMedian = ArrayList<Long>()
        val sampledSetsMetricVar = ArrayList<Double>()
        val sampledSetsMetricMean = ArrayList<Double>()
        regionLabels.forEach { regionLabel ->
            val stats = label2Stats[regionLabel]!!
            pValuesList.add(stats.pvalue(hypAlt))
            srcMetricValues.add(stats.metricValueForSrc)

            val values = stats.metricValuesForSampledChunkSets.flatMap { chunk ->
                chunk.map { it.toDouble() }
            }.toDoubleArray()
            values.sort()

            sampledSetsMetricMedian.add(values[values.size / 2].toLong())
            sampledSetsMetricMean.add(StatUtils.mean(values))
            sampledSetsMetricVar.add(StatUtils.variance(values))
        }
        val pValues = pValuesList.toDoubleArray()
        val qValues = Multiple.adjust(pValues.asF64Array()).toDoubleArray()

        return DataFrame()
                .with("name", regionLabels.toTypedArray())
                .with("test_loci_number", testLociNumber)
                .with(metric.column, srcMetricValues.toLongArray())
                .with("sampled_median_${metric.column}", sampledSetsMetricMedian.toLongArray())
                .with("sampled_mean_${metric.column}", sampledSetsMetricMean.toDoubleArray())
                .with("sampled_var_${metric.column}", sampledSetsMetricVar.toDoubleArray())
                .with("src_loci_number", sourceLociNumber)
                .with("pValue", pValues)
                .with("qValue", qValues)
                .reorder("qValue")
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStats::class.java)

        fun readLocations(path: Path, genome: Genome) = genome.toQuery().let { gq ->
            LocationsMergingList.create(
                gq,
                // parse location directly into LocationsMergingList
                readLocations(path, BedFormat.auto(path), gq)
            )

        }

        fun readLocations(path: Path, bedFormat: BedFormat, gq: GenomeQuery) =
                bedFormat.parse(path) { bedParser ->
                    bedParser.mapNotNull {
                        val chr = gq[it.chrom]
                        if (chr == null) {
                            // skip all unmapped contigs, etc
                            null
                        } else {
                            Location(it.start, it.end, chr)
                        }
                    }
                }

    }
}

data class TestedRegionStats(
    var countSetsWithMetricsAboveThr: Int = 0,
    var countSetsWithMetricsBelowThr: Int = 0,
    var simulationsNumber: Int = 0,
    var metricValueForSrc: Long = 0L,
    val metricValuesForSampledChunkSets: MutableList<LongArray> =  mutableListOf()
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


private fun Int.formatLongNumber() = String.format("%,d", this).replace(',', ' ')

enum class PermutationAltHypothesis(private val presentableString: String) {
    TWO_SIDED("two-sided"),
    GREATER("greater"),
    LESS("less");

    override fun toString() = presentableString

    companion object {
        fun converter() = object : ValueConverter<PermutationAltHypothesis> {
            override fun convert(value: String): PermutationAltHypothesis {
                for (h in values()) {
                    if (h.presentableString == value) {
                       return h
                    }
                }
                throw ValueConversionException("Unsupported hypothesis: $value")
            }

            override fun valueType() = PermutationAltHypothesis::class.java

            override fun valuePattern() = null
        }
    }
}