package org.jetbrains.bio.experiments.gsea

import joptsimple.ValueConversionException
import joptsimple.ValueConverter
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
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
        val progress = Progress {title="Loci Sampling"}.bounded(simulationsNumber.toLong())

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

            // XXX: shuffled regions not intersect by def of our shuffle procedure
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
        regionLabelAndLociToTest: List<Pair<String, LocationsList<out RangesList>>>,
        outputFolderPath: Path? = null,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER,
        aSetIsLoi: Boolean = true,
        mergeOverlapped: Boolean = true,
        intersectionFilter: LocationsList<out RangesList>? = null
    ): DataFrame {
        outputFolderPath?.createDirectories()

        val label2Stats = regionLabelAndLociToTest.map { it.first to TestedRegionStats() }.toMap()

        // TODO: keep type of list i.e mergeOverlapped
        val sourceLoci = readLocations(srcRegionsPath, genome, mergeOverlapped, intersectionFilter)

        val nChunks = ceil(simulationsNumber.toDouble() / chunkSize).toInt()

        val progress = Progress {title="Over/Under-representation check progress"}.bounded(
            nChunks.toLong() * regionLabelAndLociToTest.size.toLong()
        )
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info("Simulations: Chunk [${chunkId+1} of $nChunks], simulations ${start.formatLongNumber()}..${end.formatLongNumber() } of ${simulationsNumber.formatLongNumber()}")
            val sampledRegions: List<LocationsList<out RangesList>> =
                sampleRegions(srcRegionsPath, simulationsInChunk).let { result ->
                    if (intersectionFilter == null) {
                        result
                    } else {
                        result.map { ll ->
                            val filtered = ll.intersectRanges(intersectionFilter)
                            LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
                        }
                    }
                }

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

                val metricValueForSrc = calcMetric(sourceLoci, lociToTest, aSetIsLoi, metric)
                val metricValueForSampled = sampledRegions.stream().parallel().mapToLong { sampledLoci ->
                    calcMetric(sampledLoci, lociToTest, aSetIsLoi, metric)
                }.toArray()

                if (outputFolderPath != null) {
                    DataFrame().with("sample", metricValueForSampled)
                        .save(outputFolderPath / "${regionLabel}_${metric.column}.chunk${chunkId + 1}.tsv")
                }

                // calc 2 sided: p-value using definition:
                var countGreater = 0 // count when random regions metric value is above given value
                var countBelow = 0 // count when random regions metric value is above given value
                for (v in metricValueForSampled) {
                    if (v >= metricValueForSrc) {
                        countGreater++
                    }
                    if (v <= metricValueForSrc) {
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

    private fun calcMetric(
        sampledLoci: LocationsList<out RangesList>,
        lociToTest: LocationsList<out RangesList>,
        aSetIsLoi: Boolean,
        metric: RegionsMetric
    ): Long {
        val a = if (aSetIsLoi) sampledLoci else lociToTest
        val b = if (aSetIsLoi) lociToTest else sampledLoci

        // XXX: at the moment only 'overlap' is used here, i.e. integer metric
        return metric.calcMetric(a, b).toLong()
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStats::class.java)

        fun readLocations(
            path: Path, genome: Genome, mergeOverlapped: Boolean,
            intersectionFilter: LocationsList<out RangesList>? = null
        ): LocationsList<out RangesList> = genome.toQuery().let { gq ->
            val locations = readLocations(path, BedFormat.auto(path), gq)
            if (mergeOverlapped) {
                val mergedLocations = LocationsMergingList.create(gq, locations)
                if (locations.size != mergedLocations.size) {
                    LOG.info("$path: ${locations.size} regions merged to ${mergedLocations.size}")
                }
                if (intersectionFilter == null) {
                    mergedLocations
                } else {
                    val filtered = mergedLocations.intersectRanges(intersectionFilter)
                    LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
                }
            } else {
                val result = LocationsSortedList.create(gq, locations)
                if (intersectionFilter == null) {
                    result
                } else {
                    val filtered = result.intersectRanges(intersectionFilter)
                    LocationsSortedList.create(filtered.genomeQuery, filtered.locationIterator())
                }
            }
        }

        fun readLocations(path: Path, bedFormat: BedFormat, gq: GenomeQuery): List<Location> {
            var recordsNumber = 0
            val ignoredChrs = mutableListOf<String>()

            val loci = bedFormat.parse(path) { bedParser ->
                bedParser.mapNotNull {
                    recordsNumber++

                    val chr = gq[it.chrom]
                    if (chr == null) {
                        ignoredChrs.add(it.chrom)
                        // skip all unmapped contigs, etc
                        null
                    } else {
                        Location(it.start, it.end, chr)
                    }
                }
            }
            if (loci.size != recordsNumber) {
                val pnt = loci.size.asPercentOf(recordsNumber)
                LOG.warn("$path: Loaded $pnt % (${loci.size} of $recordsNumber) locations. Ignored chromosomes: ${ignoredChrs.size}. For more details use debug option.")
            }
            if (ignoredChrs.isNotEmpty()) {
                LOG.debug("$path: Ignored chromosomes: $ignoredChrs")
            }
            return loci
        }
    }
}

fun Int.asPercentOf(total: Int, digitsAfterDot: Int = 2) = Precision.round(100.0 * this / total, digitsAfterDot)

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