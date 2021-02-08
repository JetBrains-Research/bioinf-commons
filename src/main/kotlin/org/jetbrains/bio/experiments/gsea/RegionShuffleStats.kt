package org.jetbrains.bio.experiments.gsea

import gnu.trove.map.hash.TIntIntHashMap
import joptsimple.ValueConversionException
import joptsimple.ValueConverter
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
import org.jetbrains.bio.util.parallelismLevel
import org.jetbrains.bio.viktor.KahanSum
import org.jetbrains.bio.viktor.asF64Array
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.stream.Collectors
import java.util.stream.IntStream
import kotlin.math.ceil
import kotlin.math.min
import kotlin.math.pow
import kotlin.math.sqrt

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
        val progress = Progress { title = "Loci Sampling" }.bounded(simulationsNumber.toLong())

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
        val dumpDetails = outputFolderPath != null

        LOG.info("Regions sets to test: ${regionLabelAndLociToTest.size}")

        val label2Stats = regionLabelAndLociToTest.map { it.first to TestedRegionStats() }.toMap()
        val sourceLoci = readLocations(srcRegionsPath, genome, mergeOverlapped, intersectionFilter)
        val nChunks = ceil(simulationsNumber.toDouble() / chunkSize).toInt()

        val progress = Progress { title = "Over/Under-representation check progress (all chunks)" }.bounded(
            nChunks.toLong() * regionLabelAndLociToTest.size.toLong()
        )
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info("Simulations: Chunk [${chunkId + 1} of $nChunks], simulations " +
                    "${start.formatLongNumber()}..${end.formatLongNumber()} of ${simulationsNumber.formatLongNumber()}")

            val sampledRegions: List<List<LocationsList<out RangesList>>> =
                sampleRegions(simulationsInChunk, intersectionFilter, parallelismLevel())

            // TODO: maybe roll back, issue was other.
//            val sampledRegions: List<LocationsList<out RangesList>> =
//                sampleRegions(srcRegionsPath, simulationsInChunk).let { result ->
//                    if (intersectionFilter == null) {
//                        result
//                    } else {
//                        result.map { ll ->
//                            val filtered = ll.intersectRanges(intersectionFilter)
//                            LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
//                        }
//                    }
//                }

            /*
            // XXX: optional save sampledRegions (flatmap)
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

                val metricValueForSampled =
                    calcMetricForSampled(sampledRegions, lociToTest, aSetIsLoi, metric, metricValueForSrc)

//                val metricValueForSampled = chunked(sampledRegions.stream(), parallelismLevel()).parallel().flatMapToLong { chunk ->
//                                    chunk.stream().mapToLong {
//                                        calcMetric(it, lociToTest, aSetIsLoi, metric)
//                                    }
//                                }.toArray()

                val stats = label2Stats[regionLabel]!!
                metricValueForSampled.forEach { st ->
                    stats.countSetsWithMetricsAboveThr += st.countSetsWithMetricsAboveThr
                    stats.countSetsWithMetricsBelowThr += st.countSetsWithMetricsBelowThr
                    stats.metricHist += st.metricHist
                }
                stats.simulationsNumber = stats.simulationsNumber + simulationsInChunk
                stats.metricValueForSrc = metricValueForSrc
            }
        }

        progress.done()

        // Save results to DataFrame:
        val regionLabels = regionLabelAndLociToTest.map { it.first }
        val testLociNumber = regionLabelAndLociToTest.map { it.second.size }.toIntArray()
        val pValuesList = ArrayList<Double>()
        val srcMetricValues = ArrayList<Long>()
        val sampledSetsMetricMedian = ArrayList<Int>()
        val sampledSetsMetricVar = ArrayList<Double>()
        val sampledSetsMetricMean = ArrayList<Double>()

        regionLabels.forEach { regionLabel ->
            val stats = label2Stats[regionLabel]!!
            pValuesList.add(stats.pvalue(hypAlt))
            srcMetricValues.add(stats.metricValueForSrc)

            // For large simulation & multiple regions number we cannot store in memory
            // all metrics values - to many RAM required (e.g 7k x 10^6 simulations > 90 GB)
            // instead we could calc them from hist or use 'online' version of std and mean and skip
            // median and total hist.
            val metricHist = stats.metricHist

            if (dumpDetails) {
                val path = outputFolderPath!! / "${regionLabel}_${metric.column}.hist.tsv"
                metricHist.save(path)
            }
            require(simulationsNumber == metricHist.countValues()) {
                "Simulations number doesn't match expectations:"
            }
            val metricMean = metricHist.mean()
            val metricSd = metricHist.stdev(simulationsNumber, metricMean)
            val medianValue = metricHist.median(simulationsNumber)

            sampledSetsMetricMedian.add(medianValue)
            sampledSetsMetricMean.add(metricMean)
            sampledSetsMetricVar.add(metricSd * metricSd)
        }
        val pValues = pValuesList.toDoubleArray()
        val qValues = Multiple.adjust(pValues.asF64Array()).toDoubleArray()

        return DataFrame()
            .with("name", regionLabels.toTypedArray())
            .with("test_loci_number", testLociNumber)
            .with(metric.column, srcMetricValues.toLongArray())
            .with("sampled_median_${metric.column}", sampledSetsMetricMedian.toIntArray())
            .with("sampled_mean_${metric.column}", sampledSetsMetricMean.toDoubleArray())
            .with("sampled_var_${metric.column}", sampledSetsMetricVar.toDoubleArray())
            .with("src_loci_number", IntArray(regionLabels.size) { sourceLoci.size })
            .with("sampled_sets_number", IntArray(regionLabels.size) { simulationsNumber })
            .with("pValue", pValues)
            .with("qValue", qValues)
            .reorder("qValue")
    }

    private fun calcMetricForSampled(
        sampledRegions: List<List<LocationsList<out RangesList>>>,
        lociToTest: LocationsList<out RangesList>,
        aSetIsLoi: Boolean,
        metric: RegionsMetric,
        metricValueForSrc: Long
    ): List<PerThreadStats> = sampledRegions.parallelStream().map { chunk ->
        var countGreater = 0 // count when random regions metric value is above given value
        var countBelow = 0 // count when random regions metric value is above given value

        val metricHist = IntHistogram()
        chunk.forEachIndexed { i, ll ->
            val v = calcMetric(ll, lociToTest, aSetIsLoi, metric)
            require(v < Int.MAX_VALUE) {
                "Long values not supported, please fire a ticket. Got: $v >= ${Int.MAX_VALUE}"
            }
            metricHist.increment(v.toInt())

            // calc 2 sided: p-value using definition:
            if (v >= metricValueForSrc) {
                countGreater++
            }
            if (v <= metricValueForSrc) {
                countBelow++
            }
        }
        PerThreadStats(countGreater, countBelow, metricHist)
    }.collect(Collectors.toList())

    private fun sampleRegions(
        simulationsNumber: Int,
        intersectionFilter: LocationsList<out RangesList>?,
        threadsNum: Int
    ): List<List<LocationsMergingList>> {
        val chunked = sampleRegions(srcRegionsPath, simulationsNumber).chunked(threadsNum)
        if (intersectionFilter == null) {
            return chunked
        }

        return chunked.map { chunk ->
            chunk.map { ll ->
                val filtered = ll.intersectRanges(intersectionFilter)
                LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
            }
        }
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

data class PerThreadStats(
    val countSetsWithMetricsAboveThr: Int,
    val countSetsWithMetricsBelowThr: Int,
    val metricHist: IntHistogram = IntHistogram()
)

data class TestedRegionStats(
    var countSetsWithMetricsAboveThr: Int = 0,
    var countSetsWithMetricsBelowThr: Int = 0,
    var simulationsNumber: Int = 0,
    var metricValueForSrc: Long = 0L,
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

class IntHistogram {
    val data: TIntIntHashMap = TIntIntHashMap()

    fun increment(value: Int) {
        data.adjustOrPutValue(value, 1, 1)
    }

    fun mean(valuesNumber: Int = countValues()): Double {
        // sum(x_i) / n
        val metricMeanAcc = KahanSum();
        data.forEachEntry { metric, count ->
            metricMeanAcc += (count * metric).toDouble()
            true
        }
        return metricMeanAcc.result() / valuesNumber
    }

    fun countValues(): Int {
        var totalCount = 0
        data.forEachEntry { _, count ->
            totalCount += count
            true
        }
        return totalCount
    }

    operator fun plusAssign(other: IntHistogram) {
        other.data.forEachEntry { metric, count ->
            data.adjustOrPutValue(metric, count, count)
            true
        }
    }

    /**
     * Returns '-1' if list is empty
     */
    fun median(valuesNumber: Int = countValues()): Int {
        val medianIdx = valuesNumber / 2
        var idx = 0

        val keys = data.keys()
        keys.sort()

        for (metric in keys) {
            val count = data[metric]
            idx += count
            if (idx > medianIdx) {
                return metric
            }
        }
        return -1
    }

    fun stdev(valuesNumber: Int = countValues(), mean: Double = mean(valuesNumber)): Double {
        //listOf(1,3).stream().mapToInt { it }.average()
        //XXX: https://math.stackexchange.com/questions/857566/how-to-get-the-standard-deviation-of-a-given-histogram-image
        //XXX: online variance - https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance  and https://math.stackexchange.com/questions/198336/how-to-calculate-standard-deviation-with-streaming-inputs

        if (valuesNumber == 0) {
            return Double.NaN
        } else if (valuesNumber == 1) {
            return 0.0
        }
        //  sum((x_i - mean)^2) / (n - 1)
        val metricVarAcc = KahanSum();
        data.forEachEntry { metric, count ->
            metricVarAcc += count * (metric - mean).pow(2)
            true
        }
        return sqrt(metricVarAcc.result() / (valuesNumber - 1))
    }

    fun save(path: Path) {
        val keys = data.keys()
        keys.sort()
        DataFrame()
            .with("metric", keys)
            .with("count", keys.map { data[it] }.toIntArray())
            .save(path)
    }

    override fun toString() = data.toString()
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IntHistogram) return false

        if (data != other.data) return false

        return true
    }

    override fun hashCode(): Int = data.hashCode()

    companion object {
        fun create(data: IntArray): IntHistogram {
            val metricHist = IntHistogram()
            data.forEach { metricHist.increment(it) }
            return metricHist
        }
    }
}