package org.jetbrains.bio.gsea

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.statistics.hypothesis.BenjaminiHochberg
import org.jetbrains.bio.util.*
import org.jetbrains.bio.viktor.asF64Array
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.stream.Collectors
import java.util.stream.IntStream
import kotlin.math.ceil

/**
 * @param simulationsNumber Simulations number
 * @param chunkSize Size of chunk
 * @param maxRetries Max retries when cannot shuffle region from background
 *
 * Loci file path will be added to background
 */
abstract class RegionShuffleStats(
    protected  val simulationsNumber: Int,
    protected  val chunkSize: Int,
    protected  val maxRetries: Int
) {
    /**
     * Strand is ignore in all files
     *
     * @param loiInfos  Test sampled loci vs given regions lists (label, loci, initial #regions in loci) using given metric
     * @param outputFolderPath Results path
     * @param metric Metric will be applied to  `(sampledLoci, lociToTest)`
     */
    abstract fun calcStatistics(
        gq: GenomeQuery,
        inputRegionsPath: Path,
        backgroundPath: Path?,
        loiInfos: List<LoiInfo>,
        outputFolderPath: Path? = null,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER,
        aSetIsRegions: Boolean = true,
        mergeOverlapped: Boolean = true,
        intersectionFilter: LocationsList<out RangesList>? = null,
        genomeMaskedAreaPath: Path? = null,
        genomeAllowedAreaPath: Path? = null,
        mergeRegionsToBg: Boolean = false,
        samplingWithReplacement: Boolean = false
    ): DataFrame


    protected fun <T> doCalcStatistics(
        gq: GenomeQuery,
        loiInfos: List<LoiInfo>,
        outputFolderPath: Path? = null,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER,
        aSetIsRegions: Boolean = true,
        mergeOverlapped: Boolean = true,
        intersectionFilter: LocationsList<out RangesList>? = null,
        samplingWithReplacement: Boolean = false,
        inputRegionsAndBackgroundProvider: (GenomeQuery) -> Pair<List<Location>, T>,
        samplingFun: (GenomeQuery, List<ChromosomeRange>, T, Int, Boolean) -> List<Location>
    ): DataFrame {
        outputFolderPath?.createDirectories()
        val dumpDetails = outputFolderPath != null

        val (inputRegions, bgRegionsList) = inputRegionsAndBackgroundProvider(gq)
        require(inputRegions.isNotEmpty()) {
            "Regions file is empty or all regions were removed by filters."
        }
        val allowedInputRegionsList = when {
            mergeOverlapped -> LocationsMergingList.create(gq, inputRegions)
            else -> LocationsSortedList.create(gq, inputRegions)
        }

        LOG.info("LOI sets to test: ${loiInfos.size}")

        val label2Stats = loiInfos.associate { info ->
            info.label to TestedRegionStats()
        }

        val label2LoadedRegionsNumber = loiInfos.associate { info ->
            info.label to info.processedLoiNumber
        }
        val label2RecordsNumber = loiInfos.associate { info ->
            info.label to info.recordsNumber
        }

        val nChunks = ceil(simulationsNumber.toDouble() / chunkSize).toInt()

        val progress = Progress { title = "Over/Under-representation check progress (all chunks)" }.bounded(
            nChunks.toLong()
        )
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info(
                "Simulations: Chunk [${chunkId + 1} of $nChunks], simulations " +
                        "${start.formatLongNumber()}..${end.formatLongNumber()} of ${simulationsNumber.formatLongNumber()}"
            )

            val sampledRegions: List<List<LocationsList<out RangesList>>> =
                sampleRegions(
                    simulationsInChunk, intersectionFilter, inputRegions, bgRegionsList,
                    gq,
                    parallelismLevel(),
                    samplingFun=samplingFun,
                    withReplacement = samplingWithReplacement
                )


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
            val chunkProgress = Progress { title = "Chunk Over/Under-representation:" }.bounded(
                loiInfos.size.toLong()
            )
            for (info in loiInfos) {
                chunkProgress.report()

                val metricValueForInput = calcMetric(allowedInputRegionsList, info.loci, aSetIsRegions, metric)

                val metricValueForSampled =
                    calcMetricForSampled(sampledRegions, info.loci, aSetIsRegions, metric, metricValueForInput)

                val stats = label2Stats[info.label]!!
                metricValueForSampled.forEach { st ->
                    stats.countSetsWithMetricsAboveThr += st.countSetsWithMetricsAboveThr
                    stats.countSetsWithMetricsBelowThr += st.countSetsWithMetricsBelowThr
                    stats.metricHist += st.metricHist
                }
                stats.simulationsNumber += simulationsInChunk
                stats.metricValueForInput = metricValueForInput
            }
            chunkProgress.done()
            progress.report()
        }

        progress.done()

        // Save results to DataFrame:
        val loiLabels = loiInfos.map { it.label }
        val loiRangesNumber = loiInfos.map { it.loci.size }.toIntArray()
        val pValuesList = ArrayList<Double>()
        val metricValuesForInput = ArrayList<Long>()
        val sampledSetsMetricMedian = ArrayList<Int>()
        val sampledSetsMetricVar = ArrayList<Double>()
        val sampledSetsMetricMean = ArrayList<Double>()

        loiLabels.forEach { loiLabel ->
            val stats = label2Stats[loiLabel]!!
            pValuesList.add(stats.pvalue(hypAlt))
            metricValuesForInput.add(stats.metricValueForInput)

            // For large simulation & multiple regions number we cannot store in memory
            // all metrics values - to many RAM required (e.g 7k x 10^6 simulations > 90 GB)
            // instead we could calc them from hist or use 'online' version of std and mean and skip
            // median and total hist.
            val metricHist = stats.metricHist

            if (dumpDetails) {
                val path = outputFolderPath!! / "${loiLabel}_${metric.column}.hist.tsv"
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
        val qValues = BenjaminiHochberg.adjust(pValues.asF64Array()).toDoubleArray()

        val metricArgString = if (aSetIsRegions) "regions, loi" else "loi, regions"
        return DataFrame()
            .with("loi", loiLabels.toTypedArray())
            .with("n_loi_records", loiLabels.map { label2RecordsNumber[it]!! }.toIntArray())
            .with("n_loi_loaded", loiLabels.map { label2LoadedRegionsNumber[it]!! }.toIntArray())
            .with("n_loi_merged", loiRangesNumber)
            .with("n_regions", IntArray(loiLabels.size) { inputRegions.size })
            .with("metric", loiLabels.map { "${metric.column}($metricArgString)" }.toTypedArray())
            .with("input_${metric.column}", metricValuesForInput.toLongArray())
            .with("sampled_median_${metric.column}", sampledSetsMetricMedian.toIntArray())
            .with("sampled_mean_${metric.column}", sampledSetsMetricMean.toDoubleArray())
            .with("sampled_var_${metric.column}", sampledSetsMetricVar.toDoubleArray())
            .with("sampled_sets_n", IntArray(loiLabels.size) { simulationsNumber })
            .with("pValue", pValues)
            .with("qValue", qValues)
            .reorder("pValue")
    }

    private fun <T> sampleRegions(
        simulationsNumber: Int,
        intersectionFilter: LocationsList<out RangesList>?,
        srcLoci: List<Location>,
        background: T,
        gq: GenomeQuery,
        threadsNum: Int,
        withReplacement: Boolean,
        samplingFun: (GenomeQuery, List<ChromosomeRange>, T, Int, Boolean) -> List<Location>
    ): List<List<LocationsList<out RangesList>>> {
        // Sample:
        val progress = Progress { title = "Loci Sampling" }.bounded(simulationsNumber.toLong())
        val loci = srcLoci.map { it.toChromosomeRange() }

        // Compute different simulations in parallel, they are independent
        val sampled = IntStream.range(0, simulationsNumber).parallel().mapToObj { _ ->
            val randLoci = samplingFun(gq, loci, background, maxRetries, withReplacement)

            progress.report()

            if (withReplacement) {
                LocationsSortedList.create(gq, randLoci)
            } else {
                // XXX: shuffled regions not intersect by def of our shuffle procedure
                LocationsMergingList.create(gq, randLoci)
            }

        }.collect(Collectors.toList())
        progress.done()

        val chunked = sampled.chunked(threadsNum)
        if (intersectionFilter == null) {
            return chunked
        }

        // Apply Filters:
        return chunked.map { chunk ->
            chunk.map { ll ->
                val filtered = ll.intersectRanges(intersectionFilter)
                LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
            }
        }
    }


    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStats::class.java)

        /**
         * Reads locations from a given path, ignoring strand information.
         *
         * @param path The path to the file containing the locations.
         * @param gq The genome query.
         * @param mergeOverlapped If true,  overlapped regions will be merged
         *
         * @return A triple containing the locations list and the number of loaded locations before merge and
         *   the records number in initial file
         */
        fun readLocationsIgnoringStrand(
            path: Path, gq: GenomeQuery, mergeOverlapped: Boolean
        ): Triple<LocationsList<out RangesList>, Int, Int> {
            val (locations, recordsNumber) = readLocationsIgnoringStrand(path, gq)
            return when {
                mergeOverlapped -> {
                    val mergedLocations = LocationsMergingList.create(gq, locations)
                    if (locations.size != mergedLocations.size) {
                        LOG.info("$path: ${locations.size} regions merged to ${mergedLocations.size}")
                    }
                    Triple(mergedLocations, locations.size, recordsNumber)
                }
                else -> Triple(LocationsSortedList.create(gq, locations), locations.size, recordsNumber)
            }
        }

        fun readLocationsIgnoringStrand(
            path: Path,
            gq: GenomeQuery,
            bedFormat: BedFormat = BedFormat.auto(path)
        ): Pair<List<Location>, Int> {
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
                LOG.warn(
                    "$path: Loaded $pnt % (${loci.size} of $recordsNumber) locations." +
                            " Ignored chromosomes: ${ignoredChrs.size}. For more details use debug option."
                )
            }
            if (ignoredChrs.isNotEmpty()) {
                LOG.debug("$path: Ignored chromosomes: $ignoredChrs")
            }
            return loci to recordsNumber
        }

        /**
         * Loads complementary regions to the given masked genome regions filter.
         *
         * @param genomeMaskedLociPath The path to the file containing the masked genome loci. File is TAB-separated BED with at least 3 fields. Strand is ignored.
         * @param gq The genome query object.
         * @return A list of locations representing the complementary regions to the masked genome regions filter.
         */
        fun loadComplementaryToMaskedGenomeRegionsFilter(
            genomeMaskedLociPath: Path,
            gq: GenomeQuery
        ): LocationsMergingList {
            val maskedGenome = readLocationsIgnoringStrand(genomeMaskedLociPath, gq).first
            LOG.info("Genome masked loci: ${maskedGenome.size.formatLongNumber()} regions")
            val maskedGenomeLocations = LocationsMergingList.create(gq, maskedGenome)
            LOG.info("Genome masked loci (merged): ${maskedGenomeLocations.size.formatLongNumber()} regions")

            // complementary regions
            return maskedGenomeLocations.apply { rl, chr, _strand -> rl.complementaryRanges(chr.length) }
        }

        /**
         * Loads the genome allowed locations filter from the specified genome allowed loci file and genome query.
         *
         * @param genomeAllowedLociPath The path to the genome allowed loci file. File is TAB-separated BED with at least 3 fields. Strand is ignored.
         * @param gq The genome query.
         * @return The loaded genome allowed locations filter as a LocationsMergingList.
         */
        fun loadGenomeAllowedLocationsFilter(
            genomeAllowedLociPath: Path,
            gq: GenomeQuery
        ): LocationsMergingList {
            val allowedGenome = genomeAllowedLociPath.let {
                readLocationsIgnoringStrand(it, gq).first
            }

            LOG.info("Genome allowed loci: ${allowedGenome.size.formatLongNumber()} regions")
            val allowedGenomeLocations = LocationsMergingList.create(gq, allowedGenome)

            LOG.info("Genome allowed loci (merged): ${allowedGenomeLocations.size.formatLongNumber()} regions")
            return allowedGenomeLocations
        }

        fun makeAllowedRegionsFilter(
            genomeMaskedLociPath: Path?,
            genomeAllowedLociPath: Path?,
            gq: GenomeQuery
        ): LocationsMergingList? {
            // complementary to masked list
            val genomeMaskedComplementaryFilter = genomeMaskedLociPath?.let {
                loadComplementaryToMaskedGenomeRegionsFilter(it, gq)
            }

            // allowed list
            val genomeAllowedFilter = genomeAllowedLociPath?.let { loadGenomeAllowedLocationsFilter(it, gq) }

            return when {
                genomeAllowedFilter == null -> genomeMaskedComplementaryFilter
                genomeMaskedComplementaryFilter == null -> genomeAllowedFilter
                else -> genomeAllowedFilter.intersectRanges(genomeMaskedComplementaryFilter) as LocationsMergingList
            }
        }


        fun calcMetric(
            sampledLoci: LocationsList<out RangesList>,
            loi: LocationsList<out RangesList>,
            aSetIsSampled: Boolean,
            metric: RegionsMetric
        ): Long {
            val a = if (aSetIsSampled) sampledLoci else loi
            val b = if (aSetIsSampled) loi else sampledLoci

            // XXX: at the moment only 'overlap' is used here, i.e. integer metric
            return metric.calcMetric(a, b).toLong()
        }

        fun calcMetricForSampled(
            sampledRegions: List<List<LocationsList<out RangesList>>>,
            lociToTest: LocationsList<out RangesList>,
            aSetIsLoi: Boolean,
            metric: RegionsMetric,
            metricValueForSrc: Long
        ): List<PerThreadStats> = sampledRegions.parallelStream().map { chunk ->
            var countGreater = 0 // count when random regions metric value is above given value
            var countBelow = 0 // count when random regions metric value is above given value

            val metricHist = IntHistogram()
            chunk.forEach { ll ->
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
    }
}

data class PerThreadStats(
    val countSetsWithMetricsAboveThr: Int,
    val countSetsWithMetricsBelowThr: Int,
    val metricHist: IntHistogram = IntHistogram()
)
