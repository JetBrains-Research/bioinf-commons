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
import org.jetbrains.bio.gsea.EnrichmentInLoi.processInputRegions
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
 * @param regionSetMaxRetries Max retries when cannot shuffle whole set of regions from background
 * @param singleRegionMaxRetries Max retries when cannot shuffle one region from background
 *
 * Loci file path will be added to background
 */
abstract class RegionShuffleStats(
    protected  val simulationsNumber: Int,
    protected  val chunkSize: Int,
    protected  val regionSetMaxRetries: Int,
    protected  val singleRegionMaxRetries: Int
) {

    protected fun <T> doCalcStatistics(
        inputRegionsPath: Path,
        gq: GenomeQuery,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
        loiInfos: List<LoiInfo>,
        outputFolderPath: Path? = null,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis = PermutationAltHypothesis.GREATER,
        aSetIsRegions: Boolean = true,
        mergeOverlapped: Boolean = true,
        truncateFilter: LocationsList<out RangesList>? = null,
        samplingWithReplacement: Boolean = false,
        backgroundProvider: (GenomeQuery, List<Location>) -> T, // by genome query & input regions
        loiOverlapWithBgFun: (LocationsList<out RangesList>, T) -> Int,
        samplingFun: (GenomeQuery, List<ChromosomeRange>, T, Int, Int, Boolean) -> Pair<List<Location>, IntHistogram>
    ): DataFrame {
        outputFolderPath?.createDirectories()
        val dumpDetails = outputFolderPath != null

        val inputRegionsFiltered = processInputRegions(
            inputRegionsPath, gq,
            genomeAllowedAreaFilter = genomeAllowedAreaFilter,
            genomeMaskedAreaFilter = genomeMaskedAreaFilter
        )
        require(inputRegionsFiltered.isNotEmpty()) {
            "Regions file is empty or all regions were removed by filters."
        }
        val inputRegionsListFiltered = when {
            mergeOverlapped -> LocationsMergingList.create(gq, inputRegionsFiltered)
            else -> LocationsSortedList.create(gq, inputRegionsFiltered)
        }
        
        val bgRegionsList = backgroundProvider(gq, inputRegionsFiltered)
        

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

        LOG.info("Calc LOI overlap with BG: ${loiInfos.size}")
        val infoOverlapWithBg = loiInfos.map { info -> loiOverlapWithBgFun(info.lociFiltered, bgRegionsList) }.toIntArray()

        LOG.info("Do overrepresented check...")
        val progress = Progress { title = "Over/Under-representation check progress (all chunks)" }.bounded(
            nChunks.toLong()
        )

        val finalRegionAttemptsHist = IntHistogram()
        (0 until nChunks).forEach { chunkId ->
            val start = chunkId * chunkSize
            val end = minOf(simulationsNumber, (chunkId + 1) * chunkSize)
            val simulationsInChunk = end - start

            LOG.info(
                "Simulations: Chunk [${chunkId + 1} of $nChunks], simulations " +
                        "${start.formatLongNumber()}..${end.formatLongNumber()} of ${simulationsNumber.formatLongNumber()}"
            )

            val (sampledRegions: List<List<LocationsList<out RangesList>>>, regionAttemptsHist) =
                sampleRegions(
                    simulationsInChunk, regionSetMaxRetries, singleRegionMaxRetries, truncateFilter, inputRegionsFiltered, bgRegionsList,
                    gq,
                    parallelismLevel(),
                    samplingFun=samplingFun,
                    withReplacement = samplingWithReplacement
                )
            finalRegionAttemptsHist += regionAttemptsHist

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

                val metricValueForInput = calcMetric(inputRegionsListFiltered, info.lociFiltered, aSetIsRegions, metric)

                val metricValueForSampled =
                    calcMetricForSampled(sampledRegions, info.lociFiltered, aSetIsRegions, metric, metricValueForInput)

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

        // Save attempts Histogram
        if (outputFolderPath != null) {
            val path = outputFolderPath / "sampled_${simulationsNumber}.region_attempts.hist.tsv"
            finalRegionAttemptsHist.save(path, metric_col = "attempt")
            LOG.info("[DONE]: $path")
        }

        // Save results to DataFrame:
        val loiLabels = loiInfos.map { it.label }
        val loiFilteredAndMergedRangesNumber = loiInfos.map { it.lociFiltered.size }.toIntArray()
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
            .with("loi", loiLabels.toTypedArray()) // #1
            .with("n_loi_records", loiLabels.map { label2RecordsNumber[it]!! }.toIntArray())  // #2
            .with("n_loi_loaded", loiLabels.map { label2LoadedRegionsNumber[it]!! }.toIntArray()) // #3
            .with("n_loi_filtered_merged", loiFilteredAndMergedRangesNumber) // #4
            .with("n_regions", IntArray(loiLabels.size) { inputRegionsFiltered.size }) // #5
            .with("metric", loiLabels.map { "${metric.column}($metricArgString)" }.toTypedArray()) // #6
            .with("n_regions_filtered", IntArray(loiLabels.size) { inputRegionsListFiltered.size }) // #7
            .with("loi_filtered_merged_bg_overlap", infoOverlapWithBg) // #8
            .with("input_${metric.column}", metricValuesForInput.toLongArray()) // #9
            .with("sampled_median_${metric.column}", sampledSetsMetricMedian.toIntArray()) // #10
            .with("sampled_mean_${metric.column}", sampledSetsMetricMean.toDoubleArray()) // #11
            .with("sampled_var_${metric.column}", sampledSetsMetricVar.toDoubleArray()) // #12
            .with("sampled_sets_n", IntArray(loiLabels.size) { simulationsNumber }) // #13
            .with("test_H1", loiLabels.map { hypAlt.name }.toTypedArray()) // #14
            .with("pValue", pValues) // #15
            .with("qValue", qValues) // #16
            .with("ratio_obs_2_input",  metricValuesForInput.map { obs -> obs.toFloat() / inputRegionsListFiltered.size }.toFloatArray()) // #17
            .with("ratio_obs_2_exp",  metricValuesForInput.zip(sampledSetsMetricMedian).map { (obs, n) -> obs.toFloat() / n }.toFloatArray()) // #18
            .reorder("qValue")
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStats::class.java)

        fun <T> sampleRegions(
            simulationsNumber: Int,
            regionSetMaxRetries: Int,
            singleRegionMaxRetries: Int,
            truncateFilter: LocationsList<out RangesList>?,
            srcLoci: List<Location>,
            background: T,
            gq: GenomeQuery,
            threadsNum: Int,
            withReplacement: Boolean,
            samplingFun: (GenomeQuery, List<ChromosomeRange>, T, Int, Int, Boolean) -> Pair<List<Location>, IntHistogram>
        ): Pair<List<List<LocationsList<out RangesList>>>, IntHistogram> {
            // Sample:
            val progress = Progress { title = "Loci Sampling" }.bounded(simulationsNumber.toLong())
            val loci = srcLoci.map { it.toChromosomeRange() }

            // Compute different simulations in parallel, they are independent
            val sampled = IntStream.range(0, simulationsNumber).parallel().mapToObj { _ ->
                val res = samplingFun(gq, loci, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement)
                val randLoci = res.first
                progress.report()

                val ll = if (withReplacement) {
                    LocationsSortedList.create(gq, randLoci)
                } else {
                    // XXX: shuffled regions not intersect by def of our shuffle procedure
                    LocationsMergingList.create(gq, randLoci)
                }
                ll to res.second

            }.collect(Collectors.toList())
            progress.done()

            val chunked = sampled.chunked(threadsNum)

            val finalHist = IntHistogram()
            chunked.forEach { chunk ->
                chunk.forEach { (_, hist) ->
                    finalHist += hist
                }
            }

            if (truncateFilter == null) {
                return chunked.map { chunk ->
                    chunk.map { (ll, hist) ->
                        finalHist += hist
                        ll
                    }
                } to finalHist
            }

            // Apply Filters:
            return chunked.map { chunk ->
                chunk.map { (ll, hist) ->
                    finalHist += hist
                    val filtered = ll.intersectRanges(truncateFilter)
                    LocationsMergingList.create(filtered.genomeQuery, filtered.locationIterator())
                }
            } to finalHist
        }

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
                LOG.debug("{}: Ignored chromosomes: {}", path, ignoredChrs)
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
            return maskedGenomeLocations.makeComplementary()
        }

        /**
         * Loads the genome allowed locations filter from the specified genome allowed loci file and genome query.
         *
         * @param lociPath The path to the genome allowed loci file. File is TAB-separated BED with at least 3 fields. Strand is ignored.
         * @param gq The genome query.
         * @return The loaded genome allowed locations filter as a LocationsMergingList.
         */
        fun readGenomeAreaFilter(
            lociPath: Path,
            gq: GenomeQuery
        ): LocationsMergingList {
            val (loci, recordsNumber) = lociPath.let {
                readLocationsIgnoringStrand(it, gq)
            }

            LOG.info("Genome area filter [${lociPath.fileName}]: ${loci.size.formatLongNumber()} regions of $recordsNumber lines")
            val ll = LocationsMergingList.create(gq, loci)

            LOG.info("Genome area filter [${lociPath.fileName}]: ${ll.size.formatLongNumber()} merged regions")
            return ll
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
