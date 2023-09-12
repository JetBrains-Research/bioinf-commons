package org.jetbrains.bio.gsea

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.containers.toRangeSortedList
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream

object OverlapLoiWithEachRegion {
    private val LOG = LoggerFactory.getLogger(OverlapLoiWithEachRegion::class.java)

    @JvmStatic
    fun main(args: Array<String>) {
        with(OptionParser()) {
            // Chrom.size
            acceptsAll(
                listOf("cs", "chrom.sizes"),
                "Chromosome sizes path, can be downloaded at" +
                        " http://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes"
            ).withRequiredArg()
                .required()
                .withValuesConvertedBy(PathConverter.exists())


            // Chromosome additional mapping
            acceptsAll(
                listOf("chrmap"),
                "Comma separated list of `chrInput:chrGenome` pairs, where `chrInput` is chromosome name in " +
                        "input file and `chrGenome` chromosome name in *.chrom.sizes file."
            ).withRequiredArg()
                .withValuesSeparatedBy(",")

            acceptsAll(
                listOf("r", "regions"),
                "Regions to test file. File is TAB-separated BED with at least 3 fields. Strand is ignored." +
                        " E.g. regions could be DMRs."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 3))
                .required()

            acceptsAll(
                listOf("l", "loi"),
                "Loci of interest (LOI) to check enrichment. Could be single *.bed file or folder with *.bed" +
                        " files. Strand is ignored. E.g. loi could be CGIs file."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            acceptsAll(
                listOf("loi-filter"),
                "Filter LOI files ending with requested suffix string (not regexp)"
            )
                .withRequiredArg()

            acceptsAll(listOf("a-flanked"), "Flank 'a' ranges at both sides (non-negative dist in bp)")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(0)

            // Output
            acceptsAll(
                listOf("o", "output"),
                "Output folder path"
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.directory())

            acceptsAll(listOf("parallelism"), "parallelism level")
                .withRequiredArg()
                .ofType(Int::class.java)

//            accepts(
//                "detailed",
//                "Generated detailed stats report with metric value for each shuffle"
//            )

            acceptsAll(
                listOf("m", "merge"),
                "Merge overlapped locations (e.g regions, loi) while loading files if applicable. Will speedup computations."
            )

            // Logging level:
            acceptsAll(listOf("d", "debug"), "Print all the debug info")

            val tool = BioinfToolsCLA.Tools.OVERLAP_LOI_WITH_EACH_REGION
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                LOG.info("CHROM.SIZES: $chromSizesPath")

                val chrMapStr = options.valuesOf("chrmap").filterIsInstance<String>()
                LOG.info("CHROM MAPPING: $chrMapStr")

                val genome = Genome.get(
                    chromSizesPath,
                    chrAltName2CanonicalMapping = parseChrNamesMapping(chrMapStr)
                )
                LOG.info("GENOME BUILD: ${genome.build}")

                val inputRegions = options.valueOf("regions") as Path
                LOG.info("INPUT REGIONS: $inputRegions")

                val loiFolderPath = options.valueOf("loi") as Path
                LOG.info("LOI TO TEST: $loiFolderPath")

                val loiNameSuffix = options.valueOf("loi-filter") as String?
                LOG.info("LOI FNAME SUFFIX: ${loiNameSuffix ?: "N/A"}")

                val aSetFlankedBothSides = options.valueOf("a-flanked") as Int
                LOG.info("A FLANKED: $aSetFlankedBothSides")

                val outputFolder = options.valueOf("output") as Path
                LOG.info("OUTPUT_FOLDER: $outputFolder")

                val parallelism = options.valueOf("parallelism") as Int?
                LOG.info("THREADS: $parallelism")
                configureParallelism(parallelism)

//                val detailedReport = options.has("detailed")
//                LOG.info("DETAILED_REPORT_FLAG: $detailedReport")

                val mergeOverlapped = options.has("merge")
                LOG.info("MERGE OVERLAPPED: $mergeOverlapped")

                val metric = OverlapNumberMetric(aSetFlankedBothSides)
                LOG.info("METRIC: ${metric.column}")

                doCalculations(
                    inputRegions, loiFolderPath, genome, outputFolder, metric,
                    mergeOverlapped, loiNameSuffix
                )
            }
        }
    }

    fun doCalculations(
        inputRegionsPath: Path,
        loiFolderPath: Path,
        genome: Genome,
        outputFolderPath: Path,
        metric: RegionsMetric,
        mergeOverlapped: Boolean,
        loiNameSuffix: String?
    ) {
        outputFolderPath.createDirectories()
        val gq = genome.toQuery()

        val threadsNumber = parallelismLevel()
        val filesStream = if (loiFolderPath.isDirectory) Files.list(loiFolderPath) else Stream.of(loiFolderPath)
        val chunkedLoiInfos: List<List<LoiInfo>>  = EnrichmentInLoi.collectLoiFrom(
            filesStream, gq, mergeOverlapped, null, loiNameSuffix
        ).chunked(threadsNumber)

        val totalSetsToTest = chunkedLoiInfos.sumOf { it.size }
        LOG.info("LOI sets to test: $totalSetsToTest")
        require(totalSetsToTest > 0) {
            "No loi files passed file suffix filter."
        }

        val (inputRegions, inputRegionsBedFormat) = readNamedLocationsIgnoringStrand(inputRegionsPath, gq)
        require(inputRegions.isNotEmpty()) {
            "Regions file is empty."
        }
        LOG.info("Loaded ${inputRegions.size.formatLongNumber()} regions, bed format: $inputRegionsBedFormat")
        val regions2RangesList = inputRegions.map { it to it.location.toRangeSortedList() }

        val progress = Progress { title = "Overlap(loi set, each region) progress" }.bounded(
            totalSetsToTest.toLong()
        )
        chunkedLoiInfos.parallelStream().forEach { chunk ->
            for (info in chunk) {
                progress.report()

                val values = LongArray(regions2RangesList.size)
                regions2RangesList.forEachIndexed { i, (namedRegion, regionAsRangeList) ->
                    val metricValueForSrc = info.lociFiltered.calcAdditiveMetricDouble(
                        regionAsRangeList, namedRegion.location.chromosome, namedRegion.location.strand,
                        metric::calcMetric
                    ).toLong()
                    values[i] = metricValueForSrc
                }

                val path = outputFolderPath / "${info.label}_${metric.column}.tsv"
                // Save:
                DataFrame()
                    .with(
                        "chr",
                        regions2RangesList.map { (regin, _) -> regin.location.chromosome.name }.toTypedArray()
                    )
                    .with(
                        "startOffset",
                        regions2RangesList.map { (regin, _) -> regin.location.startOffset }.toIntArray()
                    )
                    .with(
                        "endOffset",
                        regions2RangesList.map { (regin, _) -> regin.location.endOffset }.toIntArray()
                    ).let { df ->
                        when {
                            inputRegionsBedFormat.fieldsNumber >= 4 -> {
                                df.with(
                                    "name",
                                    regions2RangesList.map { (regin, _) ->
                                        (regin as NamedLocation).name
                                    }.toTypedArray()
                                )
                            }
                            else -> df
                        }
                    }
                    .with(metric.column, values)
                    .save(path)
                LOG.info("Report saved to: $path")
            }
        }
        progress.done("Done.")
    }

    fun readNamedLocationsIgnoringStrand(
        path: Path,
        gq: GenomeQuery,
        bedFormat: BedFormat = BedFormat.auto(path)
    ): Pair<List<LocationAware>, BedFormat> {
        var recordsNumber = 0
        val ignoredChrs = mutableListOf<String>()

        val locations = bedFormat.parse(path) { bedParser ->
            bedParser.mapNotNull {
                recordsNumber++

                val chr = gq[it.chrom]
                if (chr == null) {
                    ignoredChrs.add(it.chrom)
                    // skip all unmapped contigs, etc
                    null
                } else {
                    val location = Location(it.start, it.end, chr)

                    if (bedFormat.fieldsNumber >= 4) {
                        val name = it.unpack(4, false).name
                        NamedLocation(name, location)
                    } else {
                        location
                    }
                }
            }
        }
        if (locations.size != recordsNumber) {
            val pnt = locations.size.asPercentOf(recordsNumber)
            LOG.warn(
                "$path: Loaded $pnt % (${locations.size} of $recordsNumber) locations." +
                        " Ignored chromosomes: ${ignoredChrs.size}. For more details use debug option."
            )
        }
        if (ignoredChrs.isNotEmpty()) {
            LOG.debug("{}: Ignored chromosomes: {}", path, ignoredChrs)
        }
        return locations to bedFormat
    }

}

data class NamedLocation(val name: String, override val location: Location) : LocationAware