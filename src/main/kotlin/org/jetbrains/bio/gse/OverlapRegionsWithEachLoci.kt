package org.jetbrains.bio.gse

import joptsimple.OptionParser
import org.jetbrains.bio.BioinfToolsCLA
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.containers.toRangeSortedList
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream

object OverlapRegionsWithEachLoci {
    private val LOG = LoggerFactory.getLogger(OverlapRegionsWithEachLoci::class.java)

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
                listOf("l", "loi"),
                "Loci of interest (loi) file path in TAB separated BED, BED3 or BED4 format. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.bedtoolsValidFile(minBedSpecFields = 4))
                .required()

            acceptsAll(
                listOf("r", "regions"),
                "Regions *.bed file or folder with *.bed regions files. Strand is ignored."
            )
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.noCheck())
                .required()

            acceptsAll(
                listOf("regions-filter"),
                "Filter regions files ending with requested suffix string (not regexp)"
            )
                .withRequiredArg()

            acceptsAll(listOf("regions-flanked"), "Flank regions at both sides (non-negative dist in bp)")
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
                "Merge overlapped locations while loading files if applicable. Will speedup computations."
            )

            // Logging level:
            acceptsAll(listOf("d", "debug"), "Print all the debug info")

            val tool = BioinfToolsCLA.Tools.OVERLAP_REGIONS_WITH_EACH_LOCI
            parse(args, description = tool.description) { options ->
                BioinfToolsCLA.configureLogging("quiet" in options, "debug" in options)
                LOG.info("Tool [${tool.command}]: ${tool.description} (vers: ${BioinfToolsCLA.version()})")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                LOG.info("CHROM.SIZES: $chromSizesPath")

                val chrMapStr = options.valuesOf("chrmap") as List<String>
                LOG.info("CHROM MAPPING: $chrMapStr")

                val genome = Genome.get(
                    chromSizesPath,
                    chrAltName2CanonicalMapping = parseChrNamesMapping(chrMapStr)
                )
                LOG.info("GENOME BUILD: ${genome.build}")

                val srcLoci = options.valueOf("loi") as Path
                LOG.info("SOURCE_LOCI: $srcLoci")

                val regionsFolderPath = options.valueOf("regions") as Path
                LOG.info("REGIONS: $regionsFolderPath")
                val regionsNameSuffix = options.valueOf("regions-filter") as String?
                LOG.info("REGIONS_FNAME_SUFFIX: ${regionsNameSuffix ?: "N/A"}")

                val aSetFlankedBothSides = options.valueOf("regions-flanked") as Int
                LOG.info("REGIONS FLANKED: $aSetFlankedBothSides")

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
                    srcLoci, regionsFolderPath, genome, outputFolder, metric,
                    mergeOverlapped, regionsNameSuffix
                )
            }
        }
    }

    fun doCalculations(
        loiPath: Path,
        regionsPath: Path,
        genome: Genome,
        outputFolderPath: Path,
        metric: RegionsMetric,
        mergeOverlapped: Boolean,
        regionsNameSuffix: String?
    ) {
        outputFolderPath.createDirectories()

        val threadsNumber = parallelismLevel();
        val regionLabelAndLociToTest: List<List<Pair<String, LocationsList<out RangesList>>>> =
            EnrichmentInRegions.collectRegionsFrom(
                if (regionsPath.isDirectory) Files.list(regionsPath) else Stream.of(regionsPath),
                genome, mergeOverlapped, null, regionsNameSuffix
            ).chunked(threadsNumber)

        val regionsSetsToTestNumber = regionLabelAndLociToTest.sumOf { it.size }
        LOG.info("Regions sets to test: $regionsSetsToTestNumber")
        require(regionsSetsToTestNumber > 0) {
            "No regions files passed file suffix filter."
        }
        val gq = genome.toQuery()
        val (sourceLoci, sourceLociBedFormat) = readNamedLocationsIgnoringStrand(loiPath, gq)
        require(sourceLoci.isNotEmpty()) {
            "Loci file is empty or all loci were masked."
        }
        LOG.info("Loaded ${sourceLoci.size} loci, bed format: $sourceLociBedFormat")
        val sourceLoci2RangeList = sourceLoci.map { it to it.location.toRangeSortedList() }

        val progress = Progress { title = "Overlap(region set, locus) progress" }.bounded(
            regionsSetsToTestNumber.toLong()
        )
        regionLabelAndLociToTest.parallelStream().forEach { chunk ->
            for ((regionTypeLabel, regionLociToTest) in chunk) {
                progress.report()

                val values = LongArray(sourceLoci2RangeList.size)
                sourceLoci2RangeList.forEachIndexed { i, (namedLocus, lociAsRangeList) ->
                    val metricValueForSrc = regionLociToTest.calcAdditiveMetricDouble(
                        lociAsRangeList, namedLocus.location.chromosome, namedLocus.location.strand,
                        metric::calcMetric
                    ).toLong()
                    values[i] = metricValueForSrc
                }

                val path = outputFolderPath / "${regionTypeLabel}_${metric.column}.tsv"
                // Save:
                DataFrame()
                    .with("chr", sourceLoci2RangeList.map { it.first.location.chromosome.name }.toTypedArray())
                    .with("startOffset", sourceLoci2RangeList.map { it.first.location.startOffset }.toIntArray())
                    .with("endOffset", sourceLoci2RangeList.map { it.first.location.endOffset }.toIntArray()).let {
                        if (sourceLociBedFormat.fieldsNumber >= 4) {
                            it.with(
                                "name",
                                sourceLoci2RangeList.map { (it.first as NamedLocation).name }.toTypedArray()
                            )
                        } else {
                            it
                        }
                    }.with(metric.column, values)
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

        val loci = bedFormat.parse(path) { bedParser ->
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
        return loci to bedFormat
    }

}

data class NamedLocation(val name: String, override val location: Location) : LocationAware