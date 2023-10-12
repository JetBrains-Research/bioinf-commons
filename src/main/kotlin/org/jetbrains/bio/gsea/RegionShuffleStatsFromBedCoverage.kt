package org.jetbrains.bio.gsea

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.OverlapNumberMetric
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.sampling.shuffleChromosomeRanges
import org.jetbrains.bio.util.formatLongNumber
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.ThreadLocalRandom
import kotlin.random.asKotlinRandom

class RegionShuffleStatsFromBedCoverage(
    simulationsNumber: Int,
    chunkSize: Int,
    regionSetMaxRetries: Int,
    singleRegionMaxRetries: Int
) : RegionShuffleStats(simulationsNumber, chunkSize, regionSetMaxRetries, singleRegionMaxRetries) {
    fun calcStatistics(
        gq: GenomeQuery,
        inputRegionsPath: Path,
        backgroundPath: Path?,
        loiInfos: List<LoiInfo>,
        outputFolderPath: Path?,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis,
        aSetIsRegions: Boolean,
        mergeOverlapped: Boolean,
        truncateFilter: LocationsList<out RangesList>?,
        genomeAllowedAreaFilter: LocationsMergingList?,
        genomeMaskedAreaFilter: LocationsMergingList?,
        mergeInputRegionsToBg: Boolean,
        samplingWithReplacement: Boolean,
        backgroundIsMethylome: Boolean,
        zeroBasedBackground: Boolean,
        bedBgFlnk: Int,
    ) = doCalcStatistics(
        inputRegionsPath,
        gq,
        genomeAllowedAreaFilter,
        genomeMaskedAreaFilter,
        loiInfos,
        outputFolderPath,
        metric,
        hypAlt,
        aSetIsRegions,
        mergeOverlapped,
        truncateFilter,
        samplingWithReplacement = samplingWithReplacement,
        backgroundProvider = { genomeQuery, inputRegionsFiltered ->
            backgroundProviderFun(
                inputRegionsFiltered = inputRegionsFiltered,
                backgroundPath = backgroundPath,
                mergeInputRegionsToBg = mergeInputRegionsToBg,
                genomeAllowedAreaFilter = genomeAllowedAreaFilter,
                genomeMaskedAreaFilter = genomeMaskedAreaFilter,
                gq = genomeQuery,
                backgroundIsMethylome=backgroundIsMethylome,
                zeroBasedBackground=zeroBasedBackground,
                bedBgFlnk=bedBgFlnk
            )
        },
        loiOverlapWithBgFun = { loiFiltered, background ->
            OverlapNumberMetric().calcMetric(loiFiltered, background).toInt()
        },
        samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
            val threadLocalRandom = ThreadLocalRandom.current().asKotlinRandom()
            val res = shuffleChromosomeRanges(
                genomeQuery,
                regions,
                background.asLocationSequence().map { it.toChromosomeRange() }.toList(),
                regionSetMaxRetries = regionSetMaxRetries,
                singleRegionMaxRetries = singleRegionMaxRetries,
                withReplacement = withReplacement,
                genomeMaskedAreaFilter = genomeMaskedAreaFilter,
                randomGen = threadLocalRandom
            )
            res.first.map { it.on(Strand.PLUS) } to res.second
        }
    )

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStatsFromBedCoverage::class.java)

        fun backgroundProviderFun(
            inputRegionsFiltered: List<Location>,
            backgroundPath: Path?,
            mergeInputRegionsToBg: Boolean,
            genomeAllowedAreaFilter: LocationsMergingList?,
            genomeMaskedAreaFilter: LocationsMergingList?,
            gq: GenomeQuery,
            backgroundIsMethylome: Boolean,
            zeroBasedBackground: Boolean,
            bedBgFlnk: Int,
        ): LocationsMergingList {
            // Here no need in intersecting background with `SharedSamplingOptions.truncateRangesToSpecificLocation`

            var bgRegions = if (backgroundPath != null) {
                val bgLoci = if (backgroundIsMethylome) {
                    val methylomeCov = RegionShuffleStatsFromMethylomeCoverage.loadMethylomeCovBackground(
                        backgroundPath,
                        zeroBasedBackground,
                        gq
                    )
                     SamplingMethylationValidation.makeBEDBackgroundFromMethylome(
                        gq, methylomeCov, bedBgFlnk,
                        if (mergeInputRegionsToBg) inputRegionsFiltered else emptyList()
                    )
                } else {
                    require(zeroBasedBackground) { "1-based regions not yet supported here" }
                    require(bedBgFlnk == 0) { "Flanking not yet supported here" }
                    val bgLocations = readLocationsIgnoringStrand(backgroundPath, gq).first.toMutableList()
                    if (mergeInputRegionsToBg) {
                        // merge source loci into background
                        bgLocations.addAll(inputRegionsFiltered)
                    }
                    LocationsMergingList.create(
                        gq,
                        bgLocations
                    )
                }
                LOG.info("Background regions: ${bgLoci.size} regions")

                inputRegionsFiltered.forEach {
                    require(bgLoci.includes(it)) {
                        "Background $backgroundPath regions are required to include all loci of interest, but the " +
                                "region is missing in bg: ${it.toChromosomeRange()}"
                    }
                }
                bgLoci
            } else {
                // whole genome as bg
                val bgLoci = LocationsMergingList.create(
                    gq, gq.get().map { Location(0, it.length, it) }
                )
                LOG.info("Background regions: Using whole genome as background. ${bgLoci.size} regions")

                bgLoci
            }

            if (genomeAllowedAreaFilter != null) {
                LOG.info("Background regions: Applying allowed area filter...")

                val beforeFilteringBgRegionsNum = bgRegions.size

                // In bedtools shuffle, it is important to find a random point in an allowed background, and
                // a sampled region could just intersect but not be included in the allowed area
                // this corresponds to sampling alg used for a BED-based background.
                // I.e. we need intersected BG only to place an anchor,
                // then the resulting region could go out of allowed genome and bg.

                // For the consistency with sampling results let's expect that input regions also intersect
                // the background and could be not fully included into bg:
                bgRegions = bgRegions.intersectRanges(genomeAllowedAreaFilter) as LocationsMergingList

                LOG.info("Background regions: Allowed regions: ${bgRegions.size.formatLongNumber()} regions of ${beforeFilteringBgRegionsNum.formatLongNumber()}")
            }

            if (genomeMaskedAreaFilter != null) {
                LOG.info("Background regions: Applying masked area filter...")

                val beforeFilteringBgRegionsNum = bgRegions.size

                // These points are anchor points, if they are in `genomeMaskedAreaFilter`, then resulting regions
                // will be dismissed after sampling => we could remove them now.
                // But later checks still required because region started in valid location could intersect some masked.
                bgRegions = bgRegions.intersectRanges(genomeMaskedAreaFilter.makeComplementary()) as LocationsMergingList

                LOG.info("Background regions: W/o masked regions: ${bgRegions.size.formatLongNumber()} regions of ${beforeFilteringBgRegionsNum.formatLongNumber()}")
            }
            return bgRegions
        }
    }
}