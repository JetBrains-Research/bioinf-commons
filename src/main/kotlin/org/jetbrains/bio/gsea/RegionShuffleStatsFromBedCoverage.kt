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
        intersectionFilter: LocationsList<out RangesList>?,
        genomeMaskedAreaPath: Path?,
        genomeAllowedAreaPath: Path?,
        mergeRegionsToBg: Boolean,
        samplingWithReplacement: Boolean
    ) = doCalcStatistics(
        gq,
        loiInfos,
        outputFolderPath,
        metric,
        hypAlt,
        aSetIsRegions,
        mergeOverlapped,
        intersectionFilter,
        samplingWithReplacement = samplingWithReplacement,
        inputRegionsAndBackgroundProvider = { genomeQuery ->
            loadInputRegionsAndBedLikeBackground(
                inputRegionsPath, backgroundPath, mergeRegionsToBg, genomeMaskedAreaPath, genomeAllowedAreaPath, genomeQuery
            )
        },
        loiOverlapWithBgFun = { loiFiltered, background ->
            OverlapNumberMetric().calcMetric(loiFiltered, background).toInt()
        },
        samplingFun = { genomeQuery, regions, background, regionSetMaxRetries, singleRegionMaxRetries, withReplacement ->
            val res = shuffleChromosomeRanges(
                genomeQuery,
                regions,
                background.asLocationSequence().map { it.toChromosomeRange() }.toList(),
                regionSetMaxRetries = regionSetMaxRetries,
                singleRegionMaxRetries = singleRegionMaxRetries,
                withReplacement = withReplacement
            )
            res.first.map { it.on(Strand.PLUS) } to res.second
        }
    )

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStatsFromBedCoverage::class.java)

        fun loadInputRegionsAndBedLikeBackground(
            inputRegionsPath: Path,
            backgroundRegionsPath: Path?,
            mergeRegionsToBg: Boolean,
            genomeMaskedAreaPath: Path?,
            genomeAllowedAreaPath: Path?,
            gq: GenomeQuery
        ): Pair<List<Location>, LocationsMergingList> {
            val inputRegions = readLocationsIgnoringStrand(inputRegionsPath, gq).first
            LOG.info("Input regions: ${inputRegions.size.formatLongNumber()} regions")

            val bgRegions = if (backgroundRegionsPath != null) {
                val bgLocations = readLocationsIgnoringStrand(backgroundRegionsPath, gq).first.toMutableList()
                if (mergeRegionsToBg) {
                    // merge source loci into background
                    bgLocations.addAll(inputRegions)
                }
                val bgLoci = LocationsMergingList.create(
                    gq,
                    bgLocations
                )
                LOG.info("Background regions: ${bgLoci.size} regions")

                inputRegions.forEach {
                    require(bgLoci.includes(it)) {
                        "Background $backgroundRegionsPath regions are required to include all loci of interest, but the " +
                                "region is missing in bg: ${it.toChromosomeRange()}"
                    }
                }
                bgLoci
            } else {
                // whole genome as bg
                val bgLoci = LocationsMergingList.create(
                    gq, gq.get().map { Location(0, it.length, it) }
                )
                LOG.info("Using whole genome as background. Background regions: ${bgLoci.size} regions")

                bgLoci
            }

            val allowedGenomeFilter = makeAllowedRegionsFilter(
                genomeMaskedAreaPath, genomeAllowedAreaPath, gq
            )

            @Suppress("FoldInitializerAndIfToElvis")
            if (allowedGenomeFilter == null) {
                return inputRegions to bgRegions
            }

            LOG.info("Applying allowed regions filters...")
            val allowedBg = bgRegions.intersectRanges(allowedGenomeFilter) as LocationsMergingList
            val allowedInputRegions = inputRegions.filter { allowedGenomeFilter.includes(it) }

            LOG.info("Background regions (all filters applied): ${allowedBg.size.formatLongNumber()} regions of ${bgRegions.size.formatLongNumber()}")
            LOG.info("Input regions (all filters applied): ${allowedInputRegions.size.formatLongNumber()} regions of ${inputRegions.size.formatLongNumber()}")
            return allowedInputRegions to allowedBg
        }
    }
}