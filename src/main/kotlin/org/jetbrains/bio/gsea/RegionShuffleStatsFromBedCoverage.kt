package org.jetbrains.bio.gsea

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.LocationsList
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.RangesList
import org.jetbrains.bio.genome.containers.intersection.RegionsMetric
import org.jetbrains.bio.genome.sampling.shuffleChromosomeRanges
import org.jetbrains.bio.util.formatLongNumber
import org.slf4j.LoggerFactory
import java.nio.file.Path

class RegionShuffleStatsFromBedCoverage(
    genome: Genome,
    simulationsNumber: Int,
    chunkSize: Int,
    maxRetries: Int
) : RegionShuffleStats(genome, simulationsNumber, chunkSize, maxRetries) {
    override fun calcStatistics(
        srcRegionsPath: Path,
        backgroundPath: Path?,
        regionLabelAndLociToTest: List<Pair<String, LocationsList<out RangesList>>>,
        outputFolderPath: Path?,
        metric: RegionsMetric,
        hypAlt: PermutationAltHypothesis,
        aSetIsLoi: Boolean,
        mergeOverlapped: Boolean,
        intersectionFilter: LocationsList<out RangesList>?,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        addLoiToBg: Boolean,
        samplingWithReplacement: Boolean
    ) = doCalcStatistics(
        regionLabelAndLociToTest,
        outputFolderPath,
        metric,
        hypAlt,
        aSetIsLoi,
        mergeOverlapped,
        intersectionFilter,
        srcLociAndBackgroundProvider = { gq ->
            loadSrcLociAndBedLikeBackground(
                srcRegionsPath, backgroundPath, addLoiToBg, genomeMaskedLociPath, genomeAllowedLociPath, gq
            )
        },
        samplingFun = { gq, regions, background, maxRetries, withReplacement ->
            shuffleChromosomeRanges(
                gq,
                regions,
                background.asLocationSequence().map { it.toChromosomeRange() }.toList(),
                maxRetries = maxRetries,
                withReplacement = withReplacement
            ).map { it.on(Strand.PLUS) }
        }
    )

    private fun loadSrcLociAndBedLikeBackground(
        srcRegionsPath: Path, backgroundRegionsPath: Path?,
        addLoiToBg: Boolean,
        genomeMaskedLociPath: Path?,
        genomeAllowedLociPath: Path?,
        gq: GenomeQuery
    ): Pair<List<Location>, LocationsMergingList> {
        val sourceLoci = readLocationsIgnoringStrand(srcRegionsPath, gq)
        LOG.info("Source loci: ${sourceLoci.size} regions")

        val bgLoci = if (backgroundRegionsPath != null) {
            val bgLocations = readLocationsIgnoringStrand(backgroundRegionsPath, gq).toMutableList()
            if (addLoiToBg) {
                // merge source loci into background
                bgLocations.addAll(sourceLoci)
            }
            val bgLoci = LocationsMergingList.create(
                gq,
                bgLocations
            )
            LOG.info("Background regions: ${bgLoci.size} regions")

            sourceLoci.forEach {
                require(bgLoci.includes(it)) {
                    "Background $backgroundRegionsPath regions are required to include all loci of interest, but the " +
                            "loci is missing in bg: ${it.toChromosomeRange()}"
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
            genomeMaskedLociPath, genomeAllowedLociPath, gq
        )

        val (allowedBgList, allowedSourceLoci) = if (allowedGenomeFilter == null) {
            bgLoci to sourceLoci
        } else {
            LOG.info("Applying allowed regions filters...")
            val allowedBg = bgLoci.intersectRanges(allowedGenomeFilter) as LocationsMergingList
            val allowedSourceLoci = sourceLoci.filter { allowedGenomeFilter.includes(it) }

            LOG.info("Background regions (all restrictions applied): ${allowedBg.size.formatLongNumber()} regions")
            LOG.info("Source loci (all restrictions applied): ${allowedSourceLoci.size.formatLongNumber()} regions")
            allowedBg to allowedSourceLoci
        }

        return allowedSourceLoci to allowedBgList
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(RegionShuffleStatsFromBedCoverage::class.java)
    }
}