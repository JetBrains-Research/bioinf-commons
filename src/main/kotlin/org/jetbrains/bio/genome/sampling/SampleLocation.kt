package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.big.ExtendedBedEntry
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.sequence.Nucleotide
import java.nio.file.Files
import java.nio.file.Path
import java.util.*

class SamplingError(message: String) : Exception(message)

private val RANDOM = Random()
private const val SAMPLE_ATTEMPTS_THRESHOLD = 1000

/**
 * Shuffles bed regions positions inside chromosomes. All bed regions stays in same chromosome where they were
 * and their length remains unchanged.
 *
 * @return path to newly generated Bed file.
 */
fun randomizeBedRegions(bedFilePath: Path, genomeQuery: GenomeQuery): Path {
    val randomBedPath = Files.createTempFile("randomRegions", ".bed")
    BedFormat().print(randomBedPath).use {
        BedFormat.auto(bedFilePath).parse(bedFilePath) { bedParser ->
            bedParser.distinct()
                .filter { it.chrom in genomeQuery }
                .map { it.chrom to (it.end - it.start) }
                .groupBy { it.first }
                .mapValues { it.value.map { p -> p.second } }
                .toList().parallelStream()
                .map {
                    val chr = Chromosome(genomeQuery.genome, it.first)
                    sampleLocations(chr, it.second, 0, chr.length, false, Strand.PLUS)
                }
                .sequential() // required here, otherwise flatMap().forEach() seems to look concurrently
                .flatMap { it.stream() }
                .forEach { (start, end, chr, strand) ->
                    it.print(ExtendedBedEntry(chr.name, start, end, strand = strand.char))
                }
        }
    }
    return randomBedPath
}

/**
 * Generates list of given size containing random locations (which may intersect)
 * of given length within chromosome[leftBound, rightBound)
 */
fun sampleLocations(
    chromosome: Chromosome,
    lengthList: List<Int>,
    leftBound: Int,
    rightBound: Int,
    allowNs: Boolean = false,
    predefinedStrand: Strand? = null
): List<Location> {
    if (rightBound - leftBound <= (lengthList.maxOrNull() ?: 0)) {
        throw SamplingError("Not enough space for location length ${lengthList.maxOrNull()} within [$leftBound, $rightBound)")
    }

    val locations = arrayListOf<Location>()
    lengthList.forEach { length ->
        var location: Location
        var attempts = 0
        while (true) {
            val startOffset = RANDOM.nextInt(rightBound - leftBound - length)
            val strand = predefinedStrand ?: if (RANDOM.nextBoolean()) Strand.PLUS else Strand.MINUS
            location = Location(startOffset, startOffset + length, chromosome, strand)

            if (allowNs) {
                locations.add(location)
                break
            } else {
                // check that no unknown nucleotides in sequence, otherwise
                // re-generate location
                val sequence = chromosome.sequence
                val hasNs = (location.startOffset until location.endOffset).any { sequence.charAt(it) == Nucleotide.N }
                if (!hasNs) {
                    locations.add(location)
                    break
                }
            }

            if (attempts++ > SAMPLE_ATTEMPTS_THRESHOLD) {
                throw SamplingError(
                    "Failed to sample location length $length " +
                            "within [$leftBound, $rightBound) without N"
                )
            }
        }
    }

    return locations
}