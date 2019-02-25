package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.big.ExtendedBedEntry
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.sequence.Nucleotide
import org.jetbrains.bio.io.BedFormat
import java.nio.file.Files
import java.nio.file.Path
import java.util.*
import java.util.concurrent.atomic.AtomicLong

class SamplingError(message: String) : Exception(message)

private val RANDOM = XSRandom()
private val SAMPLE_ATTEMPTS_THRESHOLD = 1000

/**
 * Shuffles bed regions positions inside chromosomes. All bed regions stays in same chromosome where they were
 * and their length remains unchanged.
 *
 * @return path to newly generated Bed file.
 */
fun randomizeBedRegions(bedFilePath: Path, genomeQuery: GenomeQuery): Path {
    val randomBedPath = Files.createTempFile("randomRegions", ".bed")
    val printer = BedFormat().print(randomBedPath)

    BedFormat.auto(bedFilePath).parse(bedFilePath) {
        it.distinct().filter { it.chrom in genomeQuery }
                .map { it.chrom to (it.end - it.start) }
                .groupBy { it.first }
                .mapValues { it.value.map { p -> p.second } }
                .toList().parallelStream()
                .map {
                    val chr = Chromosome(genomeQuery.build, it.first)
                    sampleLocations(chr, it.second, 0, chr.length, false, Strand.PLUS)
                }
                .sequential() // required here, otherwise flatMap().forEach() seems to look concurrently
                .flatMap { it.stream() }
                .forEach { (start, end, chr, strand) ->
                    printer.print(ExtendedBedEntry(chr.name, start, end, strand = strand.char))
                }
    }
    printer.close()
    return randomBedPath
}

/**
 * Generates list of given size containing random locations (which may intersect)
 * of given length within chromosome[leftBound, rightBound)
 */
fun sampleLocations(chromosome: Chromosome,
                    lengthList: List<Int>,
                    leftBound: Int,
                    rightBound: Int,
                    allowNs: Boolean = false,
                    predefinedStrand: Strand? = null): List<Location> {
    if (rightBound - leftBound <= lengthList.max() ?: 0) {
        throw SamplingError("Not enough space for location length ${lengthList.max()} within [$leftBound, $rightBound)")
    }

    val locations = arrayListOf<Location>()
    val sequence = chromosome.sequence
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
                val hasNs = (location.startOffset until location.endOffset).any { sequence.charAt(it) == Nucleotide.N }
                if (!hasNs) {
                    locations.add(location)
                    break
                }
            }

            if (attempts++ > SAMPLE_ATTEMPTS_THRESHOLD) {
                throw SamplingError("Failed to sample location length $length " +
                        "within [$leftBound, $rightBound) without N")
            }
        }
    }

    return locations
}

/**
 * See {@see http://en.wikipedia.org/wiki/Xorshift}
 * This implementation is about 30% faster than the generator from java.util.random.
 * It's output passes the Dieharder test suite with no fail and only two announced weaknesses.
 * Original code: {@see http://demesos.blogspot.ru/2011/09/replacing-java-random-generator.html}

 * @author Oleg Shpynov
 */
class XSRandom(private var seed_: Long = XSRandom.seedUniquifier() xor System.nanoTime()) : Random() {

    override fun next(nbits: Int): Int {
        var x = seed_
        x = x xor (x shl 21)
        x = x xor x.ushr(35)
        x = x xor (x shl 4)
        seed_ = x
        x = x and (1L shl nbits) - 1
        return x.toInt()
    }

    companion object {
        /**
         * The latter 2 declarations are taken from class [java.util.Random]
         * because of visibility reasons
         */
        private fun seedUniquifier(): Long {
            // L'Ecuyer, "Tables of Linear Congruential Generators of Different Sizes and Good Lattice Structure", 1999
            while (true) {
                val current = seedUniquifier.get()
                val next = current * 181783497276652981L
                if (seedUniquifier.compareAndSet(current, next))
                    return next
            }
        }

        private val seedUniquifier = AtomicLong(8682522807148012L)
    }
}