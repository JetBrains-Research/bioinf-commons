@file:Suppress("unused")

package org.jetbrains.bio.genome

import com.google.common.collect.ImmutableListMultimap
import com.google.common.collect.ListMultimap
import org.jetbrains.bio.big.BigBedFile
import java.io.BufferedReader
import java.io.FileInputStream
import java.io.InputStreamReader
import java.nio.file.Path
import java.util.zip.GZIPInputStream


object RepeatsHs1 {

    internal fun readHs1Repeats(genome: Genome, repeatsPath: Path): ListMultimap<Chromosome, Repeat> {
        val chromosomes = genome.chromosomeNamesMap

        val reader = InputStreamReader(GZIPInputStream(FileInputStream(repeatsPath.toFile())))

        val builder = ImmutableListMultimap.builder<Chromosome, Repeat>()
        BufferedReader(reader).lines()
            .forEach {
                val parsingResult = parseHs1RepeatsLine(it, chromosomes)
                if (parsingResult != null) builder.put(parsingResult.first, parsingResult.second)
            }

        return builder.build()
    }

    private fun isValidHs1RepeatsLine(line: String): Boolean {
        val substringPresentOnlyInHeader = " repeat "
        val chrMarker = "chr"

        return line.isNotBlank() && line.contains(chrMarker) && !line.contains(substringPresentOnlyInHeader)
    }

    internal fun parseHs1RepeatsLine(line: String, chromosomes: Map<String, Chromosome>): Pair<Chromosome, Repeat>? {
        if (!isValidHs1RepeatsLine(line)) {
            return null
        }

        try {
            val list = line.trim().split("\\s+".toRegex())
            if (list.size != 15) {
                return null
            }

            val chromosome = chromosomes[list[4]] ?: return null
            val strand = list[8].toStrand()
            val startOffset = list[5].toInt()
            val endOffset = list[6].toInt()
            val location = Location(startOffset, endOffset, chromosome, strand)

            val classFamily = list[10]
            val repeatClass: String
            val repeatFamily: String
            if (classFamily.contains("/")) {
                val substrings = classFamily.split("/")
                repeatClass = substrings[0].lowercase()
                repeatFamily = substrings[1].lowercase()
            } else {
                repeatClass = classFamily
                repeatFamily = classFamily
            }

            val repeat = Repeat(
                list[9], location,
                repeatClass,
                repeatFamily
            )

            return Pair(chromosome, repeat)
        } catch (e: Exception) {
            return null
        }
    }

    fun String.toStrand() = single().toStrand()

    fun Char.toStrand() = when (this) {
        '+' -> Strand.PLUS
        'C' -> Strand.MINUS
        else -> throw IllegalStateException("Unknown strand: $this")
    }
}

object CpGIslandsHs1 {
    internal fun readHs1CpGIslands(
        genome: Genome,
        islandsPath: Path,
        chrAltName2CanonicalMapping: Map<String, String>,
    ): ListMultimap<Chromosome, CpGIsland> {
        val chromosomesNames = genome.chromosomeNamesMap

        val builder = ImmutableListMultimap.builder<Chromosome, CpGIsland>()

        BigBedFile.read(islandsPath).use { bbf ->
            for (chromosome in bbf.chromosomes.valueCollection()) {
                for ((chrom, start, end, rest) in bbf.query(chromosome)) {
                    val chr = chromosomesNames[chrAltName2CanonicalMapping[chrom]]
                    val island = parseIsland(chr, start, end, rest)

                    if (chr != null && island != null) {
                        builder.put(chr, island)
                    }
                }
            }
        }

        return builder.build()
    }

    internal fun parseIsland(chr: Chromosome?, start: Int, end: Int, rest: String): CpGIsland? {
        if (chr == null) {
            return null
        }

        val otherFields = rest.split("\t")
        if (otherFields.size < 7) {
            return null
        }
        val location = Location(
            start, end,
            chr
        )

        return CpGIsland(
            otherFields[2].toInt(), otherFields[3].toInt(),
            otherFields[6].toDouble(), location
        )
    }
}