@file:Suppress("unused")

package org.jetbrains.bio.genome

import com.google.common.collect.ImmutableListMultimap
import com.google.common.collect.ListMultimap
import java.io.BufferedReader
import java.io.FileInputStream
import java.io.InputStreamReader
import java.nio.file.Path
import java.util.zip.GZIPInputStream


object RepeatsHs1 {

    internal fun readHs1(genome: Genome, repeatsPath: Path): ListMultimap<Chromosome, Repeat> {
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