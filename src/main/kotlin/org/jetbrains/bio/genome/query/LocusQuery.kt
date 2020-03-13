package org.jetbrains.bio.genome.query

import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.locationList
import org.jetbrains.bio.genome.containers.minus
import org.jetbrains.bio.util.Lexeme
import org.jetbrains.bio.util.RegexLexeme
import org.jetbrains.bio.util.Tokenizer
import org.jetbrains.bio.util.parseInt
import java.util.*


abstract class TranscriptLocusQuery(override val id: String)
    : Query<Transcript, Collection<Location>> {

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is TranscriptLocusQuery) return false

        if (id != other.id) return false

        return true
    }

    override fun hashCode(): Int {
        return id.hashCode()
    }

    fun toGeneQuery(): Query<Gene, Collection<Location>> {
        return object : Query<Gene, Collection<Location>> {
            override val id: String
                get() = this@TranscriptLocusQuery.id

            override fun apply(t: Gene): Collection<Location> {
                // Take the longest transcript
                return this@TranscriptLocusQuery.apply(
                        t.transcripts.maxBy { it.length() }!!)
            }

        }
    }

    override fun toString(): String = id

    companion object {
        /**
         * Locus for given string
         * TSS[500]
         * TSS[-2000..2000]
         */
        fun parse(text: String): TranscriptLocusQuery? {
            val lcText = text.toLowerCase()
            if (lcText in TRANSCRIPT_LOCUS_QUERIES) {
                return TRANSCRIPT_LOCUS_QUERIES[lcText]!!
            }

            val delimiter = RegexLexeme("[_;]|\\.\\.")
            val lBracket = Lexeme("[")
            val rBracket = Lexeme("]")

            // Only TSS and TES can have params
            if (!text.startsWith("tss") && !text.startsWith("tes")) {
                return null
            }

            val tokenizer = Tokenizer(lcText.substring(3), setOf(delimiter, lBracket, rBracket))
            val bracketsFound = tokenizer.fetch() == lBracket
            if (bracketsFound) {
                tokenizer.next()
            }
            return try {
                val value = tokenizer.parseInt()
                if (tokenizer.fetch() == delimiter) {
                    tokenizer.next()
                    val end = tokenizer.parseInt()
                    when (lcText.subSequence(0, 3)) {
                        "tss" -> TssQuery(value, end)
                        "tes" -> TesQuery(value, end)
                        else -> error("Unsupported query: $text")
                    }
                } else {
                    when (lcText.subSequence(0, 3)) {
                        "tss" -> TssQuery(-Math.abs(value), Math.abs(value))
                        "tes" -> TesQuery(-Math.abs(value), Math.abs(value))
                        else -> error("Unsupported query: $text")
                    }
                }
            } finally {
                if (bracketsFound) {
                    tokenizer.check(rBracket)
                }
                tokenizer.checkEnd()
            }
        }
    }
}

val TRANSCRIPT_LOCUS_QUERIES = mapOf(
        "tss" to TssQuery(),
        "tes" to TesQuery(),
        "transcript" to TranscriptQuery(),
        "cds" to CDSQuery(),
        "utr5" to UTR5Query(),
        "utr3" to UTR3Query(),
        "introns" to IntronsQuery(),
        "exons" to ExonsQuery())

class CDSQuery : TranscriptLocusQuery("cds") {
    override fun apply(input: Transcript) = input.cds
}

class ExonsQuery : TranscriptLocusQuery("exons") {
    override fun apply(input: Transcript) = input.exons
}

class TranscriptQuery : TranscriptLocusQuery("transcript") {
    override fun apply(input: Transcript) = listOf(input.location)
}

class IntronsQuery : TranscriptLocusQuery("introns") {
    override fun apply(input: Transcript) = input.introns
}

/** Transcription end site query. */
class TesQuery @JvmOverloads constructor(
        private val leftBound: Int = -2000,
        private val rightBound: Int = 2000) : TranscriptLocusQuery("tes") {

    init {
        check(leftBound < rightBound)
    }

    override fun apply(input: Transcript): Collection<Location> {
        return listOf(RelativePosition.THREE_PRIME.of(input.location, leftBound, rightBound))
    }

    override val id: String get() = "${super.id}[$leftBound..$rightBound]"

    override fun equals(other: Any?) = when {
        this === other -> true
        other == null || other !is TesQuery -> false
        else -> leftBound == other.leftBound && rightBound == other.rightBound
    }

    override fun hashCode() = Objects.hash(super.hashCode(), leftBound, rightBound)
}

/** Transcription start site query. */
class TssQuery @JvmOverloads constructor(
        private val leftBound: Int = -2000,
        private val rightBound: Int = 2000) : TranscriptLocusQuery("tss") {

    init {
        check(leftBound < rightBound)
    }

    override fun apply(input: Transcript): Collection<Location> {
        return listOf(RelativePosition.FIVE_PRIME.of(input.location, leftBound, rightBound))
    }

    override val id: String get() = "${super.id}[$leftBound..$rightBound]"

    override fun equals(other: Any?) = when {
        this === other -> true
        other == null || other !is TssQuery -> false
        else -> leftBound == other.leftBound && rightBound == other.rightBound
    }

    override fun hashCode() = Objects.hash(super.hashCode(), leftBound, rightBound)
}

class UTR3Query : TranscriptLocusQuery("utr3") {
    override fun apply(input: Transcript) = input.utr3
}

class UTR5Query : TranscriptLocusQuery("utr5") {
    override fun apply(input: Transcript) = input.utr5
}

class NearbyTranscriptQuery : TranscriptLocusQuery("nearby_transcript") {
    override fun apply(input: Transcript): List<Location> {
        return listOf(RelativePosition.ALL.of(input.location, -2000, 2000))
    }
}

abstract class ChromosomeLocusQuery(override val id: String) : Query<Chromosome, Collection<Location>> {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is ChromosomeLocusQuery) return false
        if (id != other.id) return false
        return true
    }

    override fun hashCode(): Int {
        return id.hashCode()
    }

    override fun toString(): String = id
}

val CHROMOSOME_LOCUS_QUERIES = mapOf(
        "repeats" to RepeatsQuery(),
        "non_repeats" to NonRepeatsQuery(),
        "cpg_islands" to CGIQuery()
)

class CGIQuery : ChromosomeLocusQuery("cgi") {
    override fun apply(input: Chromosome): Collection<Location> {
        return input.cpgIslands.map { it.location }
    }
}

class RepeatsQuery(private val repeatClass: String? = null) : ChromosomeLocusQuery("repeats") {
    override fun apply(input: Chromosome): Collection<Location> {
        return input.repeats.asSequence().filter { repeat ->
            repeatClass == null || repeat.repeatClass == repeatClass
        }.map { it.location }.toList()
    }

    override val id: String
        get() {
            return super.id + if (repeatClass == null) "" else "[$repeatClass]"
        }
}

class NonRepeatsQuery : ChromosomeLocusQuery("non_repeats") {
    override fun apply(input: Chromosome): Collection<Location> {
        val genomeQuery = GenomeQuery(input.genome, input.name)
        val repeats = locationList(genomeQuery, input.repeats.map { it.location })

        return Strand.values().flatMap {
            val location = input.range.on(input, it)
            location - repeats[input, Strand.PLUS]
        }
    }
}