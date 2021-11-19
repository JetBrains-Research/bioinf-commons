package org.jetbrains.bio.genome

import com.google.common.collect.ComparisonChain
import com.google.common.math.IntMath
import com.google.common.primitives.Ints
import com.google.gson.TypeAdapter
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import java.math.RoundingMode
import java.util.stream.IntStream
import java.util.stream.Stream
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

/**
 * A semi-closed interval.
 *
 * @author Oleg Shpynov
 */
data class Range(
    /** 0-based start offset (inclusive). */
    val startOffset: Int,
    /** 0-based end offset (exclusive). */
    val endOffset: Int
) : Comparable<Range> {

    init {
        require(startOffset <= endOffset) { "invalid range $this" }
    }

    fun length() = endOffset - startOffset

    fun isEmpty() = length() == 0
    fun isNotEmpty() = length() != 0

    infix fun intersects(other: Range): Boolean {
        return other.endOffset > startOffset && endOffset > other.startOffset
    }

    infix fun intersection(other: Range): Range {
        return if (this intersects other) {
            Range(
                max(startOffset, other.startOffset),
                min(endOffset, other.endOffset)
            )
        } else {
            EMPTY
        }
    }

    infix fun union(other: Range): Range {
        return Range(
            min(startOffset, other.startOffset),
            max(endOffset, other.endOffset)
        )
    }

    infix fun distanceTo(other: Range): Int {
        return if (intersects(other))
            0
        else
            min(
                abs(startOffset - other.endOffset),
                abs(endOffset - other.startOffset)
            )
    }

    fun on(chromosome: Chromosome) = ChromosomeRange(startOffset, endOffset, chromosome)

    fun on(chromosome: Chromosome, strand: Strand) = Location(startOffset, endOffset, chromosome, strand)

    operator fun contains(offset: Int) = offset in startOffset..(endOffset - 1)

    operator fun contains(other: Range) = startOffset <= other.startOffset && other.endOffset <= endOffset

    /**
     * Returns an ordered stream of sub-ranges each having a given
     * `width`, with the last sub-range possibly being an exception.
     */
    fun slice(width: Int): Stream<Range> {
        if (isEmpty()) {
            return Stream.empty()
        }

        val n = IntMath.divide(length(), width, RoundingMode.CEILING)
        return IntStream.range(0, n).mapToObj { i ->
            Range(
                startOffset + i * width,
                min(endOffset, startOffset + (i + 1) * width)
            )
        }
    }

    override fun toString() = "[$startOffset, $endOffset)"

    override fun compareTo(other: Range) = ComparisonChain.start()
        .compare(startOffset, other.startOffset)
        .compare(endOffset, other.endOffset)
        .result()

    companion object {
        /** An empty range. */
        val EMPTY = Range(0, 0)

        /**
         * A type adapter enforcing a more compact encoding for ranges.
         *
         * The notation `[24, 42]` can be confusing, because ranges are
         * semi-closed on the right. If you don't like it, feel free to
         * alter the right bound on serialization/deserialization.
         */
        internal val ADAPTER = object : TypeAdapter<Range>() {
            override fun read(`in`: JsonReader) = with(`in`) {
                beginArray()
                val startOffset = nextInt()
                val endOffset = nextInt()
                endArray()
                Range(startOffset, endOffset)
            }

            override fun write(out: JsonWriter, range: Range) {
                out.beginArray()
                    .value(range.startOffset)
                    .value(range.endOffset)
                    .endArray()
            }
        }.nullSafe()
    }
}

data class ChromosomeRange(
    val startOffset: Int,
    val endOffset: Int,
    val chromosome: Chromosome
) : LocationAware, Comparable<ChromosomeRange> {

    init {
        require(startOffset <= endOffset) { "invalid chromosome range $this" }
    }

    override val location: Location get() = on(Strand.PLUS)

    fun length() = endOffset - startOffset

    fun on(strand: Strand) = Location(startOffset, endOffset, chromosome, strand)

    fun toRange() = Range(startOffset, endOffset)

    override fun compareTo(other: ChromosomeRange) = ComparisonChain.start()
        // sort by chr name, start, end offset
        .compare(chromosome.name, other.chromosome.name)
        .compare(startOffset, other.startOffset)
        .compare(endOffset, other.endOffset)
        .result()

    override fun toString() = "${chromosome.name}:[$startOffset, $endOffset)"
}

data class Location(
    val startOffset: Int, val endOffset: Int,
    val chromosome: Chromosome,
    val strand: Strand = Strand.PLUS
) :
    LocationAware, Comparable<Location> {

    init {
        require(startOffset <= endOffset) { "invalid location $this" }
    }

    @Deprecated("", ReplaceWith("this"))
    override val location: Location
        get() = this

    val sequence: String
        get() {
            return chromosome.sequence.substring(startOffset, endOffset, strand)
        }

    operator fun contains(offset: Int) = offset in startOffset..(endOffset - 1)

    fun length() = endOffset - startOffset

    fun opposite() = Location(startOffset, endOffset, chromosome, strand.opposite())


    /**
     * Returns absolute position of offset relative to 5' bound. Differs
     * from `startOffset + relativeOffset` in case of [Strand.MINUS]
     */
    fun get5Bound(relativeOffset: Int = 0) = Companion.get5Bound(startOffset, endOffset, strand, relativeOffset)

    /**
     * Returns absolute position of offset relative to 3' bound. Differs
     * from `getEndOffset + relativeOffset` in case of [Strand.MINUS].
     */
    fun get3Bound(relativeOffset: Int = 0) = Companion.get3Bound(startOffset, endOffset, strand, relativeOffset)

    fun toRange() = Range(startOffset, endOffset)

    fun toChromosomeRange() = ChromosomeRange(startOffset, endOffset, chromosome)

    override fun toString() = "${chromosome.name}:$strand[$startOffset, $endOffset)"

    override fun compareTo(other: Location) = ComparisonChain.start()
        // sort locs by chr name, strand, start, end offset
        .compare(chromosome.name, other.chromosome.name)
        .compare(strand, other.strand)
        .compare(startOffset, other.startOffset)
        .compare(endOffset, other.endOffset)
        .result()

    companion object {
        internal val ADAPTER = object : TypeAdapter<Location>() {
            override fun read(`in`: JsonReader) = with(`in`) {
                beginArray()
                val startOffset = nextInt()
                val endOffset = nextInt()
                val chromosome = Chromosome.ADAPTER.read(this)
                val strand = nextString().toStrand()
                endArray()
                Location(startOffset, endOffset, chromosome, strand)
            }

            override fun write(out: JsonWriter, location: Location) {
                out.beginArray()
                    .value(location.startOffset)
                    .value(location.endOffset)

                Chromosome.ADAPTER.write(out, location.chromosome)

                out.value(location.strand.char.toString())
                    .endArray()
            }
        }.nullSafe()

        fun intersects(l1: Location?, l2: Location?): Boolean {
            if (l1 == null || l2 == null) {
                return false
            }
            if (l1.chromosome != l2.chromosome) {
                return false
            }
            return l1.toRange().intersects(l2.toRange())
        }

        /**
         * Returns absolute position of offset relative to 5' bound. Differs
         * from `startOffset + relativeOffset` in case of [Strand.MINUS]
         */
        fun get5Bound(
            startOffset: Int, endOffset: Int, strand: Strand,
            relativeOffset: Int = 0
        ): Int {
            return if (strand.isPlus()) {
                startOffset + relativeOffset
            } else {
                endOffset - 1 - relativeOffset
            }
        }

        /**
         * Returns absolute position of offset relative to 3' bound. Differs
         * from `getEndOffset + relativeOffset` in case of [Strand.MINUS].
         */
        fun get3Bound(
            startOffset: Int, endOffset: Int, strand: Strand,
            relativeOffset: Int = 0
        ): Int {
            return if (strand.isPlus()) {
                endOffset - 1 + relativeOffset
            } else {
                startOffset - relativeOffset
            }
        }

        fun valueOf(str: String, genome: Genome): Location? {
            return try {
                val split = str.split(":", "[", ", ", ")")
                val chr = Chromosome(genome, split[0])
                val strand = split[1].toStrand()
                val start = split[2].toInt()
                val end = split[3].toInt()
                Location(start, end, chr, strand)
            } catch (e: Exception) {
                null
            }
        }
    }
}

enum class RelativePosition {
    FIVE_PRIME {
        override fun of(location: Location, relativeStartOffset: Int, relativeEndOffset: Int): Location {
            return bracket(
                location, location.get5Bound(relativeStartOffset),
                location.get5Bound(relativeEndOffset - 1)
            )
        }
    },

    THREE_PRIME {
        override fun of(location: Location, relativeStartOffset: Int, relativeEndOffset: Int): Location {
            return bracket(
                location, location.get3Bound(relativeStartOffset),
                location.get3Bound(relativeEndOffset - 1)
            )
        }
    },

    ALL {
        override fun of(location: Location, relativeStartOffset: Int, relativeEndOffset: Int): Location {
            return if (relativeStartOffset == 0 && relativeEndOffset == 0) {
                location  // Don't allocate new location in this case
            } else {
                bracket(
                    location, location.get5Bound(relativeStartOffset),
                    location.get3Bound(relativeEndOffset - 1)
                )
            }
        }
    };

    /**
     * Computes location according given relative position and offsets.
     *
     * For more details see: [Location#get5Bound] and [Location#get3Bound]
     */
    abstract fun of(location: Location, relativeStartOffset: Int, relativeEndOffset: Int): Location

    /**
     * Ensures the resulting location is proper, i.e. has correctly ordered
     * endpoints which are bounded by [0, chromosome.length].
     */
    protected fun bracket(location: Location, newStartOffset: Int, newEndOffset: Int): Location {
        return with(location) {
            val boundedStartOffset = min(
                max(0, min(newStartOffset, newEndOffset)),
                chromosome.length
            )
            // XXX do we even need to use boundedStartOffset here?
            val boundedEndOffset = min(
                Ints.max(boundedStartOffset, newStartOffset, newEndOffset) + 1,
                chromosome.length
            )
            Location(boundedStartOffset, boundedEndOffset, chromosome, strand)
        }
    }
}

/** A simplistic interface for a "thing" with genomic location. */
interface LocationAware {
    val location: Location
}
