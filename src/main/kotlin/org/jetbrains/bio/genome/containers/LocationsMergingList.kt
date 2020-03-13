package org.jetbrains.bio.genome.containers

import com.google.common.collect.Lists
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.toBedEntry
import org.jetbrains.bio.genome.format.unpackRegularFields
import java.io.IOException
import java.nio.file.Path
import java.util.*

/**
 * A location-based friend of [RangesMergingList].
 *
 * @author Oleg Shpynov
 * @author Sergei Lebedev
 */
class LocationsMergingList private constructor(
        private val rangeLists: GenomeStrandMap<RangesMergingList>)
    : GenomeStrandMapLike<List<Location>> {

    override val genomeQuery: GenomeQuery get() = rangeLists.genomeQuery

    /**
     * Performs element-wise union on locations in the two lists.
     */
    infix fun or(other: LocationsMergingList): LocationsMergingList =
            apply(other, RangesMergingList::or)

    /**
     * Performs element-wise intersection on locations in the two lists.
     */
    infix fun and(other: LocationsMergingList): LocationsMergingList =
            apply(other, RangesMergingList::and)

    fun apply(
            other: LocationsMergingList,
            op: (RangesMergingList, RangesMergingList) -> RangesMergingList
    ) = LocationsMergingList(rangeLists.merge(other.rangeLists, op))

    override operator fun get(chromosome: Chromosome, strand: Strand): List<Location> {
        val result = Lists.newArrayList<Location>()
        for ((startOffset, endOffset) in rangeLists[chromosome, strand]) {
            result.add(Location(startOffset, endOffset,
                    chromosome, strand))
        }
        return result
    }

    /**
     * XXX: this operation takes strand into account
     */
    operator fun contains(location: Location): Boolean =
            location.toRange() in rangeLists[location.chromosome, location.strand]


    fun contains(offset: Int, chr: Chromosome, strand: Strand): Boolean = offset in rangeLists[chr, strand]

    fun intersects(location: Location): Boolean = intersect(location).isNotEmpty()

    fun intersectsBothStrands(location: Location): Boolean = intersects(location) || intersects(location.opposite())

    fun intersectionLength(location: Location): Int = intersect(location).map(Range::length).sum()

    fun intersectionLengthBothStrands(location: Location): Int =
            intersectionLength(location) + intersectionLength(location.opposite())

    fun intersect(location: Location): List<Range> {
        val rangeList = rangeLists[location.chromosome, location.strand]
        return rangeList.intersect(location.toRange())
    }

    fun intersectBothStrands(location: Location): List<Range> =
            intersect(location) + intersect(location.opposite())

    @Throws(IOException::class)
    fun save(path: Path, format: BedFormat = BedFormat()) {
        format.print(path).use { printer ->
            locationIterator().forEach { printer.print(it.toBedEntry()) }
        }
    }

    private fun asLocationSequence(): Sequence<Location> = asSequence().flatMap { it.asSequence() }

    fun locationIterator(): Iterator<Location> = asLocationSequence().iterator()

    fun toList(): List<Location> = asLocationSequence().toList()

    /** The number of locations in this list. */
    val size: Int
        get() = genomeQuery.get().sumBy {
            (rangeLists[it, Strand.PLUS].size +
                    rangeLists[it, Strand.MINUS].size)
        }

    override fun toString() = "[${joinToString(", ")}]"

    fun filter(predicate: (Location) -> Boolean): LocationsMergingList {
        val builder = builder(genomeQuery)
        locationIterator().forEach { if (predicate(it)) builder.add(it) }
        return builder.build()
    }

    class Builder(private val genomeQuery: GenomeQuery) {
        private val ranges: GenomeStrandMap<ArrayList<Range>> =
                genomeStrandMap(genomeQuery) { _, _ -> arrayListOf<Range>() }

        fun add(location: Location): Builder {
            ranges[location.chromosome, location.strand].add(location.toRange())
            return this
        }

        fun build(): LocationsMergingList {
            val rangeLists = genomeStrandMap(genomeQuery) { chromosome, strand ->
                ranges[chromosome, strand].toRangeList()
            }

            return LocationsMergingList(rangeLists)
        }
    }

    companion object {
        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        fun create(genomeQuery: GenomeQuery, locations: Iterable<Location>) =
                create(genomeQuery, locations.iterator())

        fun create(genomeQuery: GenomeQuery, locations: Iterator<Location>): LocationsMergingList {
            val builder = Builder(genomeQuery)
            locations.forEach { builder.add(it) }
            return builder.build()
        }

        @Throws(IOException::class)
        fun load(genomeQuery: GenomeQuery, path: Path, format: BedFormat = BedFormat.auto(path)): LocationsMergingList {
            val builder = builder(genomeQuery)
            format.parse(path) {
                it.forEach {
                    val chromosome = genomeQuery[it.chrom]
                    if (chromosome != null) {
                        val e = it.unpackRegularFields(format)
                        builder.add(Location(e.start, e.end, chromosome, e.strand.toStrand()))
                    }
                }
            }
            return builder.build()
        }
    }
}

/**
 * Constructs a list of complementary locations.
 *
 * Input locations may be in any order, on any strand, and are allowed
 * to overlap. Input locations on a different chromosome or on an
 * opposite strand are ignored.
 *
 * @see Range.minus for details.
 */
operator fun Location.minus(locations: List<Location>): List<Location> {
    val ranges = locations.asSequence()
            .filter { it.strand == strand && it.chromosome == chromosome }
            .map { it.toRange() }
            .toList()

    return if (ranges.isEmpty()) {
        listOf(this)
    } else {
        (toRange() - ranges).map { it.on(chromosome, strand) }
    }
}

fun locationList(genomeQuery: GenomeQuery, vararg locations: Location): LocationsMergingList =
        locationList(genomeQuery, locations.asList())

fun locationList(genomeQuery: GenomeQuery, locations: Iterable<Location>): LocationsMergingList =
        LocationsMergingList.create(genomeQuery, locations)

private fun <T> GenomeStrandMap<T>.merge(
        other: GenomeStrandMap<T>, op: (T, T) -> T): GenomeStrandMap<T> {
    require(genomeQuery == other.genomeQuery)
    return genomeStrandMap(genomeQuery) { chromosome, strand ->
        op(get(chromosome, strand), other[chromosome, strand])
    }
}
