package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.util.bufferedReader
import java.io.IOException
import java.io.Reader
import java.nio.file.Path

/**
 * A location-based friend of [RangesMergingList].
 */
class LocationsMergingList private constructor(
    override val rangeLists: GenomeStrandMap<RangesMergingList>
) : LocationsList<RangesMergingList>() {

    /**
     * Performs element-wise union on locations in the two lists.
     */
    infix fun or(other: LocationsMergingList): LocationsMergingList =
        apply(other) { ra, rb ->
            (ra as RangesMergingList).or(rb as RangesMergingList)
        }

    override fun apply(
        other: LocationsList<out RangesList>,
        op: (RangesList, RangesList) -> Iterable<Range>
    ) = LocationsMergingList(genomeStrandMap(genomeQuery) { chromosome, strand ->
        val a = rangeLists[chromosome, strand]
        val b = other.rangeLists[chromosome, strand]
        op(a, b).toRangeMergingList()
    })

    fun apply(
        op: (RangesMergingList, Chromosome, Strand) -> RangesMergingList
    ) = LocationsMergingList(genomeStrandMap(genomeQuery) { chromosome, strand ->
        op(rangeLists[chromosome, strand], chromosome, strand)
    })


    fun intersects(location: Location): Boolean = intersect(location).isNotEmpty()

    fun intersectsBothStrands(location: Location): Boolean = intersects(location) || intersects(location.opposite())

    fun intersectionLength(location: Location): Int = intersect(location).map(Range::length).sum()

    fun intersectionLengthBothStrands(location: Location): Int =
        intersectionLength(location) + intersectionLength(location.opposite())

    fun intersect(location: Location): List<Range> {
        val rangeList = rangeLists[location.chromosome, location.strand]
        return rangeList.intersectRanges(location.toRange())
    }

    fun includes(location: Location): Boolean {
        val rangeList = rangeLists[location.chromosome, location.strand]
        return rangeList.includesRange(location.toRange())
    }

    fun intersectBothStrands(location: Location): List<Range> =
        intersect(location) + intersect(location.opposite())

    override fun toString() = "[${joinToString(", ")}]"

    fun filter(predicate: (Location) -> Boolean): LocationsMergingList {
        val builder = builder(genomeQuery)
        locationIterator().forEach { if (predicate(it)) builder.add(it) }
        return builder.build()
    }

    class Builder(gq: GenomeQuery) : LocationsListBuilder<LocationsMergingList>(gq) {
        override fun build() = LocationsMergingList(genomeStrandMap(genomeQuery) { chromosome, strand ->
            ranges[chromosome, strand].toRangeMergingList()
        })
    }

    companion object {
        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        fun create(genomeQuery: GenomeQuery, locations: Iterable<Location>) =
            create(genomeQuery, locations.iterator())

        fun create(genomeQuery: GenomeQuery, locations: Iterator<Location>) =
            Builder(genomeQuery).also { builder ->
                locations.forEach { builder.add(it) }
            }.build()

        @Throws(IOException::class)
        fun load(
            gq: GenomeQuery,
            path: Path,
            format: BedFormat = BedFormat.auto(path),
            entry2LocationFun: (Chromosome, BedEntry, BedFormat) -> Location = EXTBED_2_LOC_FUN
        ) = load(
            gq, path.bufferedReader(), "${path.toAbsolutePath()}", format, entry2LocationFun
        )

        @Throws(IOException::class)
        fun load(
            gq: GenomeQuery,
            reader: Reader,
            src: String,
            format: BedFormat,
            entry2LocationFun: (Chromosome, BedEntry, BedFormat) -> Location = EXTBED_2_LOC_FUN
        ) = LocationsList.load(builder(gq), reader, src, format, entry2LocationFun)
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

private fun <T> GenomeStrandMap<T>.apply(
    other: GenomeStrandMap<T>, op: (T, T) -> T
): GenomeStrandMap<T> {
    require(genomeQuery == other.genomeQuery)
    return genomeStrandMap(genomeQuery) { chromosome, strand ->
        op(get(chromosome, strand), other[chromosome, strand])
    }
}
