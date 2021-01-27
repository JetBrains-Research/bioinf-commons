package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.util.bufferedReader
import java.io.IOException
import java.io.Reader
import java.nio.file.Path

/**
 * A location-based friend of [RangesSortedList].
 */
class LocationsSortedList private constructor(
    override val rangeLists: GenomeStrandMap<RangesSortedList>
) : LocationsList<RangesSortedList>() {

    override fun apply(
        other: LocationsList<out RangesList>,
        op: (RangesList, RangesList) -> Iterable<Range>
    ) = LocationsSortedList(genomeStrandMap(genomeQuery) { chromosome, strand ->
        val a = rangeLists[chromosome, strand]
        val b = other.rangeLists[chromosome, strand]
        op(a, b).toRangeSortedList()
    })
    
    class Builder(gq: GenomeQuery): LocationsListBuilder<LocationsSortedList>(gq) {
        override fun build() = LocationsSortedList(genomeStrandMap(genomeQuery) { chromosome, strand ->
            ranges[chromosome, strand].toRangeSortedList()
        })
    }

    override fun toString() = "[${joinToString(", ")}]"

    companion object {
        fun builder(genomeQuery: GenomeQuery) = Builder(genomeQuery)

        fun create(rangeLists: GenomeStrandMap<RangesSortedList>) = LocationsSortedList(rangeLists)

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
            entry2LocationFun: (Chromosome, BedEntry, BedFormat) -> Location = LocationsList.EXTBED_2_LOC_FUN
        ) = load(
            gq, path.bufferedReader(), "${path.toAbsolutePath()}", format, entry2LocationFun
        )

        @Throws(IOException::class)
        fun load(
            gq: GenomeQuery,
            reader: Reader,
            src: String,
            format: BedFormat,
            entry2LocationFun: (Chromosome, BedEntry, BedFormat) -> Location = LocationsList.EXTBED_2_LOC_FUN
        ) = LocationsList.load(builder(gq), reader, src, format, entry2LocationFun)
    }
}