package org.jetbrains.bio.genome.containers

import com.google.common.collect.Lists
import kotlinx.support.jdk7.use
import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.toBedEntry
import org.jetbrains.bio.genome.format.unpackRegularFields
import org.jetbrains.bio.viktor.KahanSum
import java.io.IOException
import java.io.Reader
import java.nio.file.Path
import java.util.*

abstract class LocationsList<T : RangesList> : GenomeStrandMapLike<List<Location>> {
    abstract val rangeLists: GenomeStrandMap<T>

    override val genomeQuery: GenomeQuery get() = rangeLists.genomeQuery

    /** The number of locations in this list. */
    val size: Int
        get() = genomeQuery.get().sumBy {
            (rangeLists[it, Strand.PLUS].size +
                    rangeLists[it, Strand.MINUS].size)
        }

    abstract fun apply(
        other: LocationsList<out RangesList>,
        op: (RangesList, RangesList) -> Iterable<Range>
    ): LocationsList<T>

    /**
     * Performs element-wise intersection on locations in the two lists.
     */
    infix fun intersectRanges(other: LocationsList<*>): LocationsList<*> =
        apply(other) { ra, rb -> ra.intersectRanges(rb) }

    /**
     * 'inline' is important here otherwise each method usage will
     * create new instance of anonymous 'metric' class
     */
    inline fun calcAdditiveMetric(
        other: LocationsList<out RangesList>,
        metric: (RangesList, RangesList) -> Long
    ): Long {
        require(genomeQuery == other.genomeQuery)

        var acc = 0L
        genomeQuery.get().forEach { chromosome ->
            Strand.values().forEach { strand ->
                acc += metric(
                    rangeLists[chromosome, strand],
                    other.rangeLists[chromosome, strand]
                )
            }
        }
        return acc
    }

    /**
     * 'inline' is important here otherwise each method usage will
     * create new instance of anonymous 'metric' class
     */
    inline fun calcAdditiveMetricDouble(
        other: LocationsList<out RangesList>,
        metric: (RangesList, RangesList) -> Double
    ): Double {
        require(genomeQuery == other.genomeQuery)

        val acc = KahanSum()
        genomeQuery.get().forEach { chromosome ->
            Strand.values().forEach { strand ->
                acc += metric(
                    rangeLists[chromosome, strand],
                    other.rangeLists[chromosome, strand]
                )
            }
        }

        return acc.result()
    }

    fun asLocationSequence(): Sequence<Location> = asSequence().flatMap { it.asSequence() }

    fun locationIterator(): Iterator<Location> = asLocationSequence().iterator()

    fun toList(): List<Location> = asLocationSequence().toList()

    override operator fun get(chromosome: Chromosome, strand: Strand): List<Location> {
        val result = Lists.newArrayList<Location>()
        for ((startOffset, endOffset) in rangeLists[chromosome, strand]) {
            result.add(
                Location(startOffset, endOffset, chromosome, strand)
            )
        }
        return result
    }

    /**
     * XXX: this operation takes strand into account
     */
    operator fun contains(location: Location): Boolean =
        rangeLists[location.chromosome, location.strand].includesRange(location.toRange())


    fun contains(offset: Int, chr: Chromosome, strand: Strand): Boolean = rangeLists[chr, strand].includesRange(offset)

    @Throws(IOException::class)
    fun save(path: Path, format: BedFormat = BedFormat()) {
        format.print(path).use { printer ->
            locationIterator().forEach { printer.print(it.toBedEntry()) }
        }
    }

    companion object {
        val EXTBED_2_LOC_FUN: (Chromosome, BedEntry, BedFormat) -> Location = { chr, e, fmt ->
            val ex = e.unpackRegularFields(fmt)
            Location(ex.start, ex.end, chr, ex.strand.toStrand())
        }

        @Throws(IOException::class)
        fun <T> load(
            builder: LocationsListBuilder<T>,
            reader: Reader,
            src: String,
            format: BedFormat,
            entry2LocationFun: (Chromosome, BedEntry, BedFormat) -> Location
        ): T {
            val gq = builder.genomeQuery
            format.parse(reader, src) { parser ->
                parser.forEach { entry ->
                    val chromosome = gq[entry.chrom]
                    if (chromosome != null) {
                        builder.add(entry2LocationFun(chromosome, entry, format))
                    }
                }
            }
            return builder.build()
        }
    }
}

abstract class LocationsListBuilder<T>(val genomeQuery: GenomeQuery) {
    abstract fun build(): T

    protected val ranges: GenomeStrandMap<ArrayList<Range>> =
        genomeStrandMap(genomeQuery) { _, _ -> arrayListOf<Range>() }

    fun add(location: Location): LocationsListBuilder<T> {
        ranges[location.chromosome, location.strand].add(location.toRange())
        return this
    }

}
