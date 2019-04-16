package org.jetbrains.bio.methylome

import com.google.common.base.MoreObjects
import com.google.common.primitives.Shorts
import gnu.trove.map.hash.TObjectIntHashMap
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.npy.NpzFile
import java.io.IOException
import java.nio.file.Path
import java.util.*
import java.util.concurrent.atomic.AtomicLong

/**
 * A container for WGBS data.
 *
 * Contracts:
 *
 *   * per-strand and per-chromosome counts are sorted by strand,
 *   * all cytosines are covered by at least a single read.
 *
 * @author Sergei Lebedev
 */
open class Methylome internal constructor(
        val genomeQuery: GenomeQuery,
        private val plusData: GenomeMap<MethylomeFrame>,
        private val minusData: GenomeMap<MethylomeFrame>) {

    operator fun get(chromosome: Chromosome, strand: Strand): StrandMethylomeView {
        return StrandMethylomeView(getInternal(chromosome, strand))
    }

    fun getCombined(chromosome: Chromosome): ChromosomeMethylomeView {
        return ChromosomeMethylomeView(getInternal(chromosome, Strand.PLUS),
                                       getInternal(chromosome, Strand.MINUS))
    }

    /**
     * Returns the frame corresponding to a given [chromosome] and [strand].
     */
    internal open fun getInternal(chromosome: Chromosome, strand: Strand): MethylomeFrame {
        return strand.choose(plusData, minusData)[chromosome]
    }

    val size: Int get() = genomeQuery.get().sumBy { getCombined(it).size }

    @Throws(IOException::class)
    fun save(outputPath: Path) {
        check(genomeQuery.restriction == null)
        NpzFile.write(outputPath).use { writer ->
            writer.write("version", intArrayOf(VERSION))

            val chromosomes = genomeQuery.get()
            for (chromosome in chromosomes) {
                for (strand in Strand.values()) {
                    val frame = getInternal(chromosome, strand)
                    frame.save((chromosome to strand).toKey(), writer)
                }
            }
        }
    }

    override fun toString() = MoreObjects.toStringHelper(this)
            .addValue(genomeQuery).toString()

    companion object {
        /** Binary format version.  */
        internal val VERSION: Int = 6

        fun builder(genomeQuery: GenomeQuery, ignoreDuplicatedOffsets: Boolean = false): MethylomeBuilder = MethylomeBuilder(
                genomeQuery, ignoreDuplicatedOffsets
        )

        fun lazy(genomeQuery: GenomeQuery, inputPath: Path): Methylome {
            return LazyMethylome(genomeQuery, inputPath)
        }
    }
}

/**
 * A methylome which loads its data on demand.
 *
 * The current implementation *never* invalidates the loaded entries.
 */
private class LazyMethylome(genomeQuery: GenomeQuery,
                            private val inputPath: Path)
    : Methylome(genomeQuery, frameMap(genomeQuery), frameMap(genomeQuery)) {

    private val chromosomeNamesMap = TObjectIntHashMap<String>().apply {
        genomeQuery.get().sortedBy { it.name }.forEachIndexed { i, n ->
            this.put(n.name, i)
        }
    }

    private val plusMask: BitSet = BitSet(chromosomeNamesMap.size())
    private val minusMask: BitSet = BitSet(chromosomeNamesMap.size())

    override fun getInternal(chromosome: Chromosome, strand: Strand): MethylomeFrame {
        val mask = strand.choose(plusMask, minusMask)
        val frame = super.getInternal(chromosome, strand)
        val key = chromosomeNamesMap[chromosome.name]
        if (!mask[key]) {
            NpzFile.read(inputPath).use { reader ->
                val version = reader["version"].asIntArray().single()
                require(version == Methylome.VERSION) {
                    "methylome version is $version instead of ${Methylome.VERSION}"
                }

                frame.load((chromosome to strand).toKey(), reader)
            }
            mask.set(key)
        }

        return frame
    }
}

private fun frameMap(genomeQuery: GenomeQuery): GenomeMap<MethylomeFrame> {
    return genomeMap(genomeQuery) { MethylomeFrame() }
}

private fun Pair<Chromosome, Strand>.toKey(): String {
    return first.name + '/' + second
}

/**
 * A builder for [Methylome].
 *
 * Guarantees that the constructed [Methylome] is **always** sorted.
 */
class MethylomeBuilder(
        private val genomeQuery: GenomeQuery,
        private val ignoreDuplicatedOffsets: Boolean
) {
    private val plusData = frameMap(genomeQuery)
    private val minusData = frameMap(genomeQuery)
    private val uniqueOffsets = genomeMap(genomeQuery) {
            BitterSet(it.length)
    }

    private val duplicatedOffsets = AtomicLong()
    fun duplicatedOffsets() = duplicatedOffsets.get()

    fun add(chromosome: Chromosome, strand: Strand,
            offset: Int, context: CytosineContext?,
            methylatedCount: Int, totalCount: Int): MethylomeBuilder {
        // Note(lebedev): this will silently skip possible "overflows". Not
        // sure if this is what we want.
        if (totalCount > 0) {
            val seenOffsets = uniqueOffsets[chromosome]
            val alreadySeen = seenOffsets[offset]
            if (alreadySeen) {
                duplicatedOffsets.incrementAndGet()
            }
            seenOffsets.set(offset)

            if (!ignoreDuplicatedOffsets || !alreadySeen) {
                val frame = strand.choose(plusData, minusData)[chromosome]
                // The magic constant is # of contexts + 1.
                frame.add(offset, context?.tag ?: 3,
                        Shorts.saturatedCast(methylatedCount.toLong()),
                        Shorts.saturatedCast(totalCount.toLong()))
            }
        }
        return this
    }

    /** Returns a sorted [Methylome]. */
    fun build(): Methylome {
        for (chromosome in genomeQuery.get()) {
            for (strand in Strand.values()) {
                val frame = strand.choose(plusData, minusData)[chromosome]
                frame.sort()
            }
        }

        return Methylome(genomeQuery, plusData, minusData)
    }
}