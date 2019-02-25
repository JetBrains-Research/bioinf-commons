package org.jetbrains.bio.methylome

import com.google.common.base.Joiner
import com.google.common.collect.ImmutableList
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.byByte
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.StrandFilter

/**
 * A methylome function converts a [Methylome] instance into an
 * [DataFrame].
 *
 * All of the functions are stateless and can be used simultaneously in
 * a concurrent setting.
 *
 * @author Sergei Lebedev
 */
class MethylomeToDataFrame private constructor(
        val strandFilter: StrandFilter,
        val cytosineContext: CytosineContext?) {

    fun apply(methylome: Methylome, chromosome: Chromosome): DataFrame {
        return apply(ImmutableList.of<Methylome>(methylome), chromosome)
    }

    fun apply(methylomes: List<Methylome>, chromosome: Chromosome): DataFrame {
        check(methylomes.isNotEmpty()) { "need at least one methylome" }
        val dfs = methylomes.map { filter(get(it, chromosome)) }.toTypedArray()
        if (dfs.size == 1) {
            return apply(dfs[0])
        }

        // Tag and strand for the combined methylomes should be the same.
        var res = DataFrame.mergeInner("offset", *dfs)
        for (i in dfs.indices) {
            res = res.rename("tag${i + 1}".intern(), "tag")

            val strand = "strand${i + 1}".intern()
            if (strand in res.labels) {
                res = res.rename(strand, "strand")
            }
        }

        return apply(res)
    }

    private fun DataFrame.rename(from: String, to: String): DataFrame {
        return with(rowsNumber, get(from).rename(to)).omit(from)
    }

    private fun get(methylome: Methylome, chromosome: Chromosome): MethylomeView {
        return when (strandFilter) {
            StrandFilter.PLUS  -> methylome[chromosome, Strand.PLUS]
            StrandFilter.MINUS -> methylome[chromosome, Strand.MINUS]
            StrandFilter.BOTH  -> methylome.getCombined(chromosome)
        }
    }

    private fun apply(df: DataFrame): DataFrame {
        if (df.rowsNumber == 0) {
            return df.with("d", IntArray(0))
        }

        val distances = df.sliceAsInt("offset").clone()
        for (t in distances.size - 1 downTo 1) {
            distances[t] -= distances[t - 1]
        }

        distances[0] = Integer.MAX_VALUE // should be unused.
        return df.with("d", distances)
    }

    private fun filter(methylomeView: MethylomeView): DataFrame {
        val df = methylomeView.peel()
        return if (cytosineContext == null) {
            df
        } else {
            df.filter(byByte("tag") { it == cytosineContext.tag })
        }
    }

    override fun toString(): String = Joiner.on('_')
            .join(strandFilter, cytosineContext?.toString() ?: "ANY")

    companion object {
        @JvmStatic val DEFAULT = create(StrandFilter.PLUS, CytosineContext.CG)

        @JvmStatic fun create(strandFilter: StrandFilter,
                              cytosineContext: CytosineContext? = null): MethylomeToDataFrame {
            return MethylomeToDataFrame(strandFilter, cytosineContext)
        }
    }
}
