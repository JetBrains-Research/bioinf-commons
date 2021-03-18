package org.jetbrains.bio.genome.search.sa

import org.apache.log4j.Logger
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import org.jetbrains.bio.npy.NpzFile
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.div
import org.jsuffixarrays.Algorithm
import org.jsuffixarrays.BPR
import java.nio.file.Path

/**
 * @author Roman Chernyatchik
 * @author Oleg Shpynov
 */
class SuffixArray(val sequence: NucleotideSequence, val sa: IntArray) {

    val length = sa.size

    fun index(suffix: Int): Int = sa[suffix]

    companion object {
        val LOG: Logger = Logger.getLogger(SuffixArray::class.java)
        const val VERSION = 1

        @JvmStatic
        fun create(chromosome: Chromosome) {
            val path = getPath(chromosome)
            LOG.info("Constructing SA using ${Algorithm.BPR.name}: $path")
            NpzFile.write(path).use { writer ->
                writer.write("version", intArrayOf(VERSION))
                val ints = IntArray(chromosome.length + BPR.KBS_STRING_EXTENSION_SIZE) {
                    if (it < chromosome.length) chromosome.sequence.charAt(it).toInt() else 0
                }
                val sa = Algorithm.BPR.decoratedInstance.buildSuffixArray(ints, 0, chromosome.sequence.length)
                writer.write("sa", sa)
                LOG.info("Done constructing SA $path")
            }
        }

        @JvmStatic
        fun load(chromosome: Chromosome): SuffixArray {
            val path = getPath(chromosome)
            LOG.info("Loading suffix array $path")
            NpzFile.read(path).use { reader ->
                val version = reader["version"].asIntArray().single()
                check(version == VERSION) { "index version is $version instead of $VERSION" }
                val sa = reader["sa"].asIntArray()
                val suffixArray = SuffixArray(chromosome.sequence, sa)
                LOG.info("Loaded suffix array $path")
                return suffixArray
            }
        }

        fun getPath(chromosome: Chromosome): Path {
            val indicesPath = chromosome.genome.chromSizesPath.parent / "sa"
            indicesPath.createDirectories()
            return indicesPath / "${chromosome.name}.npz"
        }
    }
}