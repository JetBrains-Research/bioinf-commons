package org.jetbrains.bio.genome.coverage

import gnu.trove.list.TIntList
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.npy.NpzFile
import java.io.IOException
import java.nio.file.Path


/**
 * [Coverage] allows to query the tags coverage of a given genomic [Location].
 * Currently, it has two implementations: [SingleEndCoverage], [PairedEndCoverage], [BasePairCoverage].
 */
interface Coverage {

    val genomeQuery: GenomeQuery

    /**
     * Returns the number of tags inside a given [location].
     */
    fun getCoverage(location: Location): Int

    /**
     * Returns the number of tags inside a given [chromosomeRange]
     * (i.e. on both strands).
     */
    fun getBothStrandsCoverage(chromosomeRange: ChromosomeRange): Int =
        getCoverage(chromosomeRange.on(Strand.PLUS)) +
                getCoverage(chromosomeRange.on(Strand.MINUS))

    val depth: Long

    enum class CoverageType {
        SINGLE_END_CHIPSEQ,
        PAIRED_END_CHIPSEQ,
        BASEPAIR,
    }
    companion object {

        /**
         * Binary storage format version. Loader will throw an [IOException]
         * if it doesn't match.
         */
        const val VERSION = 5
        const val VERSION_FIELD = "version"
        const val COV_TYPE_FIELD = "coverage_type"

        @Throws(IOException::class)
        fun load(
            inputPath: Path,
            genomeQuery: GenomeQuery,
            fragment: Fragment = AutoFragment,
            failOnMissingChromosomes: Boolean = false
        ): Coverage {
            return NpzFile.read(inputPath).use { reader ->
                val version = reader[VERSION_FIELD].asIntArray().single()
                val coverageType = if (version == 4 ) {
                    val paired = reader["paired"].asBooleanArray().single()
                    if (paired) CoverageType.PAIRED_END_CHIPSEQ else CoverageType.SINGLE_END_CHIPSEQ
                } else {
                    check(version == VERSION) {
                        "$inputPath coverage version is $version instead of $VERSION"
                    }

                    val covTypeByte = reader[COV_TYPE_FIELD].asIntArray().single()
                    require(covTypeByte >= 0)
                    require(covTypeByte < CoverageType.values().size)
                    CoverageType.values()[covTypeByte]
                }

                when (coverageType) {
                    CoverageType.SINGLE_END_CHIPSEQ -> SingleEndCoverage.load(reader, inputPath, genomeQuery, failOnMissingChromosomes)
                        .withFragment(fragment)
                    CoverageType.PAIRED_END_CHIPSEQ -> PairedEndCoverage.load(reader, inputPath, genomeQuery, failOnMissingChromosomes)
                    CoverageType.BASEPAIR -> BasePairCoverage.load(reader, inputPath, genomeQuery, failOnMissingChromosomes)
                }
            }
        }
    }
}

/**
 * Returns the insertion index of [target] into a sorted list.
 *
 * If [target] already appears in this vector, the returned
 * index is just before the leftmost occurrence of [target].
 */
fun TIntList.binarySearchLeft(target: Int): Int {
    var lo = 0
    var hi = size()
    while (lo < hi) {
        val mid = (lo + hi) ushr 1
        if (target <= this[mid]) {
            hi = mid
        } else {
            lo = mid + 1
        }
    }

    return lo
}