package org.jetbrains.bio.genome.query

import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.*
import org.jetbrains.bio.genome.format.isPaired
import org.jetbrains.bio.genome.format.processPairedReads
import org.jetbrains.bio.genome.format.processReads
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Query of tags coverage created from a BAM or BED/BED.GZ file.
 *
 * [unique] controls whether duplicate tags should be preserved ([unique] == false)
 * or squished into one tag ([unique] == true).
 *
 * Not-null [fragment] overrides the imputed fragment size of a single-end coverage.
 * If null, the imputed value is used.
 * If the coverage is paired-end, any value of [fragment] is ignored.
 *
 * [showLibraryInfo] controls whether log messages regarding fragment size are generated.
 */
class ReadsQuery(
    val genomeQuery: GenomeQuery,
    val path: Path,
    val unique: Boolean = true,
    val fragment: Fragment = AutoFragment,
    val showLibraryInfo: Boolean = true,
) : CachingInputQuery<Coverage>() {

    override fun getUncached(): Coverage = coverage()


    /**
     * Unique id is formed
     */
    override val id: String
        get() = idStem + (if (fragment is FixedFragment) "_$fragment" else "")

    override val description: String
        get() = "Path: $path, Unique: $unique, Fragment: $fragment"

    /**
     * Coverage is stored as Numpy array, with tags positions after optional fragment size shift
     */
    fun npzPath() =
        Configuration.cachePath / "coverage_${
            idStem + (if (fragment is FixedFragment) "_raw" else "")
        }${path.sha}.npz"

    /**
     * Unique identifier for file and
     */
    private val idStem = path.stemGz + (if (unique) "_unique" else "")

    /**
     * @return true when preprocessed coverage file is available or it can be recomputed by raw file
     */
    fun isAccessible(): Boolean = npzPath().isAccessible() || path.isAccessible()

    fun coverage(): Coverage {
        check(isAccessible()) { "Unable to load coverage, original reads and cache not available" }
        val npz = npzPath()
        npz.checkOrRecalculate("Coverage for ${path.name}") { (npzPath) ->
            val paired = isPaired(path)
            if (paired && fragment is AutoFragment) {
                val pairedEndCoverage = PairedEndCoverage.builder(genomeQuery).apply {
                    val unpaired = processPairedReads(genomeQuery, path) { chr, pos, pnext, len ->
                        process(chr, pos, pnext, len)
                    }
                    if (unpaired != 0) {
                        LOG.info(
                            "$unpaired single-end reads encountered when reading paired-end file $path!"
                        )
                    }
                }.build(unique)
                pairedEndCoverage.save(npzPath)
            } else {
                if (paired) {
                    LOG.info("Fragment option ($fragment) forces reading paired-end reads as single-end!")
                }
                val singleEndCoverage = SingleEndCoverage.builder(genomeQuery).apply {
                    processReads(genomeQuery, path) {
                        process(it)
                    }
                }.build(unique)
                singleEndCoverage.save(npzPath)
            }
        }
        val coverage = Coverage.load(npz, genomeQuery, fragment)
        if (showLibraryInfo) {
            showLibraryInfo(coverage)
        }
        return coverage
    }

    private fun showLibraryInfo(coverage: Coverage) {
        val libraryDepth = coverage.depth
        val information = "Library: ${path.name}, Depth: ${"%,d".format(libraryDepth)}, " +
                when (coverage) {
                    is SingleEndCoverage -> "Reads: single-ended, " + when {
                        coverage.actualFragment != coverage.detectedFragment ->
                            "Fragment size: ${coverage.actualFragment} bp " +
                                    "(overrides cross-correlation " +
                                    "estimate ${coverage.detectedFragment})"
                        fragment is AutoFragment ->
                            "Fragment size: ${coverage.detectedFragment} bp " +
                                    "(cross-correlation estimate)"
                        else ->
                            "Fragment size: ${coverage.detectedFragment} bp " +
                                    "(user input is equal to cross-correlation estimate)"
                    }
                    is PairedEndCoverage -> "Reads: paired-ended, " + if (fragment is FixedFragment) {
                        "Fragment size: ${coverage.averageInsertSize} bp (average; inferred from read pairs; " +
                                "user input $fragment is ignored)"
                    } else {
                        "Fragment size: ${coverage.averageInsertSize} bp (average; inferred from read pairs)"
                    }
                    else -> throw IllegalArgumentException("Unknown library type: ${coverage::class.java}")
                }
        LOG.info(information)
    }

    companion object {
        val LOG: Logger = LoggerFactory.getLogger(ReadsQuery::class.java)
    }
}