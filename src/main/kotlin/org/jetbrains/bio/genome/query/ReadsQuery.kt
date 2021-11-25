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
 * [logFragmentSize] controls whether log messages regarding fragment size are generated.
 */
class ReadsQuery(
    val genomeQuery: GenomeQuery,
    val path: Path,
    val unique: Boolean = true,
    val fragment: Fragment = AutoFragment,
    val logFragmentSize: Boolean = true
) : CachingInputQuery<Coverage>() {

    override fun getUncached(): Coverage = coverage()

    override val description: String
        get() = "Path: $path, Unique: $unique, Fragment: $fragment"

    fun coverage(): Coverage {
        val npz = npzPath()
        npz.checkOrRecalculate("Coverage for ${path.name}") { (npzPath) ->
            val paired = isPaired(path)
            if (paired && fragment is AutoFragment) {
                PairedEndCoverage.builder(genomeQuery).apply {
                    val unpaired = processPairedReads(genomeQuery, path) { chr, pos, pnext, len ->
                        process(chr, pos, pnext, len)
                    }
                    if (unpaired != 0) {
                        LOG.info(
                            "$unpaired single-end reads encountered when reading paired-end file $path!"
                        )
                    }
                }.build(unique).save(npzPath)
            } else {
                if (paired) {
                    LOG.info("Fragment option ($fragment) forces reading paired-end reads as single-end!")
                }
                SingleEndCoverage.builder(genomeQuery).apply {
                    processReads(genomeQuery, path) {
                        process(it)
                    }
                }.build(unique).save(npzPath)
            }
        }
        val coverage = Coverage.load(npz, genomeQuery, fragment)
        val libraryDepth = coverage.depth
        if (logFragmentSize) {
            val logMessage = "Library: ${path.name}, Depth: ${"%,d".format(libraryDepth)}, " + when (coverage) {
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
            LOG.info(logMessage)
        }
        val genomeSize = genomeQuery.get().map { it.length.toLong() }.sum()
        if (libraryDepth < genomeSize * MIN_DEPTH_THRESHOLD_PERCENT / 100.0) {
            LOG.warn(
                "Library: ${path.name}, Depth: ${"%,d".format(libraryDepth)} is less than " +
                        "$MIN_DEPTH_THRESHOLD_PERCENT% x ${"%,d".format(genomeSize)} of genome ${genomeQuery.id}"
            )
        }
        return coverage
    }

    fun npzPath() = Configuration.cachePath / "coverage_${fileId}${path.sha}.npz"

    private val idStem = path.stemGz +
            (if (unique) "_unique" else "")

    override val id: String
        get() = idStem + (if (fragment is FixedFragment) "_$fragment" else "")

    // we don't need to store fragment size in the file name
    private val fileId = idStem + (if (fragment is FixedFragment) "_raw" else "")

    companion object {
        val LOG: Logger = LoggerFactory.getLogger(ReadsQuery::class.java)

        private const val MIN_DEPTH_THRESHOLD_PERCENT = 0.1
    }
}