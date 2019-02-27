package org.jetbrains.bio.query

import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.coverage.*
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.*
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
        val fragment: Int? = null, // ignored if paired
        val logFragmentSize: Boolean = true
): CachingInputQuery<Coverage>() {

    override fun getUncached(): Coverage = coverage()

    override val description: String
        get() = "Path: $path, Unique: $unique, Fragment: $fragment"

    fun coverage(): Coverage {
        val npz = npzPath()
        npz.checkOrRecalculate("Coverage for ${path.name}") { (npzPath) ->
            val paired = isPaired(path)
            if (paired && fragment == null) {
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
        if (logFragmentSize) {
            val logMessage = "Library: ${path.name}, Depth: ${coverage.depth}, " + when (coverage) {
                is SingleEndCoverage -> "Reads: single-ended, " + when {
                    coverage.actualFragment != coverage.detectedFragment ->
                        "Fragment size: ${coverage.actualFragment} bp " +
                                "(overrides cross-correlation " +
                                "estimate ${coverage.detectedFragment})"
                    fragment == null ->
                        "Fragment size: ${coverage.detectedFragment} bp " +
                                "(cross-correlation estimate)"
                    else ->
                        "Fragment size: ${coverage.detectedFragment} bp " +
                                "(user input is equal to cross-correlation estimate)"
                }
                is PairedEndCoverage -> "Reads: paired-ended, " + if (fragment == null) {
                    "Fragment size: ${coverage.averageInsertSize} bp (average; inferred from read pairs)"
                } else {
                    "Fragment size: ${coverage.averageInsertSize} bp (average; inferred from read pairs; " +
                            "user input $fragment is ignored)"
                }
                else -> throw IllegalArgumentException("Unknown library type: ${coverage::class.java}")
            }
            LOG.info(logMessage)
        }
        return coverage
    }

    fun npzPath() = Configuration.cachePath /  "coverage_${fileId}${path.sha}.npz"

    override val id: String
        get() = fileId +  (if (fragment != null) "_$fragment" else "")

    // we don't need to store fragment size in the file name
    private val fileId: String = path.stemGz + (if (unique) "_unique" else "")

    companion object {
        val LOG: Logger = Logger.getLogger(ReadsQuery::class.java)
    }
}

/**
 * One might be tempted to use [substringBeforeLast] here.
 * However, this method can't fulfill the following requirements:
 * we want to _ignore_ the ".bed.gz" extension case but _preserve_ the stem case.
 * The current trick with string truncation does exactly that.
 */
val Path.stemGz: String get() {
    return when {
        name.toLowerCase().endsWith(".bed.gz") ->
            name.substring(0, name.length - ".bed.gz".length)
        else -> stem
    }
}
