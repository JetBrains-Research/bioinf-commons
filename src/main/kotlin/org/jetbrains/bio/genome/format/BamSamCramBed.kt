package org.jetbrains.bio.genome.format

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.util.*
import picard.sam.markduplicates.MarkDuplicates
import java.nio.file.Path
import kotlin.math.max


/**
 * Attempts to detect whether the file contains paired-end reads
 * by looking at the first valid read.
 * Currently only BAM and CRAM files are supported for paired-end processing.
 */
fun isPaired(path: Path, explicitFormat: ReadsFormat?): Boolean {
    val readsFormat = explicitFormat ?: ReadsFormat.guess(path)
    if (readsFormat == null) {
        error("unknown format $path, please specify format explicitly with --format option")
    }
    readsFormat.check(path)
    return when (readsFormat) {
        ReadsFormat.BAM, ReadsFormat.SAM, ReadsFormat.CRAM ->
            SamReaderFactory.make()
                .validationStringency(ValidationStringency.SILENT)
                .open(path.toFile()).use {
                    it.forEach { record ->
                        if (record.invalid()) {
                            return@forEach
                        }
                        return@use record.readPairedFlag
                    }
                    return@use false
                }

        else -> false
    }
}

/**
 * Extract all (valid) reads from BAM or BED[.gz] and feed them to [consumer].
 * Reads that don't belong to the provided [genomeQuery] are ignored.
 */
fun processReads(
    genomeQuery: GenomeQuery,
    path: Path,
    explicitFormat: ReadsFormat?,
    consumer: (Location) -> Unit
) {
    val readsFormat = explicitFormat ?: ReadsFormat.guess(path)
    if (readsFormat == null) {
        error("unknown format $path, please specify format explicitly with --format option")
    }
    readsFormat.check(path)
    val progress = Progress { title = "Loading reads ${path.name}" }.unbounded()
    try {
        when (readsFormat) {
            ReadsFormat.BED -> {
                processBed(path, genomeQuery, progress, consumer)
            }

            ReadsFormat.BAM, ReadsFormat.SAM, ReadsFormat.CRAM -> {
                processBamSamCram(path, genomeQuery, progress, consumer)
            }
        }
    } finally {
        progress.done()
    }
}

private fun processBamSamCram(
    bamSamCramPath: Path,
    genomeQuery: GenomeQuery,
    progress: Progress,
    consumer: (Location) -> Unit
) {
    //XXX: [for future] Previous version was processing file in parallel by chromosome
    // it required indexed bam file + provide noticeable speedup (in case of methylome parsing)
    // Parsing result is normally cached to *.npz file, so it isn't a bottle neck at the moment.
    SamReaderFactory.make().validationStringency(ValidationStringency.SILENT).open(bamSamCramPath.toFile()).use {
        it.forEach { record ->
            if (record.invalid()) {
                return@forEach
            }
            val location = record.toLocation(genomeQuery)
            progress.report()
            if (location != null) {
                consumer(location)
            }
        }
    }
}

private fun processBed(
    bedPath: Path,
    genomeQuery: GenomeQuery,
    progress: Progress,
    consumer: (Location) -> Unit
) {
    val format = BedFormat.auto(bedPath)
    format.parse(bedPath) {
        it.forEach { entry ->
            val chromosome = genomeQuery[entry.chrom]
            if (chromosome != null) {
                val e = entry.unpackRegularFields(format)
                val strand = e.strand.toStrand()
                progress.report()
                consumer(
                    Location(
                        e.start, max(e.start + 1, e.end), chromosome, strand
                    )
                )
            }
        }
    }
}

/**
 * Extract all (valid) read pairs from BAM and feed them to [consumer].
 * Only the negative strand read is passed to [consumer].
 * Pairs that don't belong to the provided [genomeQuery] are ignored.
 * Consumer accepts:
 * - the chromosome;
 * - POS (leftmost mapped position of the read);
 * - PNEXT (leftmost mapped position of the mate read);
 * - length of the read.
 *
 * Returns the number of valid unpaired reads encountered. If it's not zero,
 * something very wrong has happened.
 */
fun processPairedReads(
    genomeQuery: GenomeQuery,
    path: Path,
    explicitFormat: ReadsFormat?,
    consumer: (Chromosome, Int, Int, Int) -> Unit
): Int {
    val readsFormat = explicitFormat ?: ReadsFormat.guess(path)
    if (readsFormat == null) {
        error("unknown format $path, please specify format explicitly with --format option")
    }
    check(readsFormat in listOf(ReadsFormat.BAM, ReadsFormat.SAM, ReadsFormat.CRAM)) {
        "Only BAM supported, got: $path"
    }
    val progress = Progress { title = "Loading paired-end reads ${path.name}" }.unbounded()
    try {
        var unpairedCount = 0
        SamReaderFactory.make()
            .validationStringency(ValidationStringency.SILENT)
            .open(path.toFile()).use { iterationReader ->
                iterationReader.forEach { record ->
                    if (record.invalid()) {
                        return@forEach
                    }
                    progress.report()
                    if (!record.readPairedFlag) {
                        /* this really, really shouldn't happen,
                        * but it's nice to have a safety net */
                        unpairedCount++
                        return@forEach
                    }
                    if (record.mateUnmappedFlag) {
                        // we skip partially mapped pairs
                        return@forEach
                    }
                    if (!record.readNegativeStrandFlag || record.mateNegativeStrandFlag
                        || record.referenceName != record.mateReferenceName
                    ) {
                        // We only process negative strand reads.
                        // We skip same-strand pairs because we can't easily
                        // infer insert size for them.
                        // We also skip pairs mapped to different chromosomes.
                        return@forEach
                    }

                    val pos = record.alignmentStart
                    val pnext = record.mateAlignmentStart
                    val length = record.readLength
                    val chromosome = genomeQuery[record.referenceName]
                    if (pnext != 0 && length != 0 && chromosome != null) {
                        consumer(
                            chromosome, pos, pnext, length
                        )
                    }
                }
            }
        return unpairedCount
    } catch (e: IllegalStateException) {
        val message = Logs.getMessage(e, includeStackTrace = true)
        error("Error when processing paired-end reads: $message")
    } finally {
        progress.done()
    }
}

fun removeDuplicates(path: Path, explicitFormat: ReadsFormat?): Path {
    val readsFormat = explicitFormat ?: ReadsFormat.guess(path)
    if (readsFormat == null) {
        error("unknown format $path, please specify format explicitly with --format option")
    }
    check(readsFormat == ReadsFormat.BAM) {
        "Only BAM supported, got: $path"
    }
    val uniqueReads = path.toString().replace(".bam", "_unique.bam").toPath()
    uniqueReads.checkOrRecalculate("Unique reads") { output ->
        MarkDuplicates.main(
            arrayOf(
                "REMOVE_DUPLICATES=true",
                "INPUT=$path",
                "OUTPUT=${output.path}",
                "M=${path.toString().replace(".bam", "_metrics.txt")}"
            )
        )
    }
    return uniqueReads
}

/**
 * Checks the most common reasons to ignore a BAM record.
 */
private fun SAMRecord.invalid(): Boolean = readUnmappedFlag
        || isSecondaryOrSupplementary
        || duplicateReadFlag
        || alignmentStart == 0

private fun SAMRecord.toLocation(genomeQuery: GenomeQuery): Location? =
    // 1 based, end inclusive
    genomeQuery[referenceName]?.let { chromosome ->
        Location(
            alignmentStart - 1, alignmentEnd,
            chromosome,
            if (readNegativeStrandFlag) Strand.MINUS else Strand.PLUS
        )
    }