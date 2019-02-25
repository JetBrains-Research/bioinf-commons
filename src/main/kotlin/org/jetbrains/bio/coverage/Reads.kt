package org.jetbrains.bio.coverage

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.toStrand
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.io.unpackRegularFields
import org.jetbrains.bio.util.*
import picard.sam.markduplicates.MarkDuplicates
import java.nio.file.Path

/**
 * Attempts to detect whether the file contains paired-end reads
 * by looking at the first valid read.
 * Currently only BAM and CRAM files are supported for paired-end processing.
 */
fun isPaired(path: Path): Boolean {
    return when (path.extension) {
        "bam", "cram" ->
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
fun processReads(genomeQuery: GenomeQuery, path: Path, consumer: (Location) -> Unit) {
    val progress = Progress { title = "Loading reads ${path.name}" }.unbounded()
    try {
        when (path.extension) {
            //     vvv this is silly, yes.
            "bed", "gz", "zip" -> {
                val format = BedFormat.auto(path)
                format.parse(path) {
                    it.forEach { entry ->
                        val chromosome = genomeQuery[entry.chrom]
                        if (chromosome != null) {
                            val e = entry.unpackRegularFields(format)
                            val strand = e.strand.toStrand()
                            progress.report()
                            consumer(Location(
                                e.start, Math.max(e.start + 1, e.end), chromosome, strand
                            ))
                        }
                    }
                }
            }

            "bam", "cram" -> {
                //XXX: [for future] Previous version was processing file in parallel by chromosome
                // it required indexed bam file + provide noticeable speedup (in case of methylome parsing)
                // Parsing result is normally cached to *.npz file, so it isn't a bottle neck at the moment.
                SamReaderFactory.make().validationStringency(ValidationStringency.SILENT).open(path.toFile()).use {
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
            else -> error("unsupported file type: $path")
        }
    } finally {
        progress.done()
    }
}

/**
 * Extract all (valid) read pairs from BAM and feed them to [consumer].
 * Pairs that don't belong to the provided [genomeQuery] are ignored.
 * Consumer accepts the first read of pair and the inferred insert size.
 *
 * Returns the number of valid unpaired reads encountered. If it's not zero,
 * something very wrong has happened.
 */
fun processPairedReads(
        genomeQuery: GenomeQuery, path: Path, consumer: (Location, Int) -> Unit
): Int {
    val progress = Progress { title = "Loading paired-end reads ${path.name}" }.unbounded()
    try {
        var unpairedCount = 0
        when (path.extension) {
            "bam", "cram" -> {
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
                                if (record.secondOfPairFlag) {
                                    // only iterate each read pair once but increment progress
                                    // so that the user doesn't think we're stuck
                                    return@forEach
                                }
                                val location = record.toLocation(genomeQuery)
                                if (location != null) {
                                    consumer(location, record.inferredInsertSize)
                                }
                            }
                        }
            }
            else -> error("unsupported file type: $path")
        }
        return unpairedCount
    } catch (e: IllegalStateException) {
        error("Error when processing paired-end reads: " +
                Logs.getMessage(e, includeStackTrace = true))
    } finally {
        progress.done()
    }
}

fun removeDuplicates(path: Path): Path {
    check(path.extension == "bam") {
        "Only BAM supported, got: $path"
    }
    val uniqueReads = path.toString().replace(".bam", "_unique.bam").toPath()
    uniqueReads.checkOrRecalculate("Unique reads") { output ->
        MarkDuplicates.main(arrayOf(
                "REMOVE_DUPLICATES=true",
                "INPUT=$path",
                "OUTPUT=${output.path}",
                "M=${path.toString().replace(".bam", "_metrics.txt")}"))
    }
    return uniqueReads
}

/**
 * Checks the most common reasons to ignore a BAM record.
 */
private fun SAMRecord.invalid(): Boolean = readUnmappedFlag
        || isSecondaryOrSupplementary
        || duplicateReadFlag
        || mappingQuality == 0 // BWA multi-alignment evidence
        || alignmentStart == 0

private fun SAMRecord.toLocation(genomeQuery: GenomeQuery): Location? =
        // 1 based, end inclusive
        genomeQuery[referenceName]?.let { chromosome ->
            Location(
                alignmentStart - 1, alignmentEnd,
                chromosome,
                if (readNegativeStrandFlag) Strand.MINUS else Strand.PLUS
            ) }