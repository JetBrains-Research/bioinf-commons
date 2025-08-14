package org.jetbrains.bio.genome.format

import org.jetbrains.bio.util.extension
import org.slf4j.LoggerFactory
import java.nio.file.Path

enum class ReadsFormat {
    BAM,
    SAM,
    CRAM,
    BED;

    fun check(path: Path) {
        when (this) {
            BED -> {
                if (path.extension.lowercase() !in listOf("bed", "gz", "zip")) {
                    LOG.warn("Unexpected file extension for BED file: $path, expected bed, bed.gz or bed.zip")
                }
            }

            BAM -> {
                if (path.extension.lowercase() != "bam") {
                    LOG.warn("Unexpected file extension for BAM file: $path, expected bam")
                }
            }

            SAM -> {
                if (path.extension.lowercase() != "sam") {
                    LOG.warn("Unexpected file extension for SAM file: $path, expected sam")
                }
            }

            CRAM -> {
                if (path.extension.lowercase() != "cram") {
                    LOG.warn("Unexpected file extension for CRAM file: $path, expected cram")
                }
            }
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(ReadsFormat::class.java)


        fun guess(path: Path): ReadsFormat? = when (path.extension.lowercase()) {
            in listOf("bed", "gz", "zip") -> {
                BED
            }

            "bam" -> BAM
            "sam" -> SAM
            "cram" -> CRAM
            else -> null
        }
    }

}