package org.jetbrains.bio.io

import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools.cram.ref.ReferenceSource
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.io.File
import java.util.*

/**
 * A fake [htsjdk.samtools.cram.ref.ReferenceSource] which stores bytes
 * for a single [Chromosome].
 */
internal class WrappedReferenceSource(private val name: String,
                                      sequence: NucleotideSequence)
:
        ReferenceSource(null as File?) {

    // XXX please keep lazy to reduce memory consumption.
    private val bytes: ByteArray by lazy(LazyThreadSafetyMode.PUBLICATION) {
        sequence.toString().toByteArray()
    }

    override fun getReferenceBases(record: SAMSequenceRecord,
                                   tryNameVariants: Boolean): ByteArray? {
        if (tryNameVariants) {
            for (variant in getVariants(record.sequenceName)) {
                if (variant == name) {
                    return bytes
                }
            }
        } else {
            if (record.sequenceName == name) {
                return bytes
            }
        }

        return null
    }

    internal fun getVariants(name: String): List<String> {
        val variants = ArrayList<String>()
        when {
            name == "M"  -> variants.add("MT")
            name == "MT" -> variants.add("M")
            name.startsWith("chr") -> variants.add(name.substring(3))
            else -> variants.add("chr" + name)
        }

        if (name == "chrM") {
            // chrM case:
            variants.add("MT")
        }

        return variants
    }
}