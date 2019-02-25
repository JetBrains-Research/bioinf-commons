package org.jetbrains.bio.dataset

import com.google.common.collect.Maps

/**
 * @author Oleg Shpynov
 * @since 21/07/2017.
 */
enum class DataType(val id: String) {
    METHYLATION("methylation"),
    CHIP_SEQ("chip-seq"),
    TRANSCRIPTION("transcription"),
    MIRNA("mirna")
}

fun String.toDataType(): DataType {
    return DataType.values()
            .firstOrNull { this.toLowerCase() == it.id.toLowerCase() }
            ?: DataType.CHIP_SEQ
}

/** A target in the ChIP-seq protocol. */
@Suppress("DataClassPrivateConstructor", "unused")
data class ChipSeqTarget private constructor(val name: String) {

    override fun toString() = name

    val isInput: Boolean get() = isInput(name)

    companion object {
        private val CACHE = Maps.newConcurrentMap<String, ChipSeqTarget>()

        fun isInput(name: String): Boolean = "input" in name.toLowerCase()

        fun isWide(name: String): Boolean = (name == "H3K27me3" || name == "H3K36me3")

        init {
            // Touch classes to Load them in JVM & perform auto-registration
            HumanCells.toString()
        }

        operator fun invoke(name: String): ChipSeqTarget =
                CACHE.computeIfAbsent(name.toLowerCase()) {
                    ChipSeqTarget(name)
                }

        val H3K4me1 = ChipSeqTarget("H3K4me1")
        val H3K4me2 = ChipSeqTarget("H3K4me2")
        val H3K4me3 = ChipSeqTarget("H3K4me3")
        val H3K9ac = ChipSeqTarget("H3K9ac")
        val H3K9me3 = ChipSeqTarget("H3K9me3")
        val H3K27ac = ChipSeqTarget("H3K27ac")
        val H3K27me3 = ChipSeqTarget("H3K27me3")
        val H3K36me3 = ChipSeqTarget("H3K36me3")
        val H4K20me1 = ChipSeqTarget("H4K20me1")

        val CTCF = ChipSeqTarget("CTCF")
        val DNAse = ChipSeqTarget("DNAse")

        val Input = ChipSeqTarget("Input")
    }
}
