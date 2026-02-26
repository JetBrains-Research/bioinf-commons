package org.jetbrains.bio.genome.data

import com.google.common.collect.Maps

/** A target in the ChIP-seq protocol. */
@Suppress("DataClassPrivateConstructor", "unused") // Constructor is used in [invoke] call
data class ChipSeqTarget private constructor(val name: String) {

    override fun toString() = name

    val isInput: Boolean get() = isInput(name)

    companion object {
        private val CHIPSEQ_TARGET_CACHE = Maps.newConcurrentMap<String, ChipSeqTarget>()

        fun isInput(name: String): Boolean = "input" in name.lowercase() || "control" in name.lowercase()

        operator fun invoke(name: String): ChipSeqTarget =
            CHIPSEQ_TARGET_CACHE.computeIfAbsent(name.lowercase()) {
                ChipSeqTarget(name)
            }
    }
}
