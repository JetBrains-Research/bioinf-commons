@file:Suppress("unused")

package org.jetbrains.bio.genome.data


data class DataType(private val id: String) {

    init {
        synchronized(REGISTERED) {
            check(id !in REGISTERED) { "DataType of the type $id already registered" }
            REGISTERED[id] = this
        }
    }

    override fun toString() = id

    companion object {
        private val REGISTERED = hashMapOf<String, DataType>()

        fun fromString(id: String): DataType {
            if (id in REGISTERED) {
                return REGISTERED[id]!!
            }
            return CHIP_SEQ
        }
    }
}

// Preconfigured types
val METHYLATION = DataType("methylation")
val CHIP_SEQ = DataType("chip-seq")
val TRANSCRIPTION = DataType("transcription")
val MIRNA = DataType("mirna")

fun String.toDataType() = DataType.fromString(this)
