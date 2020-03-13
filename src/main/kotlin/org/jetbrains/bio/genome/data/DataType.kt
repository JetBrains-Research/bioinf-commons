package org.jetbrains.bio.genome.data

import com.google.common.collect.Maps


@Suppress("DataClassPrivateConstructor") // Constructor is used in [invoke] call
data class DataType private constructor(val id: String) {

    override fun toString() = id

    companion object {
        private val CACHE = Maps.newConcurrentMap<String, DataType>()

        operator fun invoke(id: String): DataType =
                CACHE.computeIfAbsent(id.toLowerCase()) {
                    DataType(id)
                }
    }
}

