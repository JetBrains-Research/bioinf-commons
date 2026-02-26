package org.jetbrains.bio.genome.data

import com.google.common.collect.Maps

@Suppress("DataClassPrivateConstructor") // Constructor is used in [invoke] call
data class Cell private constructor(val name: String, val description: String) {

    override fun toString() = name

    companion object {
        private val CELLS_CACHE = Maps.newConcurrentMap<String, Cell>()

        operator fun invoke(name: String, description: String = "N/A"): Cell =
            CELLS_CACHE.computeIfAbsent(name.lowercase()) {
                Cell(name, description)
            }
    }
}
