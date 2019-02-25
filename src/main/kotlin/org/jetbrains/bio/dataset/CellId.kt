package org.jetbrains.bio.dataset

import com.google.common.collect.Maps

@Suppress("DataClassPrivateConstructor") // Constructor is used in [invoke] call
data class CellId private constructor(val name: String, val description: String) {

    override fun toString() = name

    companion object {
        private val CACHE = Maps.newConcurrentMap<String, CellId>()

        init {
            // Touch classes to Load them in JVM & perform auto-registration
            HumanCells.toString()
        }

        operator fun invoke(name: String, description: String = "N/A"): CellId =
                CACHE.computeIfAbsent(name.toLowerCase()) {
                    CellId(name, description)
                }
    }
}

object HumanCells {
    val H1: CellId = CellId("h1", "Human H1 male embrionic stem cells")
    val H9: CellId = CellId("h9", "Human H9 male embrionic stem cells")
    val IMR90: CellId = CellId("imr90", "Human female lung fibroblasts cells")
    val AL: CellId = CellId("AL", "Human atherosclerotic lesion")
    val AO: CellId = CellId("AO", "Human aortic tissue")
}