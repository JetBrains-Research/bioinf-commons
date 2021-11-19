package org.jetbrains.bio.genome.format

import org.jetbrains.bio.util.name
import java.nio.file.Path


class FastqReads(
    val id: String,
    val reads1: List<Path>,
    val reads2: List<Path>? = null
) {
    init {
        require(reads1.isNotEmpty()) { "no data" }
        require(reads2 == null || reads1.size == reads2.size)
    }

    val isPaired: Boolean
        get() = reads2 != null

    fun toPairedArray(): Array<Path> {
        return if (reads2 != null) {
            reads1.zip(reads2).flatMap { it.toList() }.toTypedArray()
        } else {
            reads1.toTypedArray()
        }
    }

    fun getReadGroups(): List<List<Path>> {
        return if (isPaired) {
            reads1.zip(reads2!!).map { it.toList() }
        } else {
            reads1.map { listOf(it) }
        }
    }

    val path: Path get() = reads1.first().parent

}

fun createFastqReads(id: String, fastqReads: List<Path>): FastqReads {
    val fastqReads1 = fastqReads.filter { "_1.f" in it.name }.sorted()
    val fastqReads2 = fastqReads.filter { "_2.f" in it.name }.sorted()

    if (fastqReads2.isEmpty()) {
        return FastqReads(id, fastqReads, null)
    } else {
        require(fastqReads1.size == fastqReads2.size) { "Different numbers of paired reads" }
        require(fastqReads1.size + fastqReads2.size == fastqReads.size)

        return FastqReads(id, fastqReads1, fastqReads2)
    }
}
