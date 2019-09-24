package org.jetbrains.bio.io

import com.google.common.base.CharMatcher
import com.google.common.collect.Iterators
import com.google.common.collect.PeekingIterator
import com.google.common.collect.UnmodifiableIterator
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import java.io.BufferedWriter
import java.io.IOException
import java.nio.file.Path
import java.util.*
import java.util.stream.Stream
import java.util.stream.StreamSupport

/**
 * A simple reader for FASTA files.
 *
 * See http://en.wikipedia.org/wiki/FASTA_format.
 *
 * @author Sergei Lebedev
 */
class FastaReader(val readSequence: Boolean) {
    private fun read(it: PeekingIterator<String>): UnmodifiableIterator<FastaRecord> {
        return object : UnmodifiableIterator<FastaRecord>() {
            var next: FastaRecord? = null

            override fun hasNext(): Boolean {
                if (next == null && it.hasNext()) {
                    val description = it.next()
                    check(description[0] == '>')

                    val sequence = StringBuilder()
                    while (it.hasNext() && it.peek()[0] != '>') {
                        val line = it.next()
                        if (readSequence) {
                            sequence.append(line)
                        }
                    }

                    next = FastaRecord(description.substring(1), sequence.toString())
                }

                return next != null
            }

            override fun next(): FastaRecord {
                check(hasNext())
                val fastaRecord = next
                next = null
                return fastaRecord!!
            }
        }
    }

    companion object {
        @Throws(IOException::class)
        fun read(path: Path, readSequence: Boolean = true): Stream<FastaRecord> {
            val lines = path.bufferedReader().lines()
                    .map { CharMatcher.whitespace().trimFrom(it) }
                    .filter(String::isNotEmpty)
            val iterator = FastaReader(readSequence).read(Iterators.peekingIterator(lines.iterator()))
            return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator, 0), false)
        }
    }
}

@Suppress("nothing_to_inline")
private inline fun BufferedWriter.write(ch: Char) = write(ch.toInt())

@Throws(IOException::class)
fun Iterable<FastaRecord>.write(path: Path, width: Int = 80) {
    path.bufferedWriter().use { writer ->
        forEach {
            writer.write('>')
            writer.write(it.description)
            writer.write('\n')

            val length = it.sequence.length
            for (i in 0 until length) {
                writer.write(it.sequence[i].toInt())
                if (i % width == 0 && i > 0 && i < length - 1) {
                    writer.write('\n')
                }
            }

            writer.write('\n')
        }
    }
}

data class FastaRecord(val description: String, val sequence: String)

fun Genome.writeAsFasta(path: Path) {
    this.chromosomes.asSequence()
            .map { chr -> FastaRecord(chr.name, chr.sequence.toString()) }
            .asIterable().write(path)
}