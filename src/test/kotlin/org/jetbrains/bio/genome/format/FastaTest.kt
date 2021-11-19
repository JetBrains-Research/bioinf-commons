package org.jetbrains.bio.genome.format

import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import java.util.*
import java.util.stream.Collectors
import kotlin.test.assertEquals

class FastaTest {
    @Test
    fun testWriteOne() {
        withTempFile("sample", ".fa") { path ->
            val record = FastaRecord("description", "ACGT")
            listOf(record).write(path)

            path.bufferedReader().use {
                val lines = it.lineSequence().toList()
                assertEquals(2, lines.size)
                val (description, sequence) = lines
                assertEquals(">${record.description}", description)
                assertEquals(record.sequence, sequence)
            }
        }
    }

    @Test
    fun testWriteRead() {
        val r = Random()
        val records = (0..2).map {
            FastaRecord(
                "sequence$it",
                r.nextString("ACGT", r.nextInt(20 - 1) + 1)
            )
        }

        withTempFile("random", ".fa.gz") { path ->
            records.write(path)
            assertEquals(records, FastaReader.read(path).collect(Collectors.toList()))
        }
    }
}

fun Random.nextString(alphabet: String, length: Int): String {
    return ints(length.toLong(), 0, alphabet.length)
        .mapToObj { alphabet[it].toString() }
        .collect(Collectors.joining())
}