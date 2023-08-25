package org.jetbrains.bio.genome.format

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import java.util.*
import java.util.stream.Collectors
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

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

    @Test
    fun testFastaRecordWrite() {
        val record1 = FastaRecord("chr1", "a".repeat(15))
        val record2 = FastaRecord("chr2", "c".repeat(18))
        val record3 = FastaRecord("chr3", "t".repeat(24))
        val sequence = listOf(record1, record2, record3)
        withTempFile("fasta", ".fa") {
            sequence.write(it, 8)
            val lines = it.toFile().readLines()
            assertEquals(11, lines.size)
            assertEquals(">chr1", lines[0])
            assertEquals("a".repeat(8), lines[1])
            assertEquals("a".repeat(7), lines[2])
            assertEquals(">chr2", lines[3])
            assertEquals("c".repeat(8), lines[4])
            assertEquals("c".repeat(8), lines[5])
            assertEquals("c".repeat(2), lines[6])
            assertEquals(">chr3", lines[7])
            assertEquals("t".repeat(8), lines[8])
            assertEquals("t".repeat(8), lines[9])
            assertEquals("t".repeat(8), lines[10])
        }
    }

    @Test
    fun testGenomeWriteAsFasta() {
        val testGenome = Genome["to1"]

        withTempFile("fasta", ".fa") { fastaFile ->
            testGenome.writeAsFasta(fastaFile)
            val lines = fastaFile.toFile().readLines()
            val chromosomes = ArrayList<String>()
            val chromosomeLengths = ArrayList<Int>()
            var width = -1
            for (pair in lines zip lines.indices) {
                val line = pair.first
                if (line.startsWith(">")) {
                    chromosomes.add(line.substring(1, line.length))
                    assertTrue(testGenome.chromosomeNamesMap.containsKey(chromosomes.last()))
                    if (pair.second != lines.size - 1) {
                        for (i in lines[pair.second + 1].indices) {
                            if (lines[pair.second + 1][i].uppercaseChar()
                                != testGenome.chromosomeNamesMap[chromosomes.last()]!!.sequence.charAt(i).uppercaseChar()) {
                                fail("Saved file has different bases than stored in Genome")
                            }
                        }
                    }
                    chromosomeLengths.add(0)
                } else {
                    if (chromosomeLengths.isEmpty()) fail("There is no correct beginning for chromosome sequence")
                    if (width == -1) width = line.length
                    if (line.length > width) fail("Lines in fasta file does not have the same width (at least one line is longer that first one)")
                    if (line.length < width) {
                        if (pair.second != lines.size - 1
                            && !lines[pair.second + 1].startsWith(">")) {
                            fail("Lines in fasta file does not have the same width (at least one is shorter than first one)")
                        }
                    }
                    chromosomeLengths[chromosomeLengths.size - 1] = chromosomeLengths[chromosomeLengths.size - 1] + line.length
                }
            }
            for (pair in chromosomes zip testGenome.chromSizesMap.keys) {
                assertEquals(pair.second, pair.first)
            }
            var totalLength = 0
            for (pair in chromosomeLengths zip testGenome.chromSizesMap.values) {
                assertEquals(pair.second, pair.first)
                totalLength += pair.first
            }
        }
    }
}

fun Random.nextString(alphabet: String, length: Int): String {
    return ints(length.toLong(), 0, alphabet.length)
        .mapToObj { alphabet[it].toString() }
        .collect(Collectors.joining())
}