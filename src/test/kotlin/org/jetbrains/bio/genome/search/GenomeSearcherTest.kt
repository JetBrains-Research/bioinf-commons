package org.jetbrains.bio.genome.search

import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.util.SequenceUtil
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.util.div
import org.junit.Test
import java.nio.file.Path
import kotlin.test.assertTrue

class GenomeSearcherTest {
    @Test
    fun testCorrectness() {
        val genomeQuery = Genome["to1"].toQuery()
        for (chromosome in genomeQuery.get()) {
            val fastqPath = chromosome.fastqPath
//            println("Aligning reads from $fastqPath")
            for (mismatches in 0..2) {
//                println("Performing testCorrectness with $mismatches mismatch(es).")

                val genomeSearcher = GenomeSearcher(genomeQuery, mismatches)
                for (fastqRecord in FastqReader(fastqPath.toFile())) {
                    val pattern = fastqRecord.readString
                    for (l in genomeSearcher.find(pattern)) {
                        val seed = Seed(
                            if (l.strand.isPlus())
                                pattern
                            else
                                SequenceUtil.reverseComplement(pattern), 0, 0
                        )
                        assertTrue(
                            Pattern.matches(
                                l.chromosome.sequence, seed, l.startOffset,
                                mismatches
                            )
                        )
                    }
                }
            }
        }
    }

    private data class TestCase(val pos: Int, val strand: Strand, val mismatches: Int)

    private fun String.toTestCase(): TestCase {
        val (_, pos, strand, mismatches) =
            substringAfter(' ').substringBefore(' ').split(':', limit = 4)
        return TestCase(pos.toInt(), strand.toStrand(), mismatches.toInt())
    }

    @Test
    fun testSpecificity() {
        val genomeQuery = Genome["to1"].toQuery()
        for (chromosome in genomeQuery.get()) {
            val fastqPath = chromosome.fastqPath
//            println("Aligning reads from $fastqPath")
            for (mismatches in 0..2) {
//                println("Performing testSpecificity with $mismatches mismatch(es).")
                val genomeSearcher = GenomeSearcher(genomeQuery, mismatches)
                for (fastqRecord in FastqReader(fastqPath.toFile())) {
                    val header = fastqRecord.readHeader
                    val (pos, strand, expectedMismatches) = header.toTestCase()
                    if (expectedMismatches <= mismatches) {
                        val pattern = fastqRecord.readString
                        val found = genomeSearcher.find(pattern).anyMatch { l ->
                            l.startOffset == pos &&
                                    l.chromosome === chromosome &&
                                    l.strand === strand
                        }

                        assertTrue(found, "alignment not found at expected position")
                    }
                }
            }
        }
    }

    /**
     * This files are created in [TestOrganismDataGenerator]
     */
    private val Chromosome.fastqPath: Path
        get() {
            return genome.chromSizesPath.parent!! / "sa" / "$name.fastq"
        }
}
