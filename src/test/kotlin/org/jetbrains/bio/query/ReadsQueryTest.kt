package org.jetbrains.bio.query

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.coverage.PairedEndCoverage
import org.jetbrains.bio.coverage.SingleEndCoverage
import org.jetbrains.bio.coverage.processPairedReads
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.*
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse

/**
 * @author Oleg Shpynov
 * @author Aleksei Dievskii
 * @since 11/04/2018.
 *
 * The test files, "single_end.bam" and "paired_end.bam" were created in the following way:
 * 1. take a real world single- or paired-end BAM, say `file.bam`
 * 2. run `samtools view -c file.bam` to determine the number of reads
 * 3. divide the number of reads that you want in the test file
 *      (preferably about 10K for single-end and 1K for paired-end)
 *      by the total number of reads to obtain the subsampling fraction value.
 * 3. create a BED file with the single entry `chr1\t1\t10000000` and name it `chr1.bed`
 * 4. run `samtools view -b -s <subsampling fraction> -L chr1.bed file.bam > test_file.bam`.
 */
class ReadsQueryTest {

    private val TO = GenomeQuery("to1")

    @Test
    fun testStemGz() {
        withTempFile("foo", ".bam") {
            assertFalse(it.stemGz.endsWith(".bam"))
        }
        withTempFile("foo", ".bed.gz") {
            assertFalse(it.stemGz.endsWith(".bed"))
        }
        withTempFile("foo", ".BED.GZ") {
            assertFalse(it.stemGz.toLowerCase().endsWith(".bed"))
        }
    }

    @Test
    fun testCoverageId() {
        withTempDirectory("foo") {
            val treatment = (it / "OD7_k4me1_hg19.bed.gz").apply { touch() }
            assertEquals("OD7_k4me1_hg19_unique_150", ReadsQuery(TO, treatment, true, 150).id)
            assertEquals("OD7_k4me1_hg19", ReadsQuery(TO, treatment, false, null).id)
        }
    }

    @Test
    fun testLoadSingleEndBam() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val coverage = readsQuery.get()
            val chr1 = genomeQuery["chr1"]!!
            assertEquals(SingleEndCoverage::class.java, coverage::class.java)
            assertEquals(136, (coverage as SingleEndCoverage).detectedFragment)
            assertEquals(8125, coverage.getBothStrandsCoverage(chr1.range.on(chr1)))
        }
    }

    @Test
    fun testLoadPairedEndBam() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val coverage = readsQuery.get()
            val chr1 = genomeQuery["chr1"]!!
            assertEquals(PairedEndCoverage::class.java, coverage::class.java)
            assertEquals(693, coverage.getBothStrandsCoverage(chr1.range.on(chr1)))
        }
    }

    @Test
    fun testSingleEndLoggingDefault() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val (out, err) = LogsTest.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: single_end.bam", out)
            assertIn("Depth: 8125", out)
            assertIn("Reads: single-ended", out)
            assertIn("Fragment size: 136 bp (cross-correlation estimate)", out)
            assertEquals("", err)
        }
    }

    @Test
    fun testSingleEndLoggingOverride() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false, fragment = 100)
            val (out, err) = LogsTest.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: single_end.bam", out)
            assertIn("Depth: 8125", out)
            assertIn("Reads: single-ended", out)
            assertIn("Fragment size: 100 bp (overrides cross-correlation estimate 136)", out)
            assertEquals("", err)
        }
    }

    @Test
    fun testPairedEndLoggingDefault() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val (out, err) = LogsTest.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: paired_end.bam", out)
            assertIn("Depth: 693", out)
            assertIn("Reads: paired-ended", out)
            assertIn("Fragment size: 139 bp (average; inferred from read pairs)", out)
            assertEquals("", err)
        }
    }

    @Test
    fun testPairedEndLoggingOverride() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            val readsQuery = ReadsQuery(genomeQuery, path, false, fragment = 100)
            val (out, err) = LogsTest.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: paired_end.bam", out)
            assertIn("Depth: 693", out)
            assertIn("Reads: paired-ended", out)
            assertIn(
                "Fragment size: 139 bp (average; inferred from read pairs; " +
                        "user input 100 is ignored)",
                out
            )
            assertEquals("", err)
        }
    }

    @Test
    fun testLoadSingleEndBamAsPairedEnd() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery("to1")
            var matePairs = 0
            val unpairedReads = processPairedReads(genomeQuery, path) { _, _ -> matePairs++ }
            assertEquals(0, matePairs)
            assertEquals(8125, unpairedReads)
        }
    }

}