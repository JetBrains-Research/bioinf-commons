package org.jetbrains.bio.query

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertIs
import org.jetbrains.bio.coverage.PairedEndCoverage
import org.jetbrains.bio.coverage.SingleEndCoverage
import org.jetbrains.bio.coverage.processPairedReads
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.*
import org.junit.Test
import java.util.*
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

    private val TO = GenomeQuery(Genome["to1"])

    private val SINGLE_END_BAM_READS = 8125
    private val SINGLE_END_BAM_DETECTED_FRAGMENT = 136
    private val PAIRED_END_BAM_PAIRS = 693
    private val PAIRED_END_BAM_INFERRED_FRAGMENT = 139

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
            assertEquals("OD7_k4me1_hg19_unique_150", ReadsQuery(TO, treatment, true, Optional.of(150)).id)
            assertEquals("OD7_k4me1_hg19", ReadsQuery(TO, treatment, false, Optional.empty()).id)
        }
    }

    @Test
    fun testLoadSingleEndBam() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val coverage = readsQuery.get()
            val chr1 = genomeQuery["chr1"]!!
            assertEquals(SingleEndCoverage::class.java, coverage::class.java)
            assertEquals(SINGLE_END_BAM_DETECTED_FRAGMENT, (coverage as SingleEndCoverage).detectedFragment)
            assertEquals(SINGLE_END_BAM_READS, coverage.getBothStrandsCoverage(chr1.range.on(chr1)))
        }
    }

    @Test
    fun testLoadPairedEndBam() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val coverage = readsQuery.get()
            val chr1 = genomeQuery["chr1"]!!
            assertEquals(PairedEndCoverage::class.java, coverage::class.java)
            assertEquals(PAIRED_END_BAM_PAIRS, coverage.getBothStrandsCoverage(chr1.range.on(chr1)))
        }
    }

    @Test
    fun testSingleEndLoggingDefault() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val (out, err) = Logs.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: single_end.bam", out)
            assertIn("Depth: $SINGLE_END_BAM_READS", out)
            assertIn("Reads: single-ended", out)
            assertIn(
                "Fragment size: $SINGLE_END_BAM_DETECTED_FRAGMENT bp (cross-correlation estimate)",
                out
            )
            assertEquals("", err)
        }
    }

    @Test
    fun testSingleEndLoggingOverride() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false, fragment = Optional.of(100))
            val (out, err) = Logs.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: single_end.bam", out)
            assertIn("Depth: $SINGLE_END_BAM_READS", out)
            assertIn("Reads: single-ended", out)
            assertIn(
                "Fragment size: 100 bp (overrides cross-correlation estimate $SINGLE_END_BAM_DETECTED_FRAGMENT)",
                out
            )
            assertEquals("", err)
        }
    }

    @Test
    fun testPairedEndLoggingDefault() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false)
            val (out, err) = Logs.captureLoggingOutput { readsQuery.get() }
            assertIn("Library: paired_end.bam", out)
            assertIn("Depth: $PAIRED_END_BAM_PAIRS", out)
            assertIn("Reads: paired-ended", out)
            assertIn(
                "Fragment size: $PAIRED_END_BAM_INFERRED_FRAGMENT bp (average; inferred from read pairs)",
                out
            )
            assertEquals("", err)
        }
    }

    /**
     * Should be impossible to achieve via normal [ReadsQuery] invocation, so we use
     * internal [processPairedReads] method.
     */
    @Test
    fun testLoadSingleEndBamAsPairedEnd() {
        withResource(ReadsQueryTest::class.java, "single_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            var pairedReads = 0
            val unpairedReads = processPairedReads(genomeQuery, path) { _, _, _, _ -> pairedReads++ }
            assertEquals(0, pairedReads)
            assertEquals(SINGLE_END_BAM_READS, unpairedReads)
        }
    }

    /**
     * fragment=0 option should force loading as single-end
     */
    @Test
    fun testLoadPairedEndBamAsSingleEnd() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false, fragment = Optional.of(0))
            val coverage = readsQuery.get()
            val chr1 = genomeQuery["chr1"]!!
            assertIs(coverage, SingleEndCoverage::class.java)
            assertEquals(
                PAIRED_END_BAM_PAIRS * 2,
                coverage.getBothStrandsCoverage(chr1.range.on(chr1))
            )
        }
    }

    /**
     * fragment=0 option should force loading as single-end
     */
    @Test
    fun testLoadPairedEndBamAsSingleEndLogging() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQuery = ReadsQuery(genomeQuery, path, false, fragment = Optional.of(0))
            val (out, err) = Logs.captureLoggingOutput { readsQuery.get() }
            assertIn("Fragment option (0) forces reading paired-end reads as single-end!", out)
            assertEquals("", err)
        }
    }

    /**
     * We've had troubles with cache file reuse (see issue #1). This test checks that
     * the cache file is not reused when not appropriate.
     */
    @Test
    fun testLoadPairedEndBamAsPairedEndThenSingleEnd() {
        withResource(ReadsQueryTest::class.java, "paired_end.bam") { path ->
            val genomeQuery = GenomeQuery(Genome["to1"])
            val readsQueryAsPaired = ReadsQuery(genomeQuery, path, false)
            val coverageAsPaired = readsQueryAsPaired.get()
            val readsQueryAsSingle = ReadsQuery(genomeQuery, path, false, fragment = Optional.of(0))
            val coverageAsSingle = readsQueryAsSingle.get()
            assertIs(coverageAsPaired, PairedEndCoverage::class.java)
            assertIs(coverageAsSingle, SingleEndCoverage::class.java)
        }
    }

}