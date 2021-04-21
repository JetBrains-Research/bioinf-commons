package org.jetbrains.bio.genome

import org.jetbrains.bio.util.*
import org.junit.Test
import java.io.ByteArrayOutputStream
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals
import kotlin.test.assertSame
import kotlin.test.assertTrue

class GenomeTest {
    @Test
    fun equality() {
        assertEquals(Genome["to1"], Genome["to1"])
    }

    @Test
    fun presentableName() {
        assertEquals("Test Organism: to1 [test]", Genome["to1"].presentableName())
    }

    @Test
    fun get() {
        withTempDirectory("foo") { dir ->
            val chromSizesPath1 = dir / "foo1.chrom.sizes"
            chromSizesPath1.bufferedWriter().use {
                it.write("fooBarBaz\t10\n")
                it.write("unknown\t10\n")
            }
            val chromSizesPath2 = dir / "foo2.chrom.sizes"
            chromSizesPath2.bufferedWriter().use {
                it.write("chr1\t100000\n")
                it.write("chr2\t100000\n")
                it.write("chr3\t100000\n")
            }
            assertSame(Genome[chromSizesPath1, "to2"], Genome[chromSizesPath1, "to2"])
            assertEquals(Genome[chromSizesPath1, "to2"], Genome[chromSizesPath1, "to2"])
            assertNotEquals(Genome[chromSizesPath1, "to3"], Genome[chromSizesPath1, "to4"])

            // At the moment we check only build & chrom size path
            assertSame(
                    Genome.get(chromSizesPath1, "to2", genesDescriptionsPath = "vers1".toPath()),
                    Genome.get(chromSizesPath1, "to2", genesDescriptionsPath = "vers2".toPath())
            )
        }
    }

    @Test(expected = IllegalArgumentException::class)
    fun conflictingChromSizePath() {
        withTempDirectory("foo") { dir ->
            val chromSizesPath1 = dir / "foo1.chrom.sizes"
            chromSizesPath1.bufferedWriter().use {
                it.write("fooBarBaz\t10\n")
                it.write("unknown\t10\n")
            }
            val chromSizesPath2 = dir / "foo2.chrom.sizes"
            chromSizesPath2.bufferedWriter().use {
                it.write("chr1\t100000\n")
                it.write("chr2\t100000\n")
                it.write("chr3\t100000\n")
            }

            assertNotEquals(
                    Genome[chromSizesPath1, "to2"],
                    Genome[chromSizesPath2, "to2"]
            )
        }
    }

    @Test
    fun chromSizes() {
        withTempDirectory("foo") { dir ->
            val chromSizesPath = dir / "to1.chrom.sizes"
            chromSizesPath.bufferedWriter().use {
                it.write("chr1\t100000\n")
                it.write("chr10\t100000\n")
                it.write("chr100\t100000\n")
            }
            val g0 = Genome["to1"]
            assertEquals("chr1, chr2, chr3, chrX, chrM", g0.chromosomes.joinToString { it.name })

            val g1 = Genome[chromSizesPath, "to1.${chromSizesPath.name}"]
            assertEquals("chr1, chr10, chr100", g1.chromosomes.joinToString { it.name })
            assertNotEquals(g0, g1)
        }
    }

    @Test
    fun chromSizes_LoggedError() {
        withTempFile("foo", ".galaxy.dat") { path ->
            val logStream = ByteArrayOutputStream()

            //XXX: Add log appender to catch warning about suspicious genome name
            addAppender(logStream, "AppenderGenomeTest")

            try {
                val genome = Genome[path]
                val build = path.fileName.toString().substringBefore(".")
                assertEquals(genome.build, build)
                assertTrue("Unexpected chrom sizes file name: ${path.fileName}, " +
                        "expected <build>.chrom.sizes. Detected build: $build" in
                        logStream.toString())
            } finally {
                Logs.getRootLogger().detachAppender("AppenderGenomeTest")
            }
        }
    }

    @Test
    fun chromosomeNamesToAltNamesMap() {
        val genome = Genome.getCustomised("to1_gt_chromosomeNamesToAltNamesMap", "to1",
            customized = true,
            forceUpdateCache = true
        ) {
            it.copy(chrAltName2CanonicalMapping = mapOf(
                "foo" to "chr1", "foo2" to "chr1", "boo" to "chr3"
            ))
        }
        require(genome.chrAltName2CanonicalMapping.size == 3)

        val chromosomeNamesToAltNamesMap = genome.chromosomeNamesToAltNamesMap
        assertEquals(
            chromosomeNamesToAltNamesMap["chr1"]!!.sorted(),
            listOf("1",  "foo", "foo2")
        )
        assertEquals(
            chromosomeNamesToAltNamesMap["chr3"]!!.sorted(),
            listOf("3", "boo")
        )
        assertEquals(
            chromosomeNamesToAltNamesMap["chr2"]!!.sorted(),
            listOf("2")
        )
    }

    @Test
    fun chromosomeNamesMap() {
        // Use uniq name here:
        val genome = Genome.getCustomised("to1_gt_chromosomeNamesMap", "to1",
            customized = true,
            forceUpdateCache = true
        ) {
            it.copy(chrAltName2CanonicalMapping = mapOf(
                "foo" to "chr1", "foo2" to "chr1", "boo" to "chr3"
            ))
        }
        require(genome.chrAltName2CanonicalMapping.size == 3)

        val chrNamesMap = genome.chromosomeNamesMap
        assertEquals("chr1", chrNamesMap["1"]!!.name)
        assertEquals("chr2", chrNamesMap["2"]!!.name)
        assertEquals("chr1", chrNamesMap["foo"]!!.name)
        assertEquals("chr1", chrNamesMap["foo2"]!!.name)
        assertEquals("chr3", chrNamesMap["boo"]!!.name)
        assertEquals("chr1", chrNamesMap["foo2"]!!.name)
    }
}