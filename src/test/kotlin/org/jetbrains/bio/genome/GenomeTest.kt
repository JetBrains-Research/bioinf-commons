package org.jetbrains.bio.genome

import org.apache.log4j.SimpleLayout
import org.apache.log4j.WriterAppender
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
    fun testGet() {
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
    fun testConflictingChromSizePath() {
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
    fun testChromSizes() {
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
    fun testChromSizes_LoggedError() {
        withTempFile("foo", ".galaxy.dat") { path ->
            val logContent = ByteArrayOutputStream()

            //XXX: Add log appender to catch warning about suspicious genome name
            val appender = WriterAppender(SimpleLayout(), logContent).apply { name = "test appender" }
            Genome.LOG.addAppender(appender)
            
            try {
                val genome = Genome[path]
                val build = path.fileName.toString().substringBefore(".")
                assertEquals(genome.build, build)
                assertTrue("Unexpected chrom sizes file name: ${path.fileName}, " +
                        "expected <build>.chrom.sizes. Detected build: $build" in
                        logContent.toString())
            } finally {
                MultitaskProgress.LOG.removeAppender(appender)
            }
        }
    }
}