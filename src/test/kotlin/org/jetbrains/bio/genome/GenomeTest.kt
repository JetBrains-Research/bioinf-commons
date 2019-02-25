package org.jetbrains.bio.genome

import org.apache.log4j.SimpleLayout
import org.apache.log4j.WriterAppender
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempDirectory
import org.junit.Test
import java.io.ByteArrayOutputStream
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class GenomeTest {
    @Test
    fun equality() {
        assertEquals(Genome["to1"], Genome["to1"])
    }

    @Test
    fun presentableName() {
        assertEquals("Test Organism: to1", Genome["to1"].presentableName())
    }

    @Test
    fun testIgnoredChromosomes() {
        val logContent = ByteArrayOutputStream()
        val appender = WriterAppender(SimpleLayout(), logContent).apply { name = "test appender" }
        Genome.LOG.addAppender(appender)
        try {
            withTempDirectory("foo") { dir ->
                val chromSizesPath = dir / "foo.chrom.sizes"
                chromSizesPath.bufferedWriter().use {
                    it.write("fooBarBaz\t10\n")
                    it.write("unknown\t10\n")
                    it.write("chr1\t100000\n")
                    it.write("chr2\t100000\n")
                    it.write("chr3\t100000\n")
                }
                val genome = Genome["to1", chromSizesPath]
                genome.toQuery().only(listOf("chr1"))
                val log = logContent.toString()
                assertTrue("Ignored chromosomes $chromSizesPath: fooBarBaz, unknown" in log)
                assertEquals(1, "Ignored chromosomes".toRegex().findAll(log).toList().size)
            }
        } finally {
            Genome.LOG.removeAppender(appender)
        }
    }

}