package org.jetbrains.bio.genome

import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempDirectory
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals
import kotlin.test.assertSame

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
    fun testCompare() {
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
            assertSame(Genome["to1", chromSizesPath1], Genome["to1", chromSizesPath1])
            assertEquals(Genome["to1", chromSizesPath1], Genome["to1", chromSizesPath1])
            assertNotEquals(Genome["to1", chromSizesPath1], Genome["to1", chromSizesPath2])
            assertNotEquals(Genome["to1", chromSizesPath1], Genome["to2", chromSizesPath1])
            assertEquals(-1, Genome["to1", chromSizesPath1].compareTo(Genome["to2", chromSizesPath1]))
            assertEquals(chromSizesPath1.compareTo(chromSizesPath2),
                    Genome["to1", chromSizesPath1].compareTo(Genome["to1", chromSizesPath2]))
        }
    }

}