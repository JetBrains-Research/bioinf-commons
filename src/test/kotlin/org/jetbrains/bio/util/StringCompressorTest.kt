package org.jetbrains.bio.util

import org.jetbrains.bio.Tests.assertThrowsWithMessage
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class StringCompressorTest {
    @Test
    fun compressDecompressSingleString() {
        val compressor = StringCompressor()
        val s = "foobar"
        assertEquals(s, compressor.compress(s).toString())
    }

    @Test
    fun compressDecompress100KStrings() {
        val compressor = StringCompressor()
        val compressedStrings = ArrayList<CompressedString>()
        val n = 100000
        (0 until n).mapTo(compressedStrings) { compressor.compress("foo_${it}_bar") }
        (0 until n).forEach { assertEquals("foo_${it}_bar", compressedStrings[it].toString()) }
    }

    @Test
    fun addStringToFinishedChunk() {
        val chunk = StringCompressorChunk()
        chunk.compress("foo")
        chunk.finish()
        assertTrue(chunk.finished, "Chunk should be finished after finish() invocation.")

        assertThrowsWithMessage(
            IllegalStateException::class.java,
            "finished StringCompressorChunk",
            partialMessageMatch = true
        ) {
            chunk.compress("bar")
        }
    }


    @Test
    fun autoFinishChunk() {
        val chunk = StringCompressorChunk()
        val longString = "a".repeat(StringCompressorChunk.FLUSH_THRESHOLD)
        chunk.compress(longString)
        assertTrue(chunk.finished, "Chunk should be finished after reaching flush threshold.")

        assertThrowsWithMessage(
            IllegalStateException::class.java,
            "finished StringCompressorChunk",
            partialMessageMatch = true
        ) {
            chunk.compress("foo")
        }
    }
}