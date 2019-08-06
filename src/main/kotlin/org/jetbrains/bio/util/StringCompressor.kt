package org.jetbrains.bio.util

import kotlinx.support.jdk7.use
import java.io.ByteArrayOutputStream
import java.lang.ref.WeakReference
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

class StringCompressor {

    private var currentChunk: StringCompressorChunk = StringCompressorChunk()

    fun compress(s: String): CompressedString {
        synchronized(this) {
            if (currentChunk.finished) currentChunk = StringCompressorChunk()
            return currentChunk.compress(s)
        }
    }

    fun finish() {
        currentChunk.finish()
    }
}

class StringCompressorChunk {

    private var originalStringBuilder: StringBuilder? = StringBuilder(2 * FLUSH_THRESHOLD)
    private var decompressedLength: Int = 0

    private lateinit var compressedString: ByteArray
    private var weakDecompressedString: WeakReference<ByteArray> = WeakReference<ByteArray>(null)

    val decompressedString: ByteArray
        get() {
            if (!finished) {
                synchronized(this) { finish() }
            }
            val stored = weakDecompressedString.get()
            if (stored != null) return stored
            synchronized(this) {
                // check whether another thread already computed the value
                val stored2 = weakDecompressedString.get()
                if (stored2 != null) return stored2
                val computed = GZIPInputStream(compressedString.inputStream()).readBytes()
                check(computed.size == decompressedLength) { "Decompressed length doesn't match the message" }
                weakDecompressedString = WeakReference(computed)
                return computed
            }
        }

    val finished get() = originalStringBuilder == null

    fun compress(s: String): CompressedString {
        check(!finished) { "Attempting to add a String to a finished StringCompressorChunk" }
        originalStringBuilder!!.append(s)
        val offset = decompressedLength
        decompressedLength += s.length
        if (decompressedLength >= FLUSH_THRESHOLD) {
            finish()
        }
        return CompressedString(this, offset, s.length)
    }

    fun compressionRatio(): Double {
        if (!finished) return Double.NaN
        return compressedString.size * 1.0 / decompressedLength
    }

    fun finish() {
        if (finished) return
        val originalString = originalStringBuilder!!.toString().toByteArray()
        val outputStream = ByteArrayOutputStream()
        GZIPOutputStream(outputStream).use { it.write(originalString) }
        compressedString = outputStream.toByteArray()
        weakDecompressedString = WeakReference(originalString)
        originalStringBuilder = null
    }

    companion object {
        const val FLUSH_THRESHOLD = 65536
    }
}

data class CompressedString(
        private val chunk: StringCompressorChunk,
        private val offset: Int,
        private val length: Int
) {
    override fun toString() = String(chunk.decompressedString, offset, length)
}

fun main() {
    val compressor = StringCompressor()
    val list = ArrayList<CompressedString>()
    (0 until 100000).mapTo(list) { i -> compressor.compress("terribly_long_string_not_unlike_macs2_peak_name_$i") }
    compressor.finish()
    check(list.mapIndexed { index, compressedString -> compressedString.toString() == "terribly_long_string_not_unlike_macs2_peak_name_$index" }.all { it })
}