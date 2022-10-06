package org.jetbrains.bio.util

import kotlinx.support.jdk7.use
import org.jetbrains.bio.util.StringCompressorChunk.Companion.FLUSH_THRESHOLD
import java.io.ByteArrayOutputStream
import java.lang.ref.WeakReference
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

/**
 * A GZIP-based multi-string compressor.
 *
 * Designed to effectively store a large number of short and repetitive strings. Thread-safe.
 * [compress] returns a compressed representation of the provided string.
 *
 * Example:
 *
 *     val compressor = StringCompressor()
 *     val compressedStringList = ArrayList<CompressedString>()
 *     stringList.mapTo(compressedStringList) { compressor.compress(it) }
 *
 */
class StringCompressor {

    private var currentChunk: StringCompressorChunk = StringCompressorChunk()

    /**
     * Returns a compressed representation of the provided string.
     *
     * Use [CompressedString.toString] to obtain the original string. Thread-safe.
     */
    fun compress(s: String): CompressedString {
        synchronized(this) {
            if (currentChunk.finished) {
                currentChunk = StringCompressorChunk()
            }
            return currentChunk.compress(s)
        }
    }
}

/**
 * An accumulator chunk used by [StringCompressor] to compress strings.
 *
 * When not [finished], accepts strings via [compress] method. The strings are appended to a [StringBuilder].
 * Once the builder size exceeds [FLUSH_THRESHOLD], the chunk becomes [finished], compresses the builder
 * contents and releases the builder reference. A [finished] chunk doesn't accept any more strings to [compress].
 * It does, however, provide a [decompressedString] property to access the original UTF-8-encoded contents.
 */
class StringCompressorChunk {

    private var originalStringBuilder: StringBuilder? = StringBuilder(2 * FLUSH_THRESHOLD)
    private var decompressedLength: Int = 0

    private lateinit var compressedString: ByteArray
    private var weakDecompressedString: WeakReference<ByteArray> = WeakReference<ByteArray>(null)

    /**
     * A UTF-8 encoded [ByteArray] representation of the original concatenated string. Thread-safe.
     *
     * It's cached via a [WeakReference] and recalculated as necessary.
     * Accessing this property will [finish] the chunk if it's not [finished] already.
     */
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

    /**
     * If true, the chunk has been compressed and no longer accepts strings to [compress].
     */
    val finished get() = originalStringBuilder == null

    /**
     * Add a string to the chunk and return its compressed representation. Not thread-safe.
     *
     * If the chunk length exceeds [FLUSH_THRESHOLD] after the addition, [finish] the chunk.
     */
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

    /**
     * Compress the chunk contents and mark it as [finished].
     *
     * The chunk will no longer accept new strings to [compress]. The method is idempotent: calling it
     * repeatedly has the same effect as calling it once. Not thread-safe.
     */
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

/**
 * A compressed representation of a string.
 *
 * Use [toString] to get the original string. This implicitly causes the appropriate chuck
 * to [StringCompressorChunk.finish], if it hasn't done so already.
 */
data class CompressedString internal constructor(
    private val chunk: StringCompressorChunk,
    private val offset: Int,
    private val length: Int
) {
    override fun toString() = String(chunk.decompressedString, offset, length)
}