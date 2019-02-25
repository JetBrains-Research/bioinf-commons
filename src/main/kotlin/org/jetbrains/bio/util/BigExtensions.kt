package org.jetbrains.bio.util

import htsjdk.samtools.seekablestream.SeekableStream
import htsjdk.samtools.seekablestream.SeekableStreamFactory
import org.apache.log4j.Logger.getLogger
import org.jetbrains.bio.BetterSeekableBufferedStream
import org.jetbrains.bio.EndianAwareDataSeekableStream
import org.jetbrains.bio.EndianSynchronizedBufferFactory
import org.jetbrains.bio.RomBufferFactory
import org.jetbrains.bio.big.BigBedFile
import org.jetbrains.bio.big.BigFile
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.tdf.TdfFile
import org.jetbrains.bio.util.IOUsageCounter.debugMode
import org.jetbrains.bio.util.IOUsageCounter.internalCntBytesRead
import org.jetbrains.bio.util.IOUsageCounter.internalCntOpenedStreams
import java.net.URI
import java.nio.ByteOrder

/**
 * @author Roman.Chernyatchik
 */

/**
 * Opens customized BigWig file which supports parallel query operations,
 * I/O usage monitoring and works both with local and remote files.
 * Is required for JBR browser.
 *
 * @source Input path or
 * @prefetch Amount of internal indexes which will be loaded to memory
 *      If you load more indexes queries will be performed faster but file
 *      opening will be slower
 * @cancelledChecker Allows to stop current operation by throwing Exception
 *      E.g. throw cancelled exception here
 * @return BigWigFile object
 */
fun BigWigFile.Companion.readURIParallelAccess(
        source: URI,
        prefetch: Int = BigUtil.defaultPrefetchLevel,
        cancelledChecker: (() -> Unit)?
) = read(
        source.presentablePath(), prefetch,
        cancelledChecker = cancelledChecker,
        factoryProvider = BigUtil.factory
)

/**
 * Opens customized BigBed file which supports parallel query operations,
 * I/O usage monitoring and works both with local and remote files.
 * Is required for JBR browser.
 *
 * @source Input path or
 * @prefetch Amount of internal indexes which will be loaded to memory
 *      If you load more indexes queries will be performed faster but file
 *      opening will be slower
 * @cancelledChecker Allows to stop current operation by throwing Exception
 *      E.g. throw cancelled exception here
 * @return BigBedFile object
 */
fun BigBedFile.Companion.readURIParallelAccess(
        source: URI,
        prefetch: Int = BigUtil.defaultPrefetchLevel,
        cancelledChecker: (() -> Unit)?
) = read(
        source.presentablePath(), prefetch,
        cancelledChecker = cancelledChecker,
        factoryProvider = BigUtil.factory
)

/**
 * Opens customized Tdf file which supports parallel query operations,
 * I/O usage monitoring and works both with local and remote files.
 * Is required for JBR browser.
 *
 * @source Input path or
 * @prefetch Amount of internal indexes which will be loaded to memory
 *      If you load more indexes queries will be performed faster but file
 *      opening will be slower
 * @cancelledChecker Allows to stop current operation by throwing Exception
 *      E.g. throw cancelled exception here
 * @return TdfFile object
 */

fun TdfFile.Companion.readURIParallelAccess(source: URI, cancelledChecker: (() -> Unit)? = null) = read(
        source.presentablePath(),
        cancelledChecker = cancelledChecker,
        factoryProvider = BigUtil.factory
)

/**
 * Opens customized BigWig or BigBed file which supports parallel query operations,
 * I/O usage monitoring and works both with local and remote files.
 * Is required for JBR browser.
 *
 * If you know the exact type of a file, it is recommend to use specialized implementation:
 * [org.jetbrains.bio.big.BigBedFile.readURIParallelAccess] or
 * [org.jetbrains.bio.big.BigBedFile.readURIParallelAccess].
 * Specialized implementations initialize files slightly faster.
 *
 * @source Input path or
 * @prefetch Amount of internal indexes which will be loaded to memory
 *      If you load more indexes queries will be performed faster but file
 *      opening will be slower
 * @cancelledChecker Allows to stop current operation by throwing Exception
 *      E.g. throw cancelled exception here
 * @return BigFile object
 */
fun BigFile.Companion.readURIParallelAccess(
        source: URI,
        prefetch: Int = BigUtil.defaultPrefetchLevel,
        cancelledChecker: (() -> Unit)? = null
) = read(
        source.presentablePath(), prefetch,
        cancelledChecker = cancelledChecker,
        factoryProvider = BigUtil.factory
)

object BigUtil {
    private val LOG = getLogger(BigUtil::class.java)

    private const val INTERNAL_RMF_KEY = "epigenome.internal.big.rmf"
    private const val INTERNAL_BIG_PREFETCH_KEY = "epigenome.internal.big.prefetch"

    val factory: (String, ByteOrder) -> RomBufferFactory

    init {
        val defaultSize = 128000
        val romMemBufferProperty = System.getProperty(INTERNAL_RMF_KEY, "$defaultSize")
        val bufferSize = try {
            romMemBufferProperty.toInt()
        } catch (e: Exception) {
            LOG.error("Cannot parse $INTERNAL_RMF_KEY value: $romMemBufferProperty, using default " +
                    "buffer size $defaultSize bytes")
            defaultSize
        }

        // Customized big buffer factory which supports remote files, not only local ones
        factory = { p, byteOrder ->
            val stream = EndianAwareDataSeekableStream(
                    object : BetterSeekableBufferedStream(
                            // Our custom impl for better FTP access support
                            //  org.jetbrains.bio.util.seekablestream.SeekableStreamFactory
                            SeekableStreamFactory.getInstance()
                                    .getStreamFor(p).let { stream ->
                                        // Do not add overhead if not required:
                                        when {
                                            debugMode -> IOUsageCountingSeekableStream(stream)
                                            else -> stream
                                        }
                                    },
                            bufferSize
                    ) {
                        override fun fillBuffer(): Int {
                            CancellableState.current().checkCanceled()
                            return super.fillBuffer()
                        }
                    }).apply {
                order = byteOrder
            }

            EndianSynchronizedBufferFactory(stream)

        }
    }

    val defaultPrefetchLevel: Int
        get() = System.getProperty(INTERNAL_BIG_PREFETCH_KEY, "${BigFile.PREFETCH_LEVEL_DETAILED}").toInt()
}


class IOUsageCountingSeekableStream(val input: SeekableStream) : SeekableStream() {
    @Volatile var closed = false
    
    init {
        IOUsageCounter.openStreamOrReader()
    }

    override fun length() = input.length()

    override fun getSource(): String? = input.source

    override fun eof() = input.eof()

    override fun seek(position: Long) = input.seek(position)

    override fun position() = input.position()

    override fun close() {
        if (closed) {
            return
        }
        synchronized(this) {
            input.close()
            closed = true
        }
        internalCntOpenedStreams.getAndDecrement()
    }

    override fun read(buffer: ByteArray, offset: Int, length: Int): Int {
        val n = input.read(buffer, offset, length)

        internalCntBytesRead.addAndGet(maxOf(0, n).toLong())

        return n
    }

    override fun read(): Int {
        val value = input.read()

        if (value != -1) {
            internalCntBytesRead.addAndGet(1)
        }

        return value
    }
}