@file:Suppress("unused")

package org.jetbrains.bio.util

import java.io.InputStream
import java.io.Reader
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

/**
 * @author Roman.Chernyatchik
 * @date 2018-10-24
 */

object IOMonitor {
    // Used for debug
    // internal val LOG = LoggerFactory.getLogger(IOMonitor::class.java)

    @Volatile
    var debugMode: Boolean = false

    val internalCntOpenedStreams: AtomicInteger = AtomicInteger()
    val openedStreams: Int
        get() = internalCntOpenedStreams.get()

    private val internalCntMaxOpenedStreams: AtomicInteger = AtomicInteger()
    val maxOpenedStreams: Int
        get() = internalCntMaxOpenedStreams.get()

    private val internalCntTotalOpenedStreams: AtomicInteger = AtomicInteger()
    val totalOpenedStreams: Int
        get() = internalCntTotalOpenedStreams.get()

    val internalCntBytesRead: AtomicLong = AtomicLong()
    val bytesRead: Long
        get() = internalCntBytesRead.get()

    fun openStreamOrReader() {
        internalCntOpenedStreams.getAndIncrement()
        internalCntMaxOpenedStreams.set(maxOf(internalCntOpenedStreams.get(), internalCntMaxOpenedStreams.get()))
        internalCntTotalOpenedStreams.getAndIncrement()
    }
}

class IOMonitorReader(val input: Reader) : Reader() {
    @Volatile
    var closed = false

    init {
        IOMonitor.openStreamOrReader()
    }

    override fun close() {
        if (closed) {
            return
        }
        synchronized(this) {
            input.close()
            closed = true
        }
        IOMonitor.internalCntOpenedStreams.getAndDecrement()
    }

    override fun read(cbuf: CharArray?, off: Int, len: Int): Int {
        val n = input.read(cbuf, off, len)

        IOMonitor.internalCntBytesRead.addAndGet(maxOf(0, n).toLong())
        return n
    }
}

class IOMonitorInputStream(val input: InputStream) : InputStream() {
    @Volatile
    var closed = false

    init {
        IOMonitor.openStreamOrReader()
    }

    override fun close() {
        if (closed) {
            return
        }
        synchronized(this) {
            input.close()
            closed = true
        }
        IOMonitor.internalCntOpenedStreams.getAndDecrement()
    }

    override fun read(buffer: ByteArray, offset: Int, length: Int): Int {
        val n = input.read(buffer, offset, length)

        IOMonitor.internalCntBytesRead.addAndGet(maxOf(0, n).toLong())

        return n
    }

    override fun read(): Int {
        val value = input.read()

        if (value != -1) {
            IOMonitor.internalCntBytesRead.addAndGet(1)
        }

        return value
    }
}