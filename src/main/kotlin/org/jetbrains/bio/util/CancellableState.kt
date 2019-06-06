package org.jetbrains.bio.util

import java.util.concurrent.CancellationException
import kotlin.concurrent.getOrSet

/**
 * A thread-local cancellation flag.
 *
 * Contract:
 *
 * Suppose we have a thread which invokes several [CancellableTask],
 * each task has its own instance of [CancellableState], which can be cancelled from outer thread,
 * Task itself invokes [checkCanceled] which can cause exception, immediately interrupting it.
 *
 * @author Oleg Shpynov
 * @since 11/07/14
 */
class CancellableState private constructor() {
    // volatile for ThreadLocal is required here to propagate cancellation to another threads
    // Can be used in the scheme where singe instance is shared across different children threads
    @Volatile
    private var cancelled = false

    /**
     * Sets cancellation flag for the calling thread.
     */
    fun cancel() {
        cancelled = true
    }

    /**
     * Resets cancellation flag for the calling thread.
     */
    fun reset() {
        cancelled = false
    }

    /**
     * @throws [CancellationException] if cancellation flag is set.
     */
    fun checkCanceled() {
        if (cancelled) {
            throw CancellationException()
        }
    }

    companion object {
        private val THREAD_LOCALS = ThreadLocal<CancellableState>()

        fun current() = THREAD_LOCALS.getOrSet(::CancellableState)
    }
}
