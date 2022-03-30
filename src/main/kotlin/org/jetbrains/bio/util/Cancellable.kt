package org.jetbrains.bio.util

import com.google.common.annotations.VisibleForTesting
import com.google.common.base.MoreObjects
import org.slf4j.LoggerFactory
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger
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

/**
 * A Cancellable Task.
 *
 * @author Oleg Shpynov
 * @since 11/07/14
 */
class CancellableTask<T>(private val callable: Callable<T>) {
    /** Process-unique task ID. */
    val id = TASK_COUNTER.incrementAndGet()

    @Volatile
    var cancelled = false
        private set

    @Volatile
    private var cancellableState: CancellableState? = null

    @Volatile
    var task: Future<T>? = null
        internal set

    val isDone: Boolean get() = task != null && task!!.isDone

    fun cancel() {
        if (LOG.isTraceEnabled) {
            LOG.trace("Cancelled $id from ${Thread.currentThread()}")
        }
        cancelled = true
        if (cancellableState != null) {
            cancellableState!!.cancel()
        }

        if (task != null && !task!!.isDone) {
            task!!.cancel(true)
        }
    }

    fun execute() {
        if (cancelled) {
            return
        }
        if (LOG.isTraceEnabled) {
            LOG.trace("Executed $id from ${Thread.currentThread()}")
        }
        task = EXECUTOR.submit<T> {
            if (cancelled) {
                // Do not start task if it is marked as cancelled.
                return@submit null
            }

            cancellableState = CancellableState.current().apply { reset() }
            callable.call()
        }
    }

    /**
     * @return result of callable.
     * @throws CancellationException in case when process was cancelled, see [.cancel]
     */
    @Throws(CancellationException::class)
    fun get(): T {
        if (cancelled) {
            throw CancellationException()
        }
        check(task != null) { "Task not started: $id" }
        check(task!!.isDone) { "Task not ready: $id" }
        try {
            return task!!.get()
        } catch (e: InterruptedException) {
            throw CancellationException(e.message)
        } catch (e: CancellationException) {
            throw e
        } catch (e: ExecutionException) {
            // Process inner task exceptions
            val cause = e.cause
            if (cause is CancellationException) {
                throw cause
            }

            throw RuntimeException(e)
        }
    }

    fun waitAndGet(waitMillis: Long = 1000L): T? {
        while (true) {
            if (cancelled) {
                if (LOG.isTraceEnabled) {
                    LOG.trace("Cancelled $id")
                }
                return null
            }
            if (isDone) {
                if (LOG.isTraceEnabled) {
                    LOG.trace("Loaded task $id")
                }
                return get()
            }

            try {
                Thread.sleep(waitMillis)
            } catch (ignored: InterruptedException) {
                // ignored
            }
        }
    }

    override fun toString() = MoreObjects.toStringHelper(this).addValue(id).toString()

    companion object {
        private val LOG = LoggerFactory.getLogger(CancellableTask::class.java)

        private val EXECUTOR = Executors.newWorkStealingPool(parallelismLevel())

        private val TASK_COUNTER = AtomicInteger(0)

        @VisibleForTesting
        fun resetCounter() = TASK_COUNTER.set(0)
    }
}
