package org.jetbrains.bio.util

import com.google.common.annotations.VisibleForTesting
import com.google.common.base.MoreObjects
import org.slf4j.LoggerFactory
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger

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
        LOG.trace("Cancelled $id from ${Thread.currentThread()}")
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
        LOG.trace("Executed $id from ${Thread.currentThread()}")
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
                LOG.trace("Cancelled $id")
                return null
            }
            if (isDone) {
                LOG.trace("Loaded task $id")
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
