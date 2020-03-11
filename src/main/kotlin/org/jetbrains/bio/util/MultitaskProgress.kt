package org.jetbrains.bio.util

import com.google.common.annotations.VisibleForTesting
import org.jetbrains.bio.util.MultitaskProgress.addTask
import org.jetbrains.bio.util.MultitaskProgress.finishTask
import org.jetbrains.bio.util.MultitaskProgress.periodNanos
import org.jetbrains.bio.util.MultitaskProgress.reportTask
import org.slf4j.LoggerFactory
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.atomic.LongAdder

/**
 * @author Alexey Dievsky
 * @author Oleg.Shpynov
 * @date 19.05.2015
 */

/**
 * This object tracks progress of multiple tasks. Each task should have a unique identifier (of Any type),
 * which it uses to register itself [addTask], report progress [reportTask] and finish tracking [finishTask].
 * These methods are thread-safe.
 *
 * Note that reporting and finishing tasks that haven't been added isn't an error.
 */
object MultitaskProgress {
    @VisibleForTesting
    internal val LOG = LoggerFactory.getLogger(MultitaskProgress::class.java)

    // Hack to show progress in UI
    var progressBar: ProgressBar? = null

    private var startTime: Long = 0
    private val periodNanos: Long = TimeUnit.SECONDS.toNanos(10)
    private val lastDuration: AtomicLong = AtomicLong(-periodNanos)
    private var progressPercent: Long = -1

    private val totalItemsByTasks: MutableMap<Any, Long> = ConcurrentHashMap()
    private val processedItemsByTasks: MutableMap<Any, LongAdder> = ConcurrentHashMap()
    private val totalItems: AtomicLong = AtomicLong()
    private val processedItems: LongAdder = LongAdder()

    private fun tell(s: String) {
        if (totalItemsByTasks.size > 1) {
            LOG.info("Running ${totalItemsByTasks.size} tasks: $s")
        } else {
            LOG.info(s)
        }
    }

    /**
     * Restarts the progress tracker. It is assumed that no tasks are running when this method is called.
     * Not thread-safe, requires external synchronization.
     */
    fun reset() {
        startTime = System.nanoTime()
        totalItems.set(0)
        processedItems.reset()
        progressPercent = -1
        lastDuration.set(-periodNanos)
    }

    /**
     * Drops all tracked tasks and restores the item counters. Not thread-safe. Don't use outside of tests.
     */
    fun clear() {
        totalItemsByTasks.clear()
        processedItemsByTasks.clear()
        totalItems.set(0L)
        processedItems.reset()
    }

    /**
     * Add a task to be tracked. Provide a reasonable estimate of iterations needed for the task to complete;
     * this estimate doesn't need to be the upper bound. If no tasks were running prior to the invocation,
     * the tracker is reset.
     * Tries to print progress by calling [reportIfNecessary].
     */
    @JvmStatic
    fun addTask(taskID: Any, items: Long, quite: Boolean = false) {
        synchronized(taskID) {
            if (totalItemsByTasks.isEmpty()) {
                reset()
            }
            if (totalItemsByTasks.containsKey(taskID)) {
                tell("Trying to add task $taskID that is already running")
                return
            }
            totalItemsByTasks[taskID] = items
            processedItemsByTasks[taskID] = LongAdder()
            totalItems.addAndGet(items)
            if (!quite) {
                reportIfNecessary(false)
            }
        }
    }

    /**
     * Finish tracking a task. The unused iterations reserved by the task are subtracted from the total pool.
     * Tries to print progress by calling [reportIfNecessary].
     */
    @JvmStatic
    fun finishTask(taskID: Any) {
        synchronized(taskID) {
            if (!totalItemsByTasks.containsKey(taskID)) {
                return
            }
            totalItems.addAndGet(processedItemsByTasks[taskID]!!.toLong() - totalItemsByTasks[taskID]!!)
            totalItemsByTasks.remove(taskID)
            processedItemsByTasks.remove(taskID)
        }
        reportIfNecessary(true)
    }

    /**
     * Report that a task completed an iteration. If this exceeds the allotted number of iterations,
     * the number is doubled. Tries to print progress by calling [reportIfNecessary].
     */
    @JvmStatic
    fun reportTask(taskID: Any, items: Long = 1L) {
        if (!totalItemsByTasks.containsKey(taskID)) {
            return
        }
        val processedItemsForTask = processedItemsByTasks[taskID]!!
        val totalItemsForTask = totalItemsByTasks[taskID]!!
        processedItemsForTask.add(items)
        processedItems.add(items)
        if (processedItemsForTask.toLong() > totalItemsForTask) {
            synchronized(taskID) {
                // Update total items if we are out of range
                if (processedItemsForTask.toLong() > totalItemsForTask) {
                    totalItemsByTasks[taskID] = totalItemsForTask * 2
                    totalItems.addAndGet(totalItemsForTask)
                }
            }
        }
        reportIfNecessary(false)
    }

    /**
     * Reports progress iff the two conditions are met:
     * 1. no less than [periodNanos] nanoseconds have elapsed since the last report
     * OR there were no reports since the last reset
     * 2. the progress value (percentage rounded to the integer) differs from the last reported one
     * OR there were no reports since the last reset and the progress value is not zero
     *
     */
    private fun reportIfNecessary(finishEvent: Boolean) {
        val duration = Progress.getDuration(startTime, lastDuration, periodNanos, finishEvent)
        if (duration == -1L && !finishEvent) {
            return
        }

        synchronized(this) { //synchronized(LOCK) {
            val processed = processedItems.toLong()
            val total = totalItems.get()
            val processedPercent = (processed * 100.0 / total).toLong()
            if (processedPercent != progressPercent || (finishEvent && processedPercent != 100L)) {
                progressPercent = processedPercent
                if (processed < total || finishEvent) {
                    progressBar?.setState(processed, total)

                    val progressPart = "%.2f%% (%,d/%,d), Elapsed time: ${asTime(duration)}"
                            .format(processed * 100.0 / total, processed, total)

                    val throughputPart = if (processed > 1) {
                        val eta = asTime((duration.toDouble() * total / processed).toLong() - duration)
                        val avgThroughput = asThroughput(processed, duration)
                        ", Throughput: $avgThroughput, ETA: $eta"
                    } else {
                        ""
                    }

                    tell(progressPart + throughputPart)
                }
            }
        }
    }
}

interface ProgressBar {
    fun inc(amount: Long = 1)
    fun done()
    fun setTotal(total: Long)
    fun setState(amount: Long, total: Long)
}
