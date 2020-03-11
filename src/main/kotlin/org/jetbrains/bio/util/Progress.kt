package org.jetbrains.bio.util

import com.google.common.annotations.VisibleForTesting
import org.slf4j.LoggerFactory
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeUnit.*
import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.atomic.LongAccumulator
import java.util.function.LongBinaryOperator
import kotlin.math.max

/**
 * Progress tracker.
 *
 * The progress is measured in the number of processed items.
 *
 * @author Dmitry Groshev
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 */
abstract class Progress protected constructor(
        /** Title to be displayed. */
        protected val title: String,
        /** The amount of time between progress reports. */
        private val periodNanos: Long,
        /** Report into log file or stderr. Progress is often reported to stderr and stdout is kept clean */
        private val reportToStderr: Boolean
) {
    private val startTime = System.nanoTime()
    private var lastDuration = AtomicLong(-periodNanos)

    protected var accumulator: LongAccumulator = LongAccumulator(LongBinaryOperator { a, b -> a + b }, 0)

    @Volatile
    protected var isDone = false


    protected fun getDuration(): Long {
        return getDuration(startTime, lastDuration, periodNanos)
    }

    protected abstract fun format(duration: Long, processed: Long): String

    protected fun tell(message: String) = when {
        reportToStderr -> System.err.println("$title: $message")
        else -> LOG.info("$title: $message")
    }

    abstract fun report(update: Long = 1)

    /**
     * Marks progress as done. No further updates would be accumulated.
     */
    fun done(message: String = "") {
        if (!isDone) {
            isDone = true
            val duration = System.nanoTime() - startTime
            tell(format(duration, processedItems()) + " [done]")
            if (message.isNotEmpty()) {
                tell(message)
            }
        }
    }


    class Builder internal constructor() {
        var title: String? = null

        var period: Pair<Int, TimeUnit> = 10 to SECONDS

        /**
         * @param totalItems Total items number
         * @param reportToStderr If true report to [System.err] instead of LOG file
         * @param percentsOnly If true report progress only in percents and omit items number. E.g. useful
         *  when performance depends on locus length, so items number is summary loci length and progress reports
         *  locus length on update. Items number in such cases is too big and useless for user.
         */
        fun bounded(
                totalItems: Long,
                reportToStderr: Boolean=false,
                percentsOnly: Boolean=false
        ): Progress = Bounded(
                title,
                period.second.toNanos(period.first.toLong()),
                totalItems,
                reportToStderr = reportToStderr,
                percentsOnly = percentsOnly
        )

        fun unbounded(reportToStderr: Boolean=false): Progress = Unbounded(
                title,
                period.second.toNanos(period.first.toLong()),
                reportToStderr = reportToStderr
        )
    }

    abstract fun processedItems(): Long

    companion object {
        @VisibleForTesting
        internal val LOG = LoggerFactory.getLogger(Progress::class.java)

        operator fun invoke(block: Builder.() -> Unit) = Builder().apply(block)

        /**
         * This is the main "guarding" method for Progress. Its purpose is to
         * avoid excessive progress reporting (using [periodNanos]). It is
         * called for every update and therefore should be fast.
         *
         * The gist of it is an atomic CAS. It goes like this.
         *
         *   1. Calculate elapsed time from the creation of this [Progress].
         *   2. Get current value of [lastDuration] field that encodes the
         *      last reporting time. It's a *volatile* read, therefore it's
         *      visible across threads.
         *   3. Compare calculated duration with a threshold and return -1
         *      if it's too early to report progress.
         *   4. This is the tricky part. It's probable that more than one
         *      thread read lastDuration and now they are racing to set
         *      `snapshot`. We are using CAS here to gently serialize this
         *      scenario. If some other thread already changed [lastDuration],
         *      it means that other thread successfully reported and this
         *      thread shouldn't do this.
         *
         * Here is a sketchy proof that there this will work:
         * - atomics use volatiles
         * - therefore, reads and writes to [lastDuration] are strictly
         *   ordered by happens-before relation, and [System.nanoTime] calls
         *   are strictly ordered as well. It's guaranteed by JM that the
         *   ordering will be like this:
         *     - thread1 calls nanoTime;
         *     - thread1 writes [lastDuration];
         *     - thread2 reads [lastDuration];
         *     - thread2 calls nanoTime.
         * - we can show that duration is always strictly larger than `snapshot` on CAS:
         * - it's either larger because [System.nanoTime] returns larger
         *   value than previous one;
         * - or, if it returns the same value as in a different thread, duration
         *   will be equal to `snapshot` and "duration < timeoutNanos + snapshot"
         *   will be true and this method will return -1;
         * - therefore, *every* CAS will increase value of [lastDuration];
         * - therefore, there is no possibility that two threads will
         *   successfully update `snapshot` at the same time.
         *
         * Regarding return value: [System.nanoTime] is not free, therefore
         * it's best to avoid extra calls, therefore this method performs checks
         * AND returns duration.
         */
        fun getDuration(startTime: Long,
                        lastDuration: AtomicLong,
                        periodNanos: Long,
                        allowSorterPeriod: Boolean = false): Long {
            val lastDurationValue = lastDuration.get()
            val duration = System.nanoTime() - startTime
            if (duration < periodNanos + lastDurationValue && !allowSorterPeriod) {
                return -1
            }
            if (!lastDuration.compareAndSet(lastDurationValue, duration)) {
                return -1
            }
            return duration
        }

    }
}

private class Bounded(
        title: String?,
        timeOutNanos: Long,
        private val totalItems: Long,
        reportToStderr: Boolean,
        private val percentsOnly: Boolean
) : Progress(title ?: "Progress", timeOutNanos, reportToStderr = reportToStderr) {

    override fun processedItems() = totalItems

    @Volatile
    private var progressPercent: Long = -1

    override fun report(update: Long) {
        if (update == 0L) {
            // nothing to do
            return
        }

        check(!isDone) {
            "$title progress is done: ${accumulator.get()} / $totalItems. Updated requested: $update"
        }
        accumulator.accumulate(update)

        val duration = getDuration()
        if (duration == -1L) {
            return
        }

        val processed = accumulator.get()
        val processedPercent = if (totalItems > 0) (100.0 * processed / totalItems).toLong() else 100L
        if (processedPercent > progressPercent) {
            progressPercent = processedPercent
            if (processed < totalItems) {
                tell(format(duration, processed))
            } else {
                check(processed == totalItems) {
                    "$title progress overflow, expected $totalItems, got $processed"
                }
                done()
            }
        }
    }

    override fun format(duration: Long, processed: Long): String {
        val progressPart = if (totalItems > 0) {
            val itemsStr = when {
                percentsOnly -> ""
                else -> " (%,d/%,d)".format(processed, totalItems)
            }
            "%.0f%%$itemsStr, Elapsed time: ${asTime(duration)}".format(processed * 100.0 / totalItems)
        } else {
            "100%, Elapsed time: ${asTime(duration)}"
        }
        val throughputPart = if (processed > 1) {
            val eta = asTime((duration / (processed.toDouble() / totalItems)).toLong() - duration)
            val throughput = asThroughput(processed, duration)
            ", Throughput: $throughput, ETA: $eta"
        } else {
            ""
        }
        return progressPart + throughputPart
    }
}

private class Unbounded(
        title: String?,
        timeOutNanos: Long,
        reportToStderr: Boolean
) : Progress(title ?: "Processed items", timeOutNanos, reportToStderr = reportToStderr) {

    override fun processedItems() = accumulator.get()

    override fun report(update: Long) {
        check(!isDone) { "$title progress is done" }
        accumulator.accumulate(update)

        val duration = getDuration()
        if (duration == -1L) {
            return
        }

        val processed = accumulator.get()
        tell(format(duration, processed))
    }

    override fun format(duration: Long, processed: Long): String {
        return if (processed > 1) {
            val avgThroughput = asThroughput(processed, duration)
            "%,d, Elapsed time: ${asTime(duration)}, Throughput: $avgThroughput".format(processed)
        } else {
            "%,d, Elapsed time: ${asTime(duration)}".format(processed)
        }
    }
}

private fun Long.chooseUnit(): TimeUnit {
    return when {
        DAYS.convert(this, NANOSECONDS) > 0 -> DAYS
        HOURS.convert(this, NANOSECONDS) > 0 -> HOURS
        MINUTES.convert(this, NANOSECONDS) > 0 -> MINUTES
        SECONDS.convert(this, NANOSECONDS) > 0 -> SECONDS
        MILLISECONDS.convert(this, NANOSECONDS) > 0 -> MILLISECONDS
        MICROSECONDS.convert(this, NANOSECONDS) > 0 -> MICROSECONDS
        else -> NANOSECONDS
    }
}

private fun TimeUnit.abbreviate(): String {
    return when (this) {
        NANOSECONDS -> "ns"
        MICROSECONDS -> "\u03bcs" // μs
        MILLISECONDS -> "ms"
        SECONDS -> "s"
        MINUTES -> "min"
        HOURS -> "h"
        DAYS -> "d"
        else -> throw AssertionError()
    }
}

fun asThroughput(items: Long, nanos: Long): String {
    val digits = 10  // <digits> items/unit.
    val unit = (nanos * digits / items).chooseUnit()
    val amount = items * NANOSECONDS.convert(1, unit) / nanos
    return if (amount == 0L) {
        "${unit.convert(nanos, NANOSECONDS) / items} ${unit.abbreviate()}/item"
    } else {
        "$amount items/${unit.abbreviate()}"
    }
}

fun asTime(nanos: Long): String {
    val unit = nanos.chooseUnit()
    val subUnit = when (unit) {
        DAYS -> HOURS
        HOURS -> MINUTES
        else -> SECONDS // don't use subunit at all if the duration is less than 1 minute
    }
    val duration = nanos / NANOSECONDS.convert(1, unit)
    val subDuration = max(
            0,
            (nanos - duration * NANOSECONDS.convert(1, unit)) / NANOSECONDS.convert(1, subUnit)
    )
    return when {
        subDuration != 0L -> "%d %s %d %s".format(duration, unit.abbreviate(), subDuration, subUnit.abbreviate())
        else -> "%d %s".format(duration, unit.abbreviate())
    }
}