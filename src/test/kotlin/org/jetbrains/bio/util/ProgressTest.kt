package org.jetbrains.bio.util

import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.LoggerContext
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.Appender
import ch.qos.logback.core.OutputStreamAppender
import org.jetbrains.bio.Retry
import org.jetbrains.bio.RetryRule
import org.jetbrains.bio.util.ProgressPart.Companion.fromBoundedProgressString
import org.jetbrains.bio.util.ProgressPart.Companion.fromUnboundedProgressString
import org.junit.After
import org.junit.Before
import org.junit.Rule
import org.junit.Test
import org.slf4j.LoggerFactory
import java.io.ByteArrayOutputStream
import java.io.OutputStream
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.LongAdder
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ProgressFormattingTest {
    @Test
    fun testFormatTime() {
        assertEquals("1 \u03bcs", asTime(1234))
        assertEquals("10 s", asTime(TimeUnit.SECONDS.toNanos(10)))
        assertEquals("1 min", asTime(TimeUnit.SECONDS.toNanos(60)))
        assertEquals("1 min 14 s", asTime(TimeUnit.SECONDS.toNanos(74)))
        assertEquals("1 h", asTime(TimeUnit.SECONDS.toNanos(3640)))
        assertEquals("1 h 1 min", asTime(TimeUnit.SECONDS.toNanos(3680)))
    }

    @Test
    fun testFormatThroughput() {
        assertEquals("1 items/s", asThroughput(10, TimeUnit.SECONDS.toNanos(10)))
        assertEquals("2 items/s", asThroughput(20, TimeUnit.SECONDS.toNanos(10)))
        assertEquals("2 items/ms", asThroughput(25442, TimeUnit.SECONDS.toNanos(10)))
        assertEquals("2544 items/ns", asThroughput(25442, 10))
    }

    @Test
    fun fromBoundedProgressString() {
        assertEquals(
            ProgressPart(1.0, 1, 100, 18),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18 s")
        )

        assertEquals(
            ProgressPart(1.0, 1, 100, 0),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18 ns")
        )
        assertEquals(
            ProgressPart(1.0, 1, 100, 18),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18000000000 ns")
        )

        assertEquals(
            ProgressPart(1.0, 1, 100, 0),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18 ms")
        )
        assertEquals(
            ProgressPart(1.0, 1, 100, 18),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18000 ms")
        )

        assertEquals(
            ProgressPart(1.0, 1, 100, 0),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18 μs")
        )
        assertEquals(
            ProgressPart(1.0, 1, 100, 18),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18000000 μs")
        )

        assertEquals(
            ProgressPart(1.0, 1, 100, 18),
            fromBoundedProgressString("INFO - Progress: 1% (1/100), Elapsed time: 18000000 ?s")
        )
    }

    @Test
    fun fromUnBoundedProgressString() {
        assertEquals(
            ProgressPart(-1.0, 100, -1, 18),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18 s")
        )

        assertEquals(
            ProgressPart(-1.0, 100, -1, 0),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18 ns")
        )
        assertEquals(
            ProgressPart(-1.0, 100, -1, 18),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18000000000 ns")
        )

        assertEquals(
            ProgressPart(-1.0, 100, -1, 0),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18 ms")
        )
        assertEquals(
            ProgressPart(-1.0, 100, -1, 18),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18000 ms")
        )

        assertEquals(
            ProgressPart(-1.0, 100, -1, 0),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18 μs")
        )
        assertEquals(
            ProgressPart(-1.0, 100, -1, 18),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18000000 μs")
        )

        assertEquals(
            ProgressPart(-1.0, 100, -1, 18),
            fromUnboundedProgressString("INFO - Processed items: 100, Elapsed time: 18000000 ?s")
        )
    }

}

data class ProgressPart(
    val percentCompleted: Double,
    val itemsDone: Int, val itemsOverall: Int,
    val elapsedSeconds: Int
) {
    companion object {
        private val PAT_BOUNDED =
            ("Progress: (\\d+(?:\\.\\d+)?)% \\((\\d+)/(\\d+)\\), " +
                    "Elapsed time: (\\d+)(?:\\.\\d+)? ([\\wμ?]+)").toRegex()
        private val PAT_UNBOUNDED =
            ("Processed items: ((?:\\d+,)*\\d+), " +
                    "Elapsed time: (\\d+)(?:\\.\\d+)? ([\\wμ?]+)").toRegex()

        /** Converts (duration, unit) pair to seconds. */
        private fun Pair<Long, String>.toSeconds(): Int {
            val unit = when (second) {
                "\u03bcs", "?s" -> TimeUnit.MICROSECONDS
                "ms" -> TimeUnit.MILLISECONDS
                "s" -> TimeUnit.SECONDS
                "ns" -> TimeUnit.NANOSECONDS
                else -> throw AssertionError("'$second'")
            }

            return unit.toSeconds(first).toInt()
        }

        fun fromBoundedProgressString(str: String): ProgressPart {
            assertTrue(PAT_BOUNDED in str, "This should match regexp: '$str'")
            val match = PAT_BOUNDED.find(str)!!.groups
            return ProgressPart(
                match[1]!!.value.toDouble(),
                match[2]!!.value.toInt(),
                match[3]!!.value.toInt(),
                (match[4]!!.value.toLong() to match[5]!!.value).toSeconds()
            )
        }

        fun fromUnboundedProgressString(str: String): ProgressPart {
            assertTrue(PAT_UNBOUNDED in str, "This should match regexp: '$str'")
            val match = PAT_UNBOUNDED.find(str)!!.groups
            return ProgressPart(
                -1.0, match[1]!!.value.filter { it.isDigit() }.toInt(), -1,
                (match[2]!!.value.toLong() to match[3]!!.value).toSeconds()
            )
        }
    }
}

abstract class ProgressTest {
    private val APPENDER_NAME = "AppenderProgressTest"

    private lateinit var logStream: ByteArrayOutputStream

    /** A list of log output lines. */
    protected val logStrings: List<String>
        get() = logStream.toString().trim().split('\n')
            .filter { it.isNotBlank() }

    /** A list of reports from progress. */
    protected abstract val parts: List<ProgressPart>

    @Before
    fun setUp() {
        logStream = ByteArrayOutputStream()
        addAppender(logStream, APPENDER_NAME)
    }

    @After
    fun tearDown() {
        Logs.getRootLogger().detachAppender(APPENDER_NAME)
    }
}

class SequentialProgressTest : ProgressTest() {

    override val parts: List<ProgressPart>
        get() = logStrings.map {
            fromBoundedProgressString(it)
        }

    @Retry
    @Test
    fun testBasic() {
        // Test is suppose to check that we output not on each update, but only on some
        // updates depending on period setting
        // Initial test failed occasionally because 100 times * 100 ms != 10s, could be 12, 14 or 17s
        // so let's report every 1 second instead of 100ms, but 10 points instead of 1 point.
        val progress = Progress { period = 1 to TimeUnit.SECONDS }.bounded(100)
        for (i in 0 until 10) {
            progress.report(1)
            Thread.sleep(1000)
            progress.report(9)
            Thread.sleep(20)
        }

        progress.done()
        progress.done() // to ensure that this doesn't make any difference

        val uniqueSeconds = parts.asSequence().map { it.elapsedSeconds }.distinct().count()
        println("Unique seconds: $uniqueSeconds\n" +
                "Actual parts: ${parts.joinToString { "${it.elapsedSeconds} {$it}" }}"
        )

        val expPartsRange = 10..12
        val expPartsRangeStr = "[${expPartsRange.start}..${expPartsRange.endInclusive}]"
        assertTrue(
            logStrings.size in expPartsRange,
            "Timer progress lines count is expected to be $expPartsRangeStr, but was ${logStrings.size}"
        )
        assertEquals(1.0, parts[0].percentCompleted)
        assertEquals(0, parts[0].elapsedSeconds)
        assertEquals(1, parts[0].itemsDone)
        assertEquals(100.0, parts.last().percentCompleted)

        val elapsedSecs = parts.last().elapsedSeconds
        val expElapsedRange = 10..11
        val expElapsedRangeStr = "[${expElapsedRange.start}..${expElapsedRange.endInclusive}]"
        assertTrue(
            elapsedSecs in expElapsedRange,
            "Elapsed times expected to be $expElapsedRangeStr sec, " +
                    "but was $elapsedSecs"
        )
        assertEquals(100, parts.last().itemsDone)
        assertTrue("[done]" in logStrings.last())
    }

    @Test
    fun testEmpty() {
        val progress = Progress { period = 1 to TimeUnit.SECONDS }.bounded(0)
        progress.done()
        assertEquals(1, logStrings.size)
        assertTrue("[done]" in logStrings.first())
        assertTrue("Progress: 100%" in logStrings.first())
    }


    @Test
    fun testBasicIncremental() {
        // Test is suppose to check that we output not on each update, but only on some
        // updates depending on period setting
        // Initial test failed occasionally because 100 times * 100 ms != 10s, could be 12, 14 or 17s
        // so let's report every 1 second instead of 100ms, but 10 points instead of 1 point.
        // Also lets report this 10 points using 5+5 invocations.
        val progress = Progress { period = 1 to TimeUnit.SECONDS }.bounded(100)
        for (i in 0 until 10) {
            // incremental
            (0 until 5).forEach { progress.report() }
            Thread.sleep(1000)
            (0 until 5).forEach { progress.report() }
            Thread.sleep(20)
        }
        progress.done()
        progress.done() // to ensure that this doesn't make any difference

        val uniqueSeconds = parts.asSequence().map { it.elapsedSeconds }.distinct().count()
        println("Unique seconds: $uniqueSeconds\n" +
                "Actual parts: ${parts.joinToString { "${it.elapsedSeconds} {$it}" }}"
        )

        val expPartsRange = 10..12
        val expPartsRangeStr = "[${expPartsRange.start}..${expPartsRange.endInclusive}]"
        assertTrue(
            logStrings.size in expPartsRange,
            "Timer progress lines count is expected to be $expPartsRangeStr," +
                    " but was ${logStrings.size}"
        )
        assertEquals(1.0, parts[0].percentCompleted)
        // can be <1s, e.g. 258ms.
        // assertEquals(0, parts.get(0).elapsed);
        assertEquals(1, parts[0].itemsDone)
        assertEquals(100.0, parts.last().percentCompleted)
        assertEquals(100, parts.last().itemsDone)

        val elapsedSecs = parts.last().elapsedSeconds
        val expElapsedRange = 10..11
        val expElapsedRangeStr = "[${expElapsedRange.start}..${expElapsedRange.endInclusive}]"
        assertTrue(
            elapsedSecs in expElapsedRange,
            "Elapsed times expected to be $expElapsedRangeStr sec, but was $elapsedSecs"
        )
        assertTrue("[done]" in logStrings.last())
    }
}

class ParallelProgressTest : ProgressTest() {
    @get:Rule
    var rule = RetryRule(3)

    override val parts: List<ProgressPart>
        get() = logStrings.map {
            fromUnboundedProgressString(it)
        }

    @Test
    fun testParallel() {
        val progress = Progress { period = 1 to TimeUnit.SECONDS }.unbounded()
        val pool = Executors.newFixedThreadPool(THREADS)
        val adder = LongAdder()

        val endTime = System.nanoTime() + TEST_DURATION
        val callables = (0..THREADS).map {
            Callable<Void> {
                var lastTime = System.nanoTime()
                while (lastTime < endTime) {
                    // busy wait here, woo-hoo
                    while (System.nanoTime() - lastTime < WAIT_NANOS) {
                    }
                    adder.add(1)
                    progress.report()
                    lastTime = System.nanoTime()
                }

                null
            }
        }

        pool.invokeAll(callables)
        progress.done()
        progress.done() // to ensure that this doesn't make any difference

        val uniqueSeconds = parts.asSequence().map { it.elapsedSeconds }.distinct().count()
        println("Unique seconds: $uniqueSeconds\n" +
                "Actual parts: ${parts.joinToString { "${it.elapsedSeconds} {$it}" }}"
        )

        val allSeconds = parts.size
        assertEquals(
            allSeconds, uniqueSeconds + 1,
            "shouldn't report progress more often than once in second, except [done]"
        )
        assertTrue("[done]" in logStrings[logStrings.size - 1])
        assertEquals(
            adder.sum(), parts[parts.size - 1].itemsDone.toLong(),
            "should report all progress"
        )
    }

    @Retry
    @Test
    fun testParallelIncremental() {
        val progress = Progress { period = 1 to TimeUnit.SECONDS }.unbounded()
        val pool = Executors.newFixedThreadPool(THREADS)
        val adder = LongAdder()

        val endTime = System.nanoTime() + TEST_DURATION
        val callables = (0..THREADS).map {
            Callable<Void> {
                var lastTime = System.nanoTime()
                while (lastTime < endTime) {
                    // busy wait here, woo-hoo
                    while (System.nanoTime() - lastTime < WAIT_NANOS) {
                    }
                    adder.add(1)
                    progress.report()
                    lastTime = System.nanoTime()
                }

                null
            }
        }

        pool.invokeAll(callables)
        progress.done()
        progress.done() // to ensure that this doesn't make any difference

        val uniqueSeconds = parts.asSequence().map { it.elapsedSeconds }.distinct().count()
        println("Unique seconds: $uniqueSeconds\n" +
                "Actual parts: ${parts.joinToString { "${it.elapsedSeconds} {$it}" }}"
        )

        assertEquals(parts.size, uniqueSeconds + 1,
            "shouldn't report progress more often than once in second, except [done].\n" +
                    "Actual parts: ${parts.joinToString { "${it.elapsedSeconds} {$it}" }}"
        )
        assertTrue("[done]" in logStrings[logStrings.size - 1])
        assertEquals(
            adder.sum(), parts.last().itemsDone.toLong(),
            "should report all progress"
        )
    }

    companion object {
        private val THREADS = Runtime.getRuntime().availableProcessors()

        // we use 9900 ms here instead of just one second because otherwise there
        // is a race on 10th second (10th second will or will not be reported)
        private val TEST_DURATION = TimeUnit.MILLISECONDS.toNanos(9900)
        private const val WAIT_NANOS = 10
    }
}

class MultiTaskProgressTest {
    private lateinit var logStream: ByteArrayOutputStream
    private val APPENDER_NAME = "AppenderMultiTaskProgressTest"

    @Before
    fun setUp() {
        logStream = ByteArrayOutputStream()
        addAppender(logStream, APPENDER_NAME)
    }

    @After
    fun tearDown() {
        Logs.getRootLogger().detachAppender(APPENDER_NAME)
    }

    @Test
    fun testFinish() {
        val tasks = 10
        val numbers = LongArray(tasks) { (10 * it + 1).toLong() }
        (0 until tasks).forEach {
            MultitaskProgress.addTask("task_$it", numbers[it], quite = false)
        }

        (0 until tasks - 1).forEach {
            MultitaskProgress.reportTask("task_$it", numbers[it] - 1)
            MultitaskProgress.finishTask("task_$it")
            Thread.sleep(100)
        }
        MultitaskProgress.finishTask("task_${tasks - 1}")
        val log = logStream.toString()
        assertTrue("Running 9 tasks: 0.00% (0/459)" in log)
        assertTrue("Running 2 tasks: 61.95% (280/452)" in log)
        assertTrue("100.00% (360/360)" in log)
    }
}

internal fun addAppender(appenderStream: OutputStream, appenderName: String): Appender<ILoggingEvent> {
    val loggerContext = LoggerFactory.getILoggerFactory() as LoggerContext

    val logEncoder = PatternLayoutEncoder().apply {
        context = loggerContext
        pattern = "%-5level - %msg%n"
    }
    logEncoder.start()

    val appender = OutputStreamAppender<ILoggingEvent>().apply {
        context = loggerContext
        encoder = logEncoder
        outputStream = appenderStream
        name = appenderName
    }
    appender.start()

    loggerContext.getLogger(Logger.ROOT_LOGGER_NAME).addAppender(appender)
    return appender
}
