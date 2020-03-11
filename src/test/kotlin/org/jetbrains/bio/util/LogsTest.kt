package org.jetbrains.bio.util

import org.jetbrains.bio.Tests.assertIn
import org.junit.Test
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.lang.reflect.InvocationTargetException
import java.util.concurrent.ExecutionException
import kotlin.test.assertEquals

/**
 * Test for Layout format used in tools with Backlog facade
 */
class LogsTest {
    @Test
    fun getMessageRuntime() {
        val original = RuntimeException("MESSAGE")
        assertEquals(
                """ERROR MESSAGE
org.jetbrains.bio.util.LogsTest""",
                Logs.getMessage(RuntimeException(original), includeStackTrace = true)
                        .replaceAfter("org.jetbrains.bio.util.LogsTest", "")
                        .replace("\r", "")
        )
    }

    @Test
    fun getMessageExecution() {
        val original = RuntimeException()
        assertEquals(
                """ERROR
org.jetbrains.bio.util.LogsTest""",
                Logs.getMessage(ExecutionException(original), includeStackTrace = true)
                        .replaceAfter("org.jetbrains.bio.util.LogsTest", "")
                        .replace("\r", "")
        )
    }

    @Test
    fun getMessageInvocationTarget() {
        val original = RuntimeException()
        assertEquals("ERROR", Logs.getMessage(InvocationTargetException(original)))
    }

    @Test
    fun testNull() {
        assertEquals(
                """ERROR NullPointerException
org.jetbrains.bio.util.LogsTest""",
                Logs.getMessage(NullPointerException(), includeStackTrace = true)
                        .replaceAfter("org.jetbrains.bio.util.LogsTest", "")
                        .replace("\r", "")
        )
    }

    @Test
    fun testLayout() {
        val (out, _) = Logs.captureLoggingOutput {
            Logs.addConsoleAppender(Level.INFO)
            LOG.info("foobar")
        }
        assertIn(
                Regex("\\[.* \\d\\d:\\d\\d:\\d\\d] foobar\n"),
                out.replace("\r", "")
        )
    }

    @Test
    fun testDebugLayout() {
        val (out, _) = Logs.captureLoggingOutput {
            Logs.addConsoleAppender(Level.DEBUG)
            LOG.debug("foobar")
        }
        assertIn(
                Regex("\\[.* \\d\\d:\\d\\d:\\d\\d] DEBUG LogsTest foobar\n"),
                out.replace("\r", "")
        )
    }


    companion object {
        private val LOG = LoggerFactory.getLogger(LogsTest::class.java)
    }
}