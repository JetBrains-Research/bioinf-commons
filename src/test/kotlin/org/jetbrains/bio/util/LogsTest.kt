package org.jetbrains.bio.util

import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.jetbrains.bio.Tests.assertIn
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.PrintStream
import java.lang.reflect.InvocationTargetException
import java.util.concurrent.ExecutionException
import kotlin.test.assertEquals

class LogsTest {
    @Test
    fun getMessageRuntime() {
        val original = RuntimeException("MESSAGE")
        assertEquals(
            """ERROR MESSAGE
java.lang.RuntimeException: MESSAGE
${'\t'}at org.jetbrains.bio.util.LogsTest""",
            Logs.getMessage(RuntimeException(original), includeStackTrace = true)
                    .replaceAfter("at org.jetbrains.bio.util.LogsTest", "")
                    .replace("\r", "")
        )
    }

    @Test
    fun getMessageExecution() {
        val original = RuntimeException()
        assertEquals(
            """ERROR
java.lang.RuntimeException
${'\t'}at org.jetbrains.bio.util.LogsTest""",
            Logs.getMessage(ExecutionException(original), includeStackTrace = true)
                    .replaceAfter("at org.jetbrains.bio.util.LogsTest", "")
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
            """ERROR
java.lang.NullPointerException
${'\t'}at org.jetbrains.bio.util.LogsTest""",
            Logs.getMessage(NullPointerException(), includeStackTrace = true)
                    .replaceAfter("at org.jetbrains.bio.util.LogsTest", "")
                    .replace("\r", "")
        )
    }

    @Test
    fun testLayout() {
        val (out, _) = captureLoggingOutput {
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
        val (out, _) = captureLoggingOutput {
            Logs.addConsoleAppender(Level.DEBUG)
            LOG.debug("foobar")
        }
        assertIn(
            Regex("\\[.* \\d\\d:\\d\\d:\\d\\d] DEBUG LogsTest foobar\n"),
            out.replace("\r", "")
        )
    }


    companion object {
        private val LOG = Logger.getLogger(LogsTest::class.java)
        private val OUT = System.out
        private val ERR = System.err

        /**
         * Useful for checking logs.
         *
         * Runs [block] and returns the contents of stdout and stderr by redirecting them
         * and resetting the console appender.
         *
         * Clears redirection and recreates appender
         * regardless of whether [block] completed normally or threw.
         */
        fun captureLoggingOutput(block: () -> Unit): Pair<String, String> {
            val outStream = ByteArrayOutputStream()
            val errStream = ByteArrayOutputStream()
            System.setOut(PrintStream(outStream))
            System.setErr(PrintStream(errStream))
            Logs.replaceConsoleAppender(Level.INFO)
            try {
                block()
            } finally {
                System.setOut(OUT)
                System.setErr(ERR)
                Logs.replaceConsoleAppender(Level.INFO)
            }
            return String(outStream.toByteArray()) to String(errStream.toByteArray())
        }
    }
}