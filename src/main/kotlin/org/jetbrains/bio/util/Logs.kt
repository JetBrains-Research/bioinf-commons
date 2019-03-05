package org.jetbrains.bio.util

import org.apache.log4j.*
import java.io.*
import java.lang.reflect.InvocationTargetException
import java.nio.file.Path
import java.util.concurrent.ExecutionException

object Logs {
    val LAYOUT: Layout = LogsLayout()

    /**
     * Adds a new console appender with logging at [level]. If a console appender was already added,
     * removes it first. Useful when you've redirected [System.out], since the new appender
     * will write to the new stdout.
     */
    fun addConsoleAppender(level: Level) {
        with(Logger.getRootLogger()) {
            for (appender in allAppenders) {
                if (appender is ConsoleAppender) {
                    removeAppender(appender)
                }
            }
            this.level = level
            addAppender(ConsoleAppender(LAYOUT))
        }
    }

    fun quiet() {
        val nullPrintStream = PrintStream(object : OutputStream() {
            override fun write(b: Int) {}
        })
        System.setOut(nullPrintStream)
        System.setErr(nullPrintStream)
        Logs.addConsoleAppender(Level.INFO)
    }

    /**
     * Unwraps [RuntimeException], [ExecutionException], [InvocationTargetException]
     * to get original message
     */
    fun getMessage(throwable: Throwable, includeStackTrace: Boolean = false): String {
        var t: Throwable = throwable
        while (true) {
            if (t is RuntimeException || t is ExecutionException || t is InvocationTargetException) {
                val cause = t.cause
                if (cause != null) {
                    t = cause
                    continue
                }
            }
            val stackTraceWriter = StringWriter()
            t.printStackTrace(PrintWriter(stackTraceWriter))
            if (t is RuntimeException || t is ExecutionException || t is InvocationTargetException) {
                return "ERROR${if (t.message != null) " ${t.message}" else ""}" +
                        (if (includeStackTrace) "\n$stackTraceWriter" else "")
            }
            return "ERROR ${t.javaClass.simpleName}${if (t.message != null) " ${t.message}" else ""}" +
                    (if (includeStackTrace) "\n$stackTraceWriter" else "")

        }
    }

    fun addLoggingToFile(path: Path) {
        val fileAppender: FileAppender?
        try {
            // val pattern = PatternLayout("%d [%7r] %6p - %30.30c - %m\n")
            val pattern = LAYOUT
            fileAppender = FileAppender(pattern, path.toString(), false)
            Logger.getRootLogger().addAppender(fileAppender)
        } catch (e: IOException) {
            System.err.println("Failed to create log file: $path")
        }
    }


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
        val out = System.out
        val err = System.err
        val outStream = ByteArrayOutputStream()
        val errStream = ByteArrayOutputStream()
        System.setOut(PrintStream(outStream))
        System.setErr(PrintStream(errStream))
        Logs.addConsoleAppender(Level.INFO)
        try {
            block()
        } finally {
            System.setOut(out)
            System.setErr(err)
            Logs.addConsoleAppender(Level.INFO)
        }
        return String(outStream.toByteArray()) to String(errStream.toByteArray())
    }

}
