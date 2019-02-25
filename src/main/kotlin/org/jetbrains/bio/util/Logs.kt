package org.jetbrains.bio.util

import org.apache.log4j.*
import org.apache.log4j.spi.LoggingEvent
import java.io.IOException
import java.io.PrintWriter
import java.io.StringWriter
import java.lang.reflect.InvocationTargetException
import java.nio.file.Path
import java.text.DateFormat
import java.text.SimpleDateFormat
import java.util.*
import java.util.concurrent.ExecutionException

object Logs {
    const val CONSOLE_APPENDER = "CONSOLE_APPENDER"

    val LAYOUT: Layout = object : Layout() {

        override fun activateOptions() {
            // Ignore
        }

        private fun timePart(event: LoggingEvent): String {
            val date = Date(event.timeStamp)
            return "${DateFormat.getDateInstance(
                DateFormat.MEDIUM, Locale.getDefault()
            ).format(date)} ${SimpleDateFormat("HH:mm:ss").format(date)}"
        }


        private fun logLevelPart(event: LoggingEvent): String {
            return if (event.getLevel() != Level.INFO) {
                " ${event.getLevel()} ${event.loggerName.substringAfterLast('.')}"
            } else ""
        }

        private fun throwablePart(event: LoggingEvent): String {
            return if (event.throwableInformation != null) {
                "${LINE_SEP}Caused by: ${getMessage(
                    event.throwableInformation.throwable, includeStackTrace = true
                )}"
            } else
                ""
        }

        override fun format(event: LoggingEvent): String {
            return "[${timePart(event)}]${logLevelPart(event)} ${event.message}" +
                    "${throwablePart(event)}$LINE_SEP"
        }

        override fun ignoresThrowable(): Boolean {
            return false
        }
    }


    /**
     * Adds a console appender with logging at [level]. If a console appender was already added,
     * just changes the logging level.
     */
    fun addConsoleAppender(level: Level) {
        with(Logger.getRootLogger()) {
            this.level = level
            if (getAppender(CONSOLE_APPENDER) == null) {
                addAppender(ConsoleAppender(LAYOUT).apply { name = CONSOLE_APPENDER })
            }
        }
    }

    /**
     * Adds a new console appender with logging at [level]. If a console appender was already added,
     * removes it first. Useful when you've redirected [System.out], since the new appender
     * will write to the new stdout.
     */
    fun replaceConsoleAppender(level: Level) {
        with(Logger.getRootLogger()) {
            this.level = level
            if (getAppender(CONSOLE_APPENDER) != null) {
                removeAppender(CONSOLE_APPENDER)
            }
            addAppender(ConsoleAppender(LAYOUT).apply { name = CONSOLE_APPENDER })
        }
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

}
