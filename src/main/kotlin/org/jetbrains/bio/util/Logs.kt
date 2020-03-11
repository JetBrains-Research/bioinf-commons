package org.jetbrains.bio.util

import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.LoggerContext
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.spi.IThrowableProxy
import ch.qos.logback.classic.spi.ThrowableProxy
import ch.qos.logback.core.ConsoleAppender
import ch.qos.logback.core.CoreConstants
import ch.qos.logback.core.LayoutBase
import ch.qos.logback.core.encoder.LayoutWrappingEncoder
import ch.qos.logback.core.rolling.RollingFileAppender
import ch.qos.logback.core.rolling.TimeBasedRollingPolicy
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.io.ByteArrayOutputStream
import java.io.OutputStream
import java.io.PrintStream
import java.lang.reflect.InvocationTargetException
import java.nio.file.Path
import java.text.DateFormat
import java.text.SimpleDateFormat
import java.util.*
import java.util.concurrent.ExecutionException


/**
 * This is a logging facade to logback logger
 */
object Logs {
    /* Same pattern as in logback.xml configs */
    const val PATTERN = "%date %level [%thread] %logger{10} [%file:%line] %msg%n"

    /**
     * Adds a new console appender with logging at [level]. If a console appender was already added,
     * removes it first. Useful when you've redirected [System.out], since the new appender
     * will write to the new stdout.
     */
    fun addConsoleAppender(level: Level) {
        val loggerContext = LoggerFactory.getILoggerFactory() as LoggerContext
        val rootLogger = loggerContext.getLogger(Logger.ROOT_LOGGER_NAME)
        rootLogger.level = ch.qos.logback.classic.Level.toLevel(level.name)

        // Remove other Console appenders
        for (appender in rootLogger.iteratorForAppenders().asSequence().toList()) {
            if (appender is ConsoleAppender) {
                rootLogger.detachAppender(appender)
            }
        }
        val logEncoder = LayoutWrappingEncoder<ILoggingEvent>().apply {
            this.layout = PRETTY_PRINT_LAYOUT
        }
        logEncoder.start()

        val logConsoleAppender: ConsoleAppender<ILoggingEvent> = ConsoleAppender<ILoggingEvent>().apply {
            context = loggerContext
            name = "console"
            setEncoder(logEncoder)
        }
        logConsoleAppender.start()

        rootLogger.addAppender(logConsoleAppender)
    }

    fun getRootLogger(): Logger {
        val loggerContext = LoggerFactory.getILoggerFactory() as LoggerContext
        return loggerContext.getLogger(Logger.ROOT_LOGGER_NAME)
    }

    fun quiet() {
        val nullPrintStream = PrintStream(object : OutputStream() {
            override fun write(b: Int) {}
        })
        System.setOut(nullPrintStream)
        System.setErr(nullPrintStream)
        addConsoleAppender(Level.INFO)
    }

    /**
     * Unwraps [RuntimeException], [ExecutionException], [InvocationTargetException]
     * to get original message
     */
    fun getMessage(throwable: Throwable, includeStackTrace: Boolean = false): String {
        return getMessage(ThrowableProxy(throwable), includeStackTrace)
    }

    fun getMessage(throwable: IThrowableProxy, includeStackTrace: Boolean = false): String {
        var t: IThrowableProxy = throwable
        while (true) {
            if (isAux(t)) {
                val cause = t.cause
                if (cause != null) {
                    t = cause
                    continue
                }
            }
            val stackTrace = t.stackTraceElementProxyArray.iterator().asSequence().map {
                it.stackTraceElement.toString()
            }.joinToString("\n")
            if (isAux(t)) {
                return "ERROR${if (t.message != null) " ${t.message}" else ""}" +
                        (if (includeStackTrace) "\n$stackTrace" else "")
            }
            return "ERROR ${t.className.substringAfterLast('.')}" +
                    (if (t.message != null) " ${t.message}" else "") +
                    (if (includeStackTrace) "\n$stackTrace" else "")

        }
    }

    private fun isAux(t: IThrowableProxy) =
            t.className in setOf(
                    RuntimeException::class.qualifiedName,
                    ExecutionException::class.qualifiedName,
                    InvocationTargetException::class.qualifiedName)

    fun addLoggingToFile(path: Path) {
        val loggerContext = LoggerFactory.getILoggerFactory() as LoggerContext

        val rootLogger = loggerContext.getLogger(Logger.ROOT_LOGGER_NAME)
        val logEncoder = PatternLayoutEncoder().apply {
            context = loggerContext
            pattern = PATTERN
        }
        logEncoder.start()

        val logFileAppender = RollingFileAppender<ILoggingEvent>().apply {
            context = loggerContext
            name = "logFile ${path.fileName}"
            encoder = logEncoder
            file = path.toString()
            // Recreate log file on start
            setAppend(false)
        }

        val logFilePolicy = TimeBasedRollingPolicy<ILoggingEvent>()
        logFilePolicy.context = loggerContext
        logFilePolicy.setParent(logFileAppender)
        logFilePolicy.fileNamePattern = "${path.parent}/${path.name}-%d{yyyy-MM-dd_HH}.${path.extension}"
        logFilePolicy.start()

        logFileAppender.rollingPolicy = logFilePolicy
        logFileAppender.start()

        rootLogger.addAppender(logFileAppender)
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
        addConsoleAppender(Level.INFO)
        try {
            block()
        } finally {
            System.setOut(out)
            System.setErr(err)
            addConsoleAppender(Level.INFO)
        }
        return String(outStream.toByteArray()) to String(errStream.toByteArray())
    }

    val PRETTY_PRINT_LAYOUT = object : LayoutBase<ILoggingEvent>() {

        override fun doLayout(event: ILoggingEvent): String {
            return "[${timePart(event)}]${logLevelPart(event)} ${event.message}" +
                    "${throwablePart(event)}${CoreConstants.LINE_SEPARATOR}"
        }

        private fun timePart(event: ILoggingEvent): String {
            val date = Date(event.timeStamp)
            return "${DateFormat.getDateInstance(
                    DateFormat.MEDIUM, Locale.getDefault()
            ).format(date)} ${SimpleDateFormat("HH:mm:ss").format(date)}"
        }


        private fun logLevelPart(event: ILoggingEvent): String {
            return if (event.level != ch.qos.logback.classic.Level.INFO) {
                " ${event.level} ${event.loggerName.substringAfterLast('.')}"
            } else ""
        }

        private fun throwablePart(event: ILoggingEvent): String {
            return if (event.throwableProxy != null) {
                "${CoreConstants.LINE_SEPARATOR}Caused by: ${
                getMessage(event.throwableProxy, includeStackTrace = true)}"
            } else
                ""
        }
    }

}
