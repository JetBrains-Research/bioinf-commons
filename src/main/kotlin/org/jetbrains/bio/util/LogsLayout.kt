package org.jetbrains.bio.util

import org.apache.log4j.Layout
import org.apache.log4j.Level
import org.apache.log4j.spi.LoggingEvent
import java.text.DateFormat
import java.text.SimpleDateFormat
import java.util.*

/**
 * Layout used in log4j.xml configs
 */
@Suppress("unused")
class LogsLayout : Layout() {

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
            "${LINE_SEP}Caused by: ${Logs.getMessage(
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
