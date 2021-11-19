package org.jetbrains.bio.util

import com.google.common.base.Stopwatch
import org.slf4j.Logger
import org.slf4j.event.Level
import java.util.concurrent.CancellationException

/**
 *  Measures the running time of a given possibly impure [block].
 */
fun <R> Logger.time(
    level: Level = Level.DEBUG,
    message: String = "",
    block: () -> R
): R {
    log(level, "$message...")
    val stopwatch = Stopwatch.createStarted()

    val res = try {
        block()
    } catch (e: CancellationException) {
        stopwatch.stop()
        log(level, "$message: [CANCELED] after $stopwatch")
        throw e
    } catch (e: Exception) {
        stopwatch.stop()
        error("$message: [FAILED] after $stopwatch", e)
        throw e
    }
    stopwatch.stop()
    log(level, "$message: done in $stopwatch")
    return res
}

fun Logger.log(level: Level, msg: String) {
    when (level) {
        Level.INFO -> info(msg)
        Level.DEBUG -> debug(msg)
        Level.ERROR -> error(msg)
        else -> throw IllegalArgumentException("Unknown log level $level")
    }
}