package org.jetbrains.bio.util

import com.google.common.base.Stopwatch
import org.apache.log4j.Level
import org.apache.log4j.Logger
import java.util.concurrent.CancellationException

/**
 *  Measures the running time of a given possibly impure [block].
 */
inline fun <R> Logger.time(level: Level = Level.DEBUG,
                           message: String = "",
                           block: () -> R): R {
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
        log(level, "$message: [FAILED] after $stopwatch", e)
        throw e
    }
    stopwatch.stop()
    log(level, "$message: done in $stopwatch")
    return res
}