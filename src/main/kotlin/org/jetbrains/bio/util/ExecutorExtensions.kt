package org.jetbrains.bio.util

import java.util.concurrent.Callable
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import kotlin.math.min

const val PARALLELISM = "parallelism.level"

/**
 * Required parallelism level configured by application
 */
fun parallelismLevel(): Int {
    val p = System.getProperty(PARALLELISM)
    if (p != null) {
        return min(Runtime.getRuntime().availableProcessors(), p.toInt())
    }
    return Runtime.getRuntime().availableProcessors()
}

fun configureParallelism(parallelism: Int?) {
    if (parallelism != null) {
        require(parallelism.toInt() > 0) {
            "Parallelism level should be at least 1, got: $parallelism"
        }
        System.setProperty(PARALLELISM, parallelism.toString())
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", parallelism.toString())
    }
}

/**
 * Executes tasks re-throwing any exception occurred.
 */
fun ExecutorService.awaitAll(tasks: Iterable<Callable<*>>, ignoreInterrupted: Boolean = false) {
    tasks.map { submit(it) }.toList().forEach {
        try {
            if (ignoreInterrupted) {
                try {
                    it.get()
                } catch (e: InterruptedException) {
                    // Ignore
                }
            } else {
                it.get()
            }
        } catch (e: ExecutionException) {
            shutdownNow()
            throw e.cause ?: e
        }
    }
}

/**
 * We limit number of threads in the pool to obey requested parallelism level
 * IMPORTANT: limited thread pool can cause deadlocks, e.g. when using countdown latch
 * with number of tasks bigger than thread pool size.
 */
fun <T> List<Callable<T>>.await(parallel: Boolean) {
    if (parallel) {
        val executor = Executors.newWorkStealingPool(min(size, parallelismLevel()))
        executor.awaitAll(this)
        check(executor.shutdownNow().isEmpty())
    } else {
        forEach { it.call() }
    }
}