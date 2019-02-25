package org.jetbrains.bio.util

import org.junit.Test
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import kotlin.test.assertEquals

private object OhNoesException : Exception()

class ExecutorExtensionsTest {
    @Test(expected = OhNoesException::class)
    fun testAwaitAll() {
        val tasks = (1..4).map {
            Callable {
                if (it % 2 == 0) {
                    throw OhNoesException
                }
            }
        }

        val executor = Executors.newCachedThreadPool()
        executor.awaitAll(tasks)
        check(executor.shutdownNow().isEmpty())
    }

    @Test
    fun testSingleLaunch() {
        val list = Collections.synchronizedList(arrayListOf<Int>())
        val tasks = (1..3).map { n ->
            Callable {
                list.add(n)
            }
        }
        val executor = Executors.newWorkStealingPool(parallelismLevel())
        executor.awaitAll(tasks)
        assertEquals(3, list.size)
    }

}