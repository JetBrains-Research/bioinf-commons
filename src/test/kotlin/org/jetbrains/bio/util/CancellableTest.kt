package org.jetbrains.bio.util

import org.junit.Test
import java.util.concurrent.Callable
import java.util.concurrent.CancellationException
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNull

/**
 * @author Oleg Shpynov
 */
class CancellableTest {
    private val intCallable = Callable<Int> {
        for (i in 0..9) {
            CancellableState.current().checkCanceled()
            Thread.sleep(100)
        }

        42
    }

    private val incorrectCallable = Callable<Int> { throw RuntimeException("WOOT") }

    @Test fun notReady() {
        assertFalse(CancellableTask(intCallable).isDone)
    }

    @Test(expected = CancellationException::class) fun getAfterCancel() {
        val task = CancellableTask(incorrectCallable)
        task.cancel()
        task.get()
    }

    @Test fun cancelledWaitAndGet() {
        val task = CancellableTask(intCallable)
        Thread { task.cancel() }.start()
        assertNull(task.waitAndGet())
    }

    @Test fun testWaitAndGet() {
        assertEquals(42, CancellableTask(intCallable).apply { execute() }.waitAndGet()!!.toInt())
    }

    @Test(expected = RuntimeException::class) fun runtimeExceptionGet() {
        val cancellableTask = CancellableTask(incorrectCallable)
        Thread.sleep(200)
        cancellableTask.get()
    }

    @Test(expected = RuntimeException::class) fun runtimeExceptionWait() {
        CancellableTask(incorrectCallable).apply { execute() }.waitAndGet()
    }
}
