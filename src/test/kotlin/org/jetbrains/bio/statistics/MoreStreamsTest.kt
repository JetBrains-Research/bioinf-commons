package org.jetbrains.bio.statistics

import com.google.common.collect.Sets
import com.google.common.math.IntMath
import org.apache.commons.math3.distribution.PoissonDistribution
import org.jetbrains.bio.Retry
import org.jetbrains.bio.RetryRule
import org.junit.Rule
import org.junit.Test
import java.util.concurrent.ForkJoinPool
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class MoreStreamsTest {

    @get:Rule
    var rule = RetryRule(3)

    // the test checks empty range handling, so the empty range is intended
    @Suppress("EmptyRange")
    @Test
    fun chunkedStreamEmpty() {
        assertEquals(
            emptyList(),
            (1..0).chunked().iterator().asSequence().toList()
        )
        assertEquals(
            emptyList(),
            (1..0).chunked(8).iterator().asSequence().toList()
        )
    }

    @Test
    fun chunkedStreamSingleton() {
        assertEquals(
            listOf(Chunk(0, 1)),
            (0..0).chunked().iterator().asSequence().toList()
        )
        assertEquals(
            listOf(Chunk(0, 1)),
            (0..0).chunked(8).iterator().asSequence().toList()
        )
    }

    @Test
    fun chunkedStreamMultiple() {
        val chunks = (0..2).chunked().iterator().asSequence().toList()
        when (ForkJoinPool.getCommonPoolParallelism()) {
            1 -> assertEquals(listOf(Chunk(0, 3)), chunks)
            2 -> assertEquals(listOf(Chunk(0, 2), Chunk(2, 3)), chunks)
            else -> assertEquals(listOf(Chunk(0, 1), Chunk(1, 2), Chunk(2, 3)), chunks)
        }
    }

    @Test
    fun chunkedStreamSum() {
        val numObservations = IntMath.pow(2, 16)
        val values = PoissonDistribution(42.0).sample(numObservations)
        val sum = (0 until numObservations).chunked(2).mapToLong { chunk ->
            var acc: Long = 0
            for (i in chunk.lo until chunk.hi) {
                acc += values[i].toLong()
            }

            acc
        }.sum()

        var acc: Long = 0
        for (value in values) {
            acc += value.toLong()
        }

        assertEquals(acc, sum)
    }

    @Retry
    @Test
    fun chunkedStreamIsParallel() {
        val numObservations = IntMath.pow(2, 16)
        val chunked = (0 until numObservations).chunked(4)
        val threadsInvolved = chunked.map {
            setOf(Thread.currentThread().id)
        }.reduce(emptySet()) { a, b -> Sets.union(a, b) }
        // We cannot be 100% sure that all the available threads will be used,
        // check range to increase chances
        if (ForkJoinPool.getCommonPoolParallelism() > 1) {
            assertTrue(threadsInvolved.size in 2..ForkJoinPool.getCommonPoolParallelism() + 1) // + main thread
        }
    }
}
