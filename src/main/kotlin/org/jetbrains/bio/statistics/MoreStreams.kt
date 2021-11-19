package org.jetbrains.bio.statistics

import com.google.common.primitives.Ints
import java.util.*
import java.util.concurrent.ForkJoinPool
import java.util.concurrent.ForkJoinTask
import java.util.function.Consumer
import java.util.stream.IntStream
import java.util.stream.Stream
import java.util.stream.StreamSupport

data class Chunk(val lo: Int, val hi: Int)

private data class ChunkSpliterator(
    /** Start index (inclusive). */
    private var from: Int,
    /** End index (exclusive). */
    private val to: Int,
    /** Maximum number of elements in a chunk. */
    private val chunkSize: Int
) : Spliterator<Chunk> {

    override fun tryAdvance(action: Consumer<in Chunk>): Boolean {
        if (from >= to) {
            return false
        }

        val chunkTo = Math.min(from + chunkSize, to)
        action.accept(Chunk(from, chunkTo))
        from = chunkTo
        return true
    }

    override fun trySplit(): Spliterator<Chunk>? {
        return if (to - from <= chunkSize) {
            null
        } else {
            val newFrom = from + (to - from) / 2
            val split = ChunkSpliterator(from, newFrom, chunkSize)
            from = newFrom
            split
        }
    }

    override fun estimateSize(): Long {
        return ((to - from + chunkSize - 1) / chunkSize).toLong()
    }

    override fun characteristics(): Int {
        return (Spliterator.ORDERED or Spliterator.SORTED or
                Spliterator.SIZED or Spliterator.SUBSIZED or
                Spliterator.NONNULL or Spliterator.IMMUTABLE or
                Spliterator.DISTINCT)
    }

    override fun getComparator(): Comparator<in Chunk> {
        return Comparator { o1, o2 -> Ints.compare(o1.lo, o2.lo) }
    }
}

/**
 * Creates a parallel stream which operates over chunks of a given range.
 * Think a chunk-parallel for-loop.
 *
 * Internally we use `ChunkSpliterator` which unlike the spliterator
 * used in `IntStream` does not assume low computational load on
 * each int.
 *
 * @param perCore number of chunks to be delegated to each available core.
 */
fun IntRange.chunked(perCore: Int = 1): Stream<Chunk> {
    if (start > endInclusive) {
        return Stream.empty()  // Just like a forloop.
    }

    // I'm sure there's a better way to handle #chunks > range length
    // situation, but at least this code doesn't throw.
    val chunkCount = ForkJoinPool.getCommonPoolParallelism() * perCore
    val chunkSize = Math.max((endInclusive + 1 - start + chunkCount - 1) / chunkCount, 1)
    return StreamSupport.stream(ChunkSpliterator(start, endInclusive + 1, chunkSize), true)
}

/**
 * Converts an [IntRange] to a parallel [IntStream].
 */
@Suppress("nothing_to_inline")
inline fun IntRange.parallel(): IntStream = IntStream.range(start, endInclusive + 1).parallel()

/**
 * A parallel for-loop over ForkJoin pool.
 *
 * Unlike [java.util.stream.IntStream] it doesn't assume a single
 * iteration is cheap and creates a [ForkJoinTask] per iteration.
 *
 * Avoid using [forking] for long blocking operations, e.g. I/O.
 */
inline fun IntRange.forking(crossinline block: (Int) -> Unit) {
    val tasks = Array(endInclusive - start + 1) {
        ForkJoinTask.adapt { block(start + it) }
    }

    ForkJoinTask.invokeAll(*tasks)
}