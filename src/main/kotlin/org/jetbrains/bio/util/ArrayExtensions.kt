package org.jetbrains.bio.dataframe

/**
 * Manually copy-pasted routines for primitive arrays.
 *
 * @author Sergei Lebedev
 */

import com.google.common.collect.ComparisonChain
import java.util.stream.IntStream

/**
 * Returns the ordering of the original array which makes it sorted.
 *
 * If `res` is the result of method call, then for any appropriate
 *
 *   i < j => values[res[i]] <= values[res[j]]
 *
 * E.g., `values[res[0]]` is the minimum of the array.
 *
 * @param reverse sort into descending order.
 */
fun ByteArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { IntIntPair(this[it].toInt(), it) }
        .sorted(IntIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

fun ShortArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { IntIntPair(this[it].toInt(), it) }
        .sorted(IntIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

fun IntArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { IntIntPair(this[it], it) }
        .sorted(IntIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

fun LongArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { LongIntPair(this[it], it) }
        .sorted(LongIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

fun FloatArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { DoubleIntPair(this[it].toDouble(), it) }
        .sorted(DoubleIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

fun DoubleArray.argSort(reverse: Boolean = false): IntArray {
    return IntStream.range(0, size)
        .mapToObj { DoubleIntPair(this[it], it) }
        .sorted(DoubleIntPair.comparator(reverse))
        .mapToInt { it.second }
        .toArray()
}

private class IntIntPair(val first: Int, val second: Int) {
    companion object {
        fun comparator(reverse: Boolean): Comparator<IntIntPair> {
            return if (reverse) COMPARATOR.reversed() else COMPARATOR
        }

        private val COMPARATOR: Comparator<IntIntPair> = Comparator { p1, p2 ->
            ComparisonChain.start()
                .compare(p1.first, p2.first)
                .compare(p1.second, p2.second)
                .result()
        }
    }
}

private class LongIntPair(val first: Long, val second: Int) {
    companion object {
        fun comparator(reverse: Boolean): Comparator<LongIntPair>? {
            return if (reverse) COMPARATOR.reversed() else COMPARATOR
        }

        private val COMPARATOR: Comparator<LongIntPair> = Comparator { p1, p2 ->
            ComparisonChain.start()
                .compare(p1.first, p2.first)
                .compare(p1.second, p2.second)
                .result()
        }
    }
}

private class DoubleIntPair(val first: Double, val second: Int) {
    companion object {
        fun comparator(reverse: Boolean): Comparator<DoubleIntPair>? {
            return if (reverse) COMPARATOR.reversed() else COMPARATOR
        }

        private val COMPARATOR: Comparator<DoubleIntPair> = Comparator { p1, p2 ->
            ComparisonChain.start()
                .compare(p1.first, p2.first)
                .compare(p1.second, p2.second)
                .result()
        }
    }
}

/**
 * Reorders the array in place using the given ordering, so that
 * `values[i]` after the method call is equal to `values[indices[i]]`
 * before the method call.
 *
 * When `indices` are the result of `values.argSort()`, the method
 * makes the array sorted.
 */
fun ByteArray.reorder(indices: IntArray) {
    require(size == indices.size) { "non-conformable arrays" }
    val copy = indices.clone()
    for (i in indices) {
        val value = this[i]
        var j = i
        while (true) {
            val k = copy[j]
            copy[j] = j
            if (k == i) {
                this[j] = value
                break
            } else {
                this[j] = this[k]
                j = k
            }
        }
    }
}

fun ShortArray.reorder(indices: IntArray) {
    require(size == indices.size) { "non-conformable arrays" }
    val copy = indices.clone()
    for (i in indices) {
        val value = this[i]
        var j = i
        while (true) {
            val k = copy[j]
            copy[j] = j
            if (k == i) {
                this[j] = value
                break
            } else {
                this[j] = this[k]
                j = k
            }
        }
    }
}

fun IntArray.reorder(indices: IntArray) {
    require(size == indices.size) { "non-conformable arrays" }
    val copy = indices.clone()
    for (i in indices) {
        val value = this[i]
        var j = i
        while (true) {
            val k = copy[j]
            copy[j] = j
            if (k == i) {
                this[j] = value
                break
            } else {
                this[j] = this[k]
                j = k
            }
        }
    }
}

fun LongArray.reorder(indices: IntArray) {
    require(size == indices.size) { "non-conformable arrays" }
    val copy = indices.clone()
    for (i in indices) {
        val value = this[i]
        var j = i
        while (true) {
            val k = copy[j]
            copy[j] = j
            if (k == i) {
                this[j] = value
                break
            } else {
                this[j] = this[k]
                j = k
            }
        }
    }
}

fun DoubleArray.reorder(indices: IntArray) {
    require(size == indices.size) { "non-conformable arrays" }
    val copy = indices.clone()
    for (i in indices) {
        val value = this[i]
        var j = i
        while (true) {
            val k = copy[j]
            copy[j] = j
            if (k == i) {
                this[j] = value
                break
            } else {
                this[j] = this[k]
                j = k
            }
        }
    }
}

operator fun ShortArray.div(other: ShortArray): FloatArray {
    require(size == other.size) { "non-conformable arrays" }
    val res = FloatArray(size)
    for (i in res.indices) {
        res[i] = this[i].toFloat() / other[i]
    }

    return res
}