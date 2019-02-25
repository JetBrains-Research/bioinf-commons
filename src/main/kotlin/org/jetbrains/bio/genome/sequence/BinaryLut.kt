package org.jetbrains.bio.genome.sequence

import java.util.*

/**
 * A lookup table reducing the search space for binary search
 * over a sorted array of non-negative integers.
 *
 * See https://geidav.wordpress.com/2013/12/29/optimizing-binary-search.
 *
 * @author Sergei Lebedev
 */
open class BinaryLut(private val index: IntArray,
                     /** The number of bits to use for LUT indexing. */
                     private val bits: Int,
                     /**
                      * The last LUT boundary. Keys larger than this should
                      * have the `toIndex` equal to the total number of
                      * elements indexed.
                      */
                     private val end: Int) {

    open fun binarySearch(data: IntArray, key: Int): Int {
        if (key < 0) return -1 // a negative key is strictly less than any data element
        val idx = key ushr Integer.SIZE - bits
        val from = index[idx]
        val to = if (idx + 1 > end) data.size else index[idx + 1]
        return Arrays.binarySearch(data, from, to, key)
    }

    /**
     * Returns an element having minimal distance to "key", or -1 when "data" is empty.
     * If there are several such elements, an arbitrary one is returned.
     */
    fun nearestElemIdx(data: IntArray, key: Int): Int {
        val i = binarySearch(data, key)
        if (i >= 0) {
            return i
        }
        val idx = i.inv()
        return when (idx) {
            data.size -> idx - 1
            0 -> idx
            else -> if (key - data[idx - 1] <= data[idx] - key) idx - 1 else idx
        }
    }

    /**
     * Returns an inclusive range of indices of elements having minimal distance to "key",
     * or (-1,-1) when "data" is empty. There are two reasons for this range being non-singular:
     * 1. equidistant elements, e.g. searching for an element nearest to 2 in an array {0, 1, 3, 4} will return (1,2)
     * 2. multiple elements, e.g. searching for an element nearest to 2 in an array {0, 2, 2, 4} will return (1,2)
     * We expect to have a few of these cases, so after we find a match the linear search is implemented.
     */
    fun nearestElemDist(data: IntArray, key: Int): Pair<Int, Int> {
        val idx = nearestElemIdx(data, key)
        if (idx == -1) return Pair(-1, -1)
        val dist = Math.abs(data[idx] - key)
        val lowVal = key - dist
        val highVal = key + dist
        val low = run {
            var i = idx
            while (i > 0 && (data[i - 1] == lowVal || data[i - 1] == highVal)) {
                i--
            }
            i
        }
        val high = run {
            var i = idx
            while (i < data.size - 1 && (data[i + 1] == lowVal || data[i + 1] == highVal)) {
                i++
            }
            i
        }
        return Pair(low, high)
    }

    /**
     * Returns an inclusive range of indices of elements framing "key" both on the left and on the right side,
     * or (-1,-1) when "data" is empty. There can be multiple framing elements on each side, since they can be equal.
     * Also, the exact matches are included in the answer, but aren't counted as framing
     * (so even when "key" is present in "data", we look to the sides all the same).
     */
    fun nearestElemLR(data: IntArray, key: Int): Pair<Int, Int> {
        if (data.isEmpty()) return -1 to -1
        val idx = binarySearch(data, key)
        var low: Int
        var high: Int
        if (idx >= 0) {
            low = run {
                var i = idx - 1
                while (i >= 0 && data[i] == key) { i-- }
                i
            }
            high = run {
                var i = idx + 1
                while (i < data.size && data[i] == key) { i++ }
                i
            }
        } else {
            low = idx.inv() - 1
            high = idx.inv()
        }
        if (low < 0) {
            low = 0
        } else {
            while (low > 0 && data[low] == data[low - 1]) { low-- }
        }
        if (high == data.size) {
            high = data.size - 1
        } else {
            while (high < data.size - 1 && data[high] == data[high + 1]) { high++ }
        }
        return low to high
    }

    /**
     * Returns an inclusive range of indices of elements framing "key" both on the left and on the right side
     * and matching the predicate, or (-1,-1) when no elements match the predicate.
     * There can be multiple framing elements on each side, since they can be equal.
     * If an exact match (also matching the predicate) is present, the resulting range will contain only exact matches.
     * The resulting range can contain some elements not matching the predicate.
     */
    fun framingElem(data: IntArray, key: Int, predicate: (Int) -> Boolean): List<Int> {
        val idx = binarySearch(data, key)
        val left = run {
            var i = if (idx >= 0) idx else idx.inv() - 1
            while (i >= 0 && !predicate(i)) { i-- }
            i
        }
        val right = run {
            var i = if (idx >= 0) idx else idx.inv()
            while (i < data.size && !predicate(i)) { i++ }
            i
        }
        if (left == -1 && right == data.size) return emptyList()
        var low: Int
        var high: Int
        if (left >= 0) {
            low = left
            val leftVal = data[left]
            var i = left - 1
            while (i >= 0 && data[i] == leftVal) {
                if (predicate(i)) low = i
                i--
            }
        } else {
            low = right
        }
        if (right != data.size) {
            high = right
            val rightVal = data[right]
            var i = right + 1
            while (i < data.size && data[i] == rightVal) {
                if (predicate(i)) high = i
                i++
            }
        } else {
            high = left
        }
        return (low..high).filter(predicate)
    }

    /**
     * Returns an inclusive range of indices of elements that are at a given distance or less from "key"
     * both on the left and on the right side, or (-1,-1) when there are no such elements.
     */
    fun elemWithinDist(data: IntArray, key: Int, left: Int, right: Int = left): Pair<Int, Int> {
        if (data.isEmpty()) return -1 to -1
        val idxLeft = binarySearch(data, key - left)
        val idxRight = binarySearch(data, key + right)
        val low = run {
            var i = if (idxLeft >= 0) idxLeft else idxLeft.inv()
            while (i > 0 && data[i - 1] == key - left) { i-- }
            i
        }
        val high = run {
            var i = if (idxRight >= 0) idxRight else idxRight.inv() - 1
            while (i < data.size - 1 && data[i + 1] == key + right) { i++ }
            i
        }
        return if (low <= high) {
            low to high
        } else {
            -1 to -1
        }
    }


    companion object {
        fun of(data: IntArray, bits: Int): BinaryLut {
            require(bits < Integer.SIZE) { "bits must be <${Integer.SIZE}" }
            if (data.isEmpty()) {
                return EmptyBinaryLut(bits)
            }

            // Invariant: index[key(x)] = i, s.t. data[i] <= x
            val index = IntArray((1 shl bits) + 1)
            var bound = 0
            var ptr = 0
            for (i in 0..data.size - 1) {
                val nextBound = data[i] ushr Integer.SIZE - bits
                index[bound] = ptr

                if (nextBound > bound) {
                    ptr = i
                    index.fill(ptr, bound + 1, nextBound)
                }

                bound = nextBound
            }

            index.fill(ptr, bound, index.size)
            return BinaryLut(index, bits, bound)
        }
    }
}

private class EmptyBinaryLut(bits: Int) : BinaryLut(IntArray(0), bits, 0) {
    override fun binarySearch(data: IntArray, key: Int) = -1
}
