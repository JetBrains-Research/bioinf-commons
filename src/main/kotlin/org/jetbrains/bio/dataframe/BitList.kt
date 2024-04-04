package org.jetbrains.bio.dataframe

import org.jetbrains.bio.genome.Range
import java.util.*
import kotlin.math.min

/**
 * A sibling of [BitSet] which is aware of the universe cardinality.
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 */
class BitList(private val universe: Int) : BitSet() {
    fun copy() = wrap(universe, this)

    /**
     * Iterates over set bits
     */
    operator fun iterator(): IntIterator {
        var ptr: Int
        var nextPtr = nextSetBit(0)

        return object : IntIterator() {
            override fun nextInt(): Int {
                ptr = nextPtr
                nextPtr = nextSetBit(ptr + 1)
                return ptr
            }

            override fun hasNext() = nextPtr >= 0
        }
    }

    override fun set(bitIndex: Int) {
        if (bitIndex >= universe) throw IndexOutOfBoundsException("bitIndex >= size: $bitIndex >= $universe")
        super.set(bitIndex)
    }

    override fun set(fromIndex: Int, toIndex: Int) {
        if (fromIndex >= universe) throw IndexOutOfBoundsException("fromIndex >= size: $fromIndex >= $universe")
        if (toIndex > universe) throw IndexOutOfBoundsException("toIndex > size: $fromIndex > $universe")
        super.set(fromIndex, toIndex)
    }

    override fun clear(bitIndex: Int) {
        if (bitIndex >= universe) throw IndexOutOfBoundsException("bitIndex >= size: $bitIndex >= $universe")
        super.clear(bitIndex)
    }

    override fun clear(fromIndex: Int, toIndex: Int) {
        if (fromIndex >= universe) throw IndexOutOfBoundsException("fromIndex >= size: $fromIndex >= $universe")
        if (toIndex > universe) throw IndexOutOfBoundsException("toIndex > size: $fromIndex > $universe")
        super.clear(fromIndex, toIndex)
    }

    operator fun plus(other: BitList): BitList {
        val acc = copy()
        val prevUniverse = universe
        return wrap(prevUniverse + other.universe, acc).apply {
            for (i in other) {
                set(prevUniverse + i)
            }
        }
    }

    /**
     * Unlike its relative [BitList] returns the universe cardinality.
     *
     *     val bs = BitSet()
     *     bs.set(1)
     *     bs.size()   # == 64 == Long.SIZE
     *
     *     val bts = BitterSet.of(3, bs)
     *     bts.size()  # == 3
     *     bts.set(1)
     *     bts.size()  # == 3
     *
     * Why? Because, concatenation.
     */
    override fun size() = universe

    /**
     * Returns a list of 1-bit runs in this set using gap.
     *
     * For example the following bit set has 2 runs with gap = 0, and 1 run with gap = 1
     *
     *     110111100
     *
     * The first run is [0, 2) and the second [3, 8).
     */
    fun aggregate(gap: Int = 0): List<Range> {
        val ranges = arrayListOf<Range>()
        var offset = 0
        while (offset < size()) {
            val left = nextSetBit(offset)
            if (left == -1) {
                break
            }
            var right = min(nextClearBit(left + 1), size())
            while (gap > 0 && right < size()) {
                val nextSet = nextSetBit(right + 1)
                if (nextSet != -1 && nextSet - right <= gap) {
                    right = min(nextClearBit(nextSet + 1), size())
                } else {
                    break
                }
            }
            ranges.add(Range(left, right))
            offset = right
        }

        return ranges
    }
    override fun equals(other: Any?) = when {
        this === other -> true
        other !is BitList -> false
        else -> universe == other.universe && super.equals(other)
    }

    override fun hashCode() = Objects.hash(super.hashCode(), universe)

    override fun toString() = "$universe@${super.toString()}"

    companion object {
        internal fun wrap(universe: Int, wrapped: BitSet): BitList {
            require(wrapped.cardinality() <= universe)
            return BitList(universe).apply {
                or(wrapped)
            }
        }

        inline operator fun invoke(universe: Int, block: (Int) -> Boolean): BitList {
            return BitList(universe).apply {
                0.until(universe)
                    .filter { block(it) }
                    .forEach(::set)
            }
        }
    }
}

fun BooleanArray.toBitList() = let { arr ->
    BitList(arr.size) { i -> arr[i] }
}

fun IntArray.toBitList() =
    BitList((maxOrNull() ?: 0) + 1) { it in this }
