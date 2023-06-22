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
class BitterSet(private val universe: Int) : BitSet() {
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

    operator fun plus(other: BitterSet): BitterSet {
        val acc = copy()
        for (i in other) {
            acc.set(universe + i)
        }
        return wrap(universe + other.universe, acc)
    }

    /**
     * Unlike its relative [BitterSet] returns the universe cardinality.
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
     * Returns a list of 1-bit consequent blocks.
     *
     * For example the following bit set has 2 blocks
     *
     *     110111100
     *
     * The first run is [0, 2) and the second [3, 8).
     */
    fun aggregate(): List<Range> {
        val ranges = arrayListOf<Range>()
        var offset = 0
        while (offset < size()) {
            val left = nextSetBit(offset)
            if (left == -1) {
                break
            }
            val right = min(nextClearBit(left + 1), size())
            ranges.add(Range(left, right))
            offset = right
        }

        return ranges
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is BitterSet -> false
        else -> universe == other.universe && super.equals(other)
    }

    override fun hashCode() = Objects.hash(super.hashCode(), universe)

    override fun toString() = "$universe@${super.toString()}"

    companion object {
        internal fun wrap(universe: Int, wrapped: BitSet): BitterSet {
            require(wrapped.cardinality() <= universe)
            return BitterSet(universe).apply {
                or(wrapped)
            }
        }

        inline operator fun invoke(universe: Int, block: (Int) -> Boolean): BitterSet {
            return BitterSet(universe).apply {
                0.until(universe)
                    .filter { block(it) }
                    .forEach(::set)
            }
        }
    }
}

fun BooleanArray.toBitterSet() = let { arr ->
    BitterSet(arr.size) { i -> arr[i] }
}

fun IntArray.toBitterSet() =
    BitterSet((maxOrNull() ?: 0) + 1) { it in this }
