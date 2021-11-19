package org.jetbrains.bio.dataframe

import java.util.*

/**
 * A container for typed single-column data.
 *
 * @author Sergei Lebedev
 */
abstract class Column<T> protected constructor(
    val label: String, val data: T, val boxedType: Class<out Any>
) {

    // Note(lebedev): I'm not sure if this is the best way to perform
    // T->double coercion. Any suggestions -- let me know.
    abstract fun getAsDouble(row: Int): Double

    /**
     * Constructs a column with the supplied `data`.
     */
    abstract fun wrap(newData: T): Column<T>

    /**
     * Constructs a column holding the same `data` but with the supplied
     * `label`.
     */
    abstract fun rename(newLabel: String): Column<T>

    /**
     * Shuffles column `data` according to a given permutation of indices.
     */
    abstract fun reorder(indices: IntArray): Column<T>

    /**
     * Alters column size.
     *
     * @see Arrays.copyOf methods for semantics.
     */
    abstract fun resize(newSize: Int): Column<T>

    /** Returns a new column with a copy of the underlying array. */
    fun copy() = resize(size)

    /**
     * Concatenates two columns of the same type and returns a new column.
     */
    abstract operator fun plus(other: Column<*>): Column<T>

    /**
     * Similar to [plus], but only adds the values not yet present in `data`.
     */
    abstract fun merge(other: Column<*>): Column<T>

    /**
     * Filters rows based on a given `mask`. The resulting column
     * only contains the rows corresponding to non-zero bits.
     */
    abstract fun filter(mask: BitSet): Column<T>

    /**
     * Intersects a column with another column.
     *
     * @param c column.
     * @return a predicate which for a given array and index outputs
     *         `true` if an element belongs to the intersection
     *         and `false` otherwise.
     */
    internal abstract fun intersect(c: Column<*>): ObjIntPredicate<T>

    /**
     * Returns a sorted permutation of indices for a given column.
     */
    abstract fun sorted(reverse: Boolean = false): IntArray

    abstract val size: Int

    internal abstract fun dump(row: Int): String

    // Note(lebedev): I don't like this method, because it mutates
    // the column internal state.
    abstract fun load(row: Int, value: String)

    internal open fun typeName() = boxedType.simpleName

    @Suppress("unchecked_cast")
    internal fun test(p: ObjIntPredicate<*>): BitSet {
        val bs = BitSet(size)
        for (i in 0 until size) {
            if ((p as ObjIntPredicate<T>).test(data, i)) {
                bs.set(i)
            }
        }

        return bs
    }
}

interface ObjIntPredicate<T> {
    /**
     * Evaluates this predicate on the given arguments.
     *
     * @param t the first input argument
     * @param value the second input argument
     * @return `true` if the input arguments match the predicate,
     *         otherwise `false`
     */
    fun test(t: T, value: Int): Boolean

    infix fun and(other: ObjIntPredicate<in T>): ObjIntPredicate<T> {
        return invoke { t, value -> test(t, value) && other.test(t, value) }
    }

    companion object {
        // Syntax sugar for instantiation like: ObjIntPredicate { data, i -> ... }
        inline operator fun <T> invoke(crossinline f: (T, Int) -> Boolean): ObjIntPredicate<T> {
            return object : ObjIntPredicate<T> {
                override fun test(t: T, value: Int) = f(t, value)
            }
        }
    }
}
