package org.jetbrains.bio.dataframe

import java.util.function.IntPredicate

/**
 * @author Roman Chernyatchik
 */

interface RowPredicateFactory {
    operator fun invoke(df: DataFrame): IntPredicate

    companion object {
        inline operator fun invoke(crossinline f: (DataFrame) -> IntPredicate) = object : RowPredicateFactory {
            override fun invoke(df: DataFrame) = f(df)
        }
    }
}

inline fun byByte(label: String, crossinline p: (Byte) -> Boolean) = RowPredicateFactory { df ->
    val data = df.sliceAsByte(label)
    IntPredicate { row -> p(data[row]) }
}

inline fun byShort(label: String, crossinline p: (Short) -> Boolean) = RowPredicateFactory { df ->
    val data = df.sliceAsShort(label)
    IntPredicate { row -> p(data[row]) }
}

inline fun byInt(label: String, crossinline p: (Int) -> Boolean) = RowPredicateFactory { df ->
    val col = df.sliceAsInt(label)
    IntPredicate { row -> p(col[row]) }
}

inline fun byLong(label: String, crossinline p: (Long) -> Boolean) = RowPredicateFactory { df ->
    val col = df.sliceAsLong(label)
    IntPredicate { row -> p(col[row]) }
}

inline fun byFloat(label: String, crossinline p: (Float) -> Boolean) = RowPredicateFactory { df ->
    val data = df.sliceAsFloat(label)
    IntPredicate { row -> p(data[row]) }
}

inline fun byDouble(label: String, crossinline p: (Double) -> Boolean) = RowPredicateFactory { df ->
    val data = df.sliceAsDouble(label)
    IntPredicate { row -> p(data[row]) }
}

@Suppress("nothing_to_inline")
inline fun byBool(label: String) = RowPredicateFactory { df ->
    val data = df.sliceAsBool(label)
    IntPredicate { row -> data[row] }
}

@Suppress("nothing_to_inline")
inline fun byNotBool(label: String) = RowPredicateFactory { df ->
    val data = df.sliceAsBool(label)
    IntPredicate { row -> !data[row] }
}

inline fun <T> byObj(label: String, crossinline p: (T) -> Boolean) = RowPredicateFactory { df ->
    val data = df.sliceAsObj<T>(label)
    IntPredicate { row -> p(data[row]) }
}

inline fun byString(label: String, crossinline p: (String) -> Boolean) = byObj(label, p)

inline fun predicate(crossinline rowPredicate: DataFrame.(Int) -> Boolean) = RowPredicateFactory { df -> IntPredicate { row -> df.rowPredicate(row) } }

fun all(pf: RowPredicateFactory, vararg rest: RowPredicateFactory): RowPredicateFactory {
    val pfs = arrayOf(pf) + rest
    return when (pfs.size) {
        1 -> pfs.single()
        2 -> RowPredicateFactory { df ->
            val p0 = pfs[0](df)
            val p1 = pfs[1](df)
            IntPredicate { p0.test(it) && p1.test(it) }
        }
        3 -> RowPredicateFactory { df ->
            val p0 = pfs[0](df)
            val p1 = pfs[1](df)
            val p2 = pfs[2](df)
            IntPredicate { p0.test(it) && p1.test(it) && p2.test(it) }
        }
        else -> RowPredicateFactory { df ->
            // XXX we could optimize it by processing pfs in threes.
            val ps = pfs.map { it(df) }
            IntPredicate { row ->
                for (p in ps) {
                    if (!p.test(row)) {
                        return@IntPredicate false
                    }
                }

                true
            }
        }
    }
}
