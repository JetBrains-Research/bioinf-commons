@file:Suppress("platform_class_mapped_to_kotlin")

package org.jetbrains.bio.dataframe

import com.google.common.base.MoreObjects
import com.google.common.primitives.*
import gnu.trove.list.array.*
import gnu.trove.set.hash.*
import org.jetbrains.bio.util.Logs
import org.slf4j.LoggerFactory
import java.util.*

private val LOG = LoggerFactory.getLogger("org.jetbrains.bio.dataframe.PrimitiveColumns")

// XXX this is done because kotlin.Byte (and other primitives)
// resolve to a _boxed_ java.lang.Byte, instead of a primitive
// class.
class ByteColumn(
    label: String, data: ByteArray,
    val valueFormatter: ((Byte) -> String)? = null
) : Column<ByteArray>(label, data, java.lang.Byte::class.java) {

    override fun getAsDouble(row: Int) = data[row].toDouble()

    override fun rename(newLabel: String) = ByteColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: ByteArray) = ByteColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<ByteArray> {
        return wrap(Bytes.concat(data, other.data as ByteArray))
    }

    override fun merge(other: Column<*>): Column<ByteArray> {
        val seen = TByteHashSet(data)
        check(seen.size() == data.size) {
            check(seen.size() == data.size) {
                val seenElements = TByteHashSet(data.size)
                val duplicates = arrayListOf<Pair<Int, Byte>>()
                for ((idx, e) in data.withIndex()) {
                    if (e in seenElements) {
                        duplicates.add(idx to e)
                        if (duplicates.size >= 10) {
                            break
                        }
                    } else {
                        seenElements.add(e)
                    }
                }
                "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                        " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
            }
        }

        val merged = TByteArrayList(data)
        for (value in other.data as ByteArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<ByteArray> {
        val copy = ByteArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    override fun intersect(c: Column<*>): ObjIntPredicate<ByteArray> {
        val other = c as ByteColumn
        val seen = TByteHashSet(data)
        seen.retainAll(other.data.clone())
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<ByteArray> {
        require(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int) = wrap(data.copyOf(newSize))

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toByte()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is ByteColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class ShortColumn(
    label: String,
    data: ShortArray,
    val valueFormatter: ((Short) -> String)? = null
) : Column<ShortArray>(label, data, java.lang.Short::class.java) {

    override fun getAsDouble(row: Int) = data[row].toDouble()

    override fun rename(newLabel: String) = ShortColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: ShortArray) = ShortColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<ShortArray> {
        return wrap(Shorts.concat(data, other.data as ShortArray))
    }

    override fun merge(other: Column<*>): Column<ShortArray> {
        val seen = TShortHashSet(data)
        check(seen.size() == data.size) {
            check(seen.size() == data.size) {
                val seenElements = TShortHashSet(data.size)
                val duplicates = arrayListOf<Pair<Int, Short>>()
                for ((idx, e) in data.withIndex()) {
                    if (e in seenElements) {
                        duplicates.add(idx to e)
                        if (duplicates.size >= 10) {
                            break
                        }
                    } else {
                        seenElements.add(e)
                    }
                }
                "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                        " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
            }
        }

        val merged = TShortArrayList(data)
        for (value in other.data as ShortArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<ShortArray> {
        val copy = ShortArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<ShortArray> {
        require(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int): Column<ShortArray> = wrap(data.copyOf(newSize))

    override fun intersect(c: Column<*>): ObjIntPredicate<ShortArray> {
        val other = c as ShortColumn
        val seen = TShortHashSet(data)
        seen.retainAll(other.data)
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toShort()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is ShortColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class IntColumn(
    label: String,
    data: IntArray,
    val valueFormatter: ((Int) -> String)? = null
) : Column<IntArray>(label, data, java.lang.Integer::class.java) {

    override fun getAsDouble(row: Int) = data[row].toDouble()

    override fun rename(newLabel: String) = IntColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: IntArray) = IntColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<IntArray> {
        return wrap(Ints.concat(data, other.data as IntArray))
    }

    override fun merge(other: Column<*>): Column<IntArray> {
        val seen = TIntHashSet(data)
        check(seen.size() == data.size) {
            val seenElements = TIntHashSet(data.size)
            val duplicates = arrayListOf<Pair<Int, Int>>()
            for ((idx, e) in data.withIndex()) {
                if (e in seenElements) {
                    duplicates.add(idx to e)
                    if (duplicates.size >= 10) {
                        break
                    }
                } else {
                    seenElements.add(e)
                }
            }
            "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                    " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
        }

        val merged = TIntArrayList(data)
        for (value in other.data as IntArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<IntArray> {
        val copy = IntArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<IntArray> {
        check(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int): Column<IntArray> {
        if (newSize < 0) {
            Logs.getRootLogger()
            LOG.error("Current size: $size, new size = $newSize")
        }
        return wrap(data.copyOf(newSize))
    }

    override fun intersect(c: Column<*>): ObjIntPredicate<IntArray> {
        val other = c as IntColumn
        val seen = TIntHashSet(data)
        seen.retainAll(other.data.clone())
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toInt()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is IntColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class LongColumn(
    label: String, data: LongArray,
    val valueFormatter: ((Long) -> String)? = null
) : Column<LongArray>(label, data, java.lang.Long::class.java) {

    override fun getAsDouble(row: Int) = data[row].toDouble()

    override fun rename(newLabel: String) = LongColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: LongArray) = LongColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<LongArray> {
        return wrap(Longs.concat(data, other.data as LongArray))
    }

    override fun merge(other: Column<*>): Column<LongArray> {
        val seen = TLongHashSet(data)
        check(seen.size() == data.size) {
            check(seen.size() == data.size) {
                val seenElements = TLongHashSet(data.size)
                val duplicates = arrayListOf<Pair<Int, Long>>()
                for ((idx, e) in data.withIndex()) {
                    if (e in seenElements) {
                        duplicates.add(idx to e)
                        if (duplicates.size >= 10) {
                            break
                        }
                    } else {
                        seenElements.add(e)
                    }
                }
                "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                        " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
            }
        }

        val merged = TLongArrayList(data)
        for (value in other.data as LongArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<LongArray> {
        val copy = LongArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    @Deprecated("")
    override fun intersect(c: Column<*>): ObjIntPredicate<LongArray> {
        // BitSet doesn't allow long indices.
        throw UnsupportedOperationException()
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<LongArray> {
        require(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int): Column<LongArray> = wrap(data.copyOf(newSize))

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toLong()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is LongColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class DoubleColumn(
    label: String, data: DoubleArray,
    val valueFormatter: ((Double) -> String)? = null
) : Column<DoubleArray>(label, data, java.lang.Double::class.java) {

    override fun getAsDouble(row: Int) = data[row]

    override fun rename(newLabel: String) = DoubleColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: DoubleArray) = DoubleColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<DoubleArray> {
        return wrap(Doubles.concat(data, other.data as DoubleArray))
    }

    override fun merge(other: Column<*>): Column<DoubleArray> {
        val seen = TDoubleHashSet(data)
        check(seen.size() == data.size) {
            check(seen.size() == data.size) {
                val seenElements = TDoubleHashSet(data.size)
                val duplicates = arrayListOf<Pair<Int, Double>>()
                for ((idx, e) in data.withIndex()) {
                    if (e in seenElements) {
                        duplicates.add(idx to e)
                        if (duplicates.size >= 10) {
                            break
                        }
                    } else {
                        seenElements.add(e)
                    }
                }
                "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                        " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
            }
        }

        val merged = TDoubleArrayList(data)
        for (value in other.data as DoubleArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<DoubleArray> {
        val copy = DoubleArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<DoubleArray> {
        require(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int) = wrap(data.copyOf(newSize))

    override fun intersect(c: Column<*>): ObjIntPredicate<DoubleArray> {
        val other = c as DoubleColumn
        val seen = TDoubleHashSet(data)
        seen.retainAll(other.data.clone())
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toDouble()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is DoubleColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class FloatColumn(
    label: String,
    data: FloatArray,
    val valueFormatter: ((Float) -> String)? = null
) : Column<FloatArray>(label, data, java.lang.Float::class.java) {

    override fun getAsDouble(row: Int) = data[row].toDouble()

    override fun rename(newLabel: String) = FloatColumn(newLabel, data, valueFormatter)

    override fun wrap(newData: FloatArray) = FloatColumn(label, newData, valueFormatter)

    override fun plus(other: Column<*>): Column<FloatArray> {
        return wrap(Floats.concat(data, other.data as FloatArray))
    }

    override fun merge(other: Column<*>): Column<FloatArray> {
        val seen = TFloatHashSet(data)
        check(seen.size() == data.size) {
            check(seen.size() == data.size) {
                val seenElements = TFloatHashSet(data.size)
                val duplicates = arrayListOf<Pair<Int, Float>>()
                for ((idx, e) in data.withIndex()) {
                    if (e in seenElements) {
                        duplicates.add(idx to e)
                        if (duplicates.size >= 10) {
                            break
                        }
                    } else {
                        seenElements.add(e)
                    }
                }
                "duplicates found in '$label' column: unique items: ${seenElements.size()}, column size: ${data.size}" +
                        " Duplicates: ${duplicates.joinToString { (i, e) -> "($i: $e)" }}"
            }
        }

        val merged = TFloatArrayList(data)
        for (value in other.data as FloatArray) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toArray())
    }

    override fun filter(mask: BitSet): Column<FloatArray> {
        val copy = FloatArray(mask.cardinality())
        var offset = 0
        for (i in data.indices) {
            if (mask[i]) {
                copy[offset++] = data[i]
            }
        }

        return wrap(copy)
    }

    override fun sorted(reverse: Boolean) = data.argSort(reverse)

    override fun reorder(indices: IntArray): Column<FloatArray> {
        require(data.size == indices.size)
        val clone = data.clone()
        for (i in 0 until size) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int): Column<FloatArray> = wrap(data.copyOf(newSize))

    override fun intersect(c: Column<*>): ObjIntPredicate<FloatArray> {
        val other = c as FloatColumn
        val seen = TFloatHashSet(data)
        seen.retainAll(other.data.clone())
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override val size: Int get() = data.size

    override fun load(row: Int, value: String) {
        data[row] = value.toFloat()
    }

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: value.toString()
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is FloatColumn -> false
        else -> label == other.label && data.contentEquals(other.data)
    }

    override fun hashCode() = arrayOf(label, data).contentDeepHashCode()

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data.contentToString())
        .toString()
}

class BooleanColumn(
    label: String,
    data: BitList,
    val valueFormatter: ((Boolean) -> String)? = null
) : Column<BitList>(label, data, java.lang.Boolean::class.java) {

    override fun load(row: Int, value: String) {
        data[row] = value == "1" || value.equals("true", ignoreCase = true)
    }

    override fun getAsDouble(row: Int) = if (data[row]) 1.0 else 0.0

    override fun wrap(newData: BitList) = BooleanColumn(label, newData, valueFormatter)

    override fun rename(newLabel: String) = BooleanColumn(newLabel, data, valueFormatter)

    override fun reorder(indices: IntArray): Column<BitList> {
        require(size == indices.size)
        return wrap(BitList(size) { data[indices[it]] })
    }

    override fun resize(newSize: Int): Column<BitList> {
        val copy = data.copy()
        if (copy.size() > newSize) {
            copy.clear(newSize, copy.size())
        }

        return wrap(BitList.wrap(newSize, copy))
    }

    override fun plus(other: Column<*>) = wrap(data + (other as BooleanColumn).data)

    @Deprecated("")
    override fun merge(other: Column<*>): Column<BitList> {
        throw UnsupportedOperationException()
    }

    override fun filter(mask: BitSet): Column<BitList> {
        val numRows = mask.cardinality()

        val filtered = BitList(numRows)
        var offset = 0
        for (i in 0 until size) {
            if (mask.get(i)) {
                filtered.set(offset++, data.get(i))
            }
        }
        return wrap(filtered)
    }

    override fun sorted(reverse: Boolean): IntArray {
        val indices = IntArray(size)
        val firstValue = reverse

        var offset = 0
        // fill first values
        for (row in 0 until size) {
            if (data[row] == firstValue) {
                indices[offset++] = row
            }
        }
        // fill second values
        for (row in 0 until size) {
            if (data[row] != firstValue) {
                indices[offset++] = row
            }
        }

        return indices
    }

    @Deprecated("")
    override fun intersect(c: Column<*>): ObjIntPredicate<BitList> {
        throw UnsupportedOperationException()
    }

    override val size: Int get() = data.size()

    override fun dump(row: Int): String {
        val value = data[row]
        return valueFormatter?.let { it(value) } ?: (if (data[row]) "1" else "0")
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is BooleanColumn -> false
        else -> label == other.label && data == other.data
    }

    override fun hashCode() = Objects.hash(label, data)

    override fun toString() = MoreObjects.toStringHelper(this)
        .add("label", label)
        .add("data", data)
        .toString()
}
