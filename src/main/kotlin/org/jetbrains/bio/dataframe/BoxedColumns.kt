package org.jetbrains.bio.dataframe

import com.google.common.base.MoreObjects
import com.google.common.collect.ComparisonChain
import com.google.common.collect.ObjectArrays
import java.util.*

@Suppress("unchecked_cast")
abstract class ObjColumn<T: Any> protected constructor(
        label: String, valueType: Class<T>, data: Array<T>)
:
        Column<Array<T>>(label, data, valueType) {

    override fun getAsDouble(row: Int): Double {
        check(false) { "not a numeric column" }
        return Double.NaN
    }

    override fun plus(other: Column<*>): Column<Array<T>> {
        return wrap(ObjectArrays.concat(data, other.data as Array<T>,
                                        boxedType as Class<T>))
    }

    @Deprecated("")
    override fun merge(other: Column<*>): Column<Array<T>> {
        throw UnsupportedOperationException()
    }

    override fun filter(mask: BitSet): Column<Array<T>> {
        val newData = ObjectArrays.newArray(boxedType as Class<T>, mask.cardinality())
        var offset = 0
        for ((i, value) in data.withIndex()) {
            if (mask[i]) {
                newData[offset++] = value
            }
        }

        return wrap(newData)
    }

    override fun sorted(reverse: Boolean): IntArray {
        throw UnsupportedOperationException()
    }

    override fun reorder(indices: IntArray): Column<Array<T>> {
        val clone = data.clone()
        for (i in 0..size - 1) {
            clone[i] = data[indices[i]]
        }

        return wrap(clone)
    }

    override fun resize(newSize: Int) = wrap(Arrays.copyOf(data, newSize))

    override fun intersect(c: Column<*>): ObjIntPredicate<Array<T>> {
        val seen = data.intersect((c as ObjColumn<T>).data.asList())
        return ObjIntPredicate { data, i -> data[i] in seen }
    }

    override val size: Int get() = data.size

    override fun typeName() = boxedType.name

    override fun dump(row: Int) = data[row].toString()

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is ObjColumn<*> -> false
        else -> label == other.label && Arrays.equals(data, other.data)
    }

    override fun hashCode() = Arrays.deepHashCode(arrayOf(label, data))

    override fun toString() = MoreObjects.toStringHelper(this)
            .add("label", label)
            .add("data", Arrays.toString(data))
            .toString()
}

class StringColumn(label: String, data: Array<String>) :
        ObjColumn<String>(label, String::class.java, data) {

    override fun rename(newLabel: String) = StringColumn(newLabel, data)

    override fun wrap(newData: Array<String>) = StringColumn(label, newData)

    @Suppress("unchecked_cast")
    override fun merge(other: Column<*>): Column<Array<String>> {
        val seen = data.toMutableSet()
        check(seen.size == data.size) { "duplicates found" }

        val merged = data.toMutableList()
        for (value in other.data as Array<String>) {
            if (value !in seen) {
                merged.add(value)
                seen.add(value)
            }
        }

        return wrap(merged.toTypedArray())
    }

    override fun sorted(reverse: Boolean): IntArray {
        val comparator = Comparator<IndexedValue<String>> { p1, p2 ->
            ComparisonChain.start()
                    .compare(p1.value, p2.value)
                    .compare(p1.index, p2.index)
                    .result()
        }

        return data.withIndex()
                .sortedWith(if (reverse) comparator.reversed() else comparator)
                .map { it.index }
                .toIntArray()
    }

    override fun load(row: Int, value: String) {
        data[row] = value
    }

    override fun typeName() = boxedType.simpleName
}

@Suppress("unchecked_cast")
class EnumColumn<T : Enum<T>>(label: String, val enumType: Class<T>,
                              data: Array<T>)
:
        ObjColumn<T>(label, enumType, data) {

    override fun rename(newLabel: String): Column<Array<T>> {
        return EnumColumn(newLabel, boxedType as Class<T>, data)
    }

    override fun wrap(newData: Array<T>): EnumColumn<T> {
        return EnumColumn(label, boxedType as Class<T>, newData)
    }

    override fun sorted(reverse: Boolean): IntArray {
        val ordinals = Arrays.stream(data).mapToInt { it.ordinal }.toArray()
        return ordinals.argSort(reverse)
    }

    override fun load(row: Int, value: String) {
        data[row] = java.lang.Enum.valueOf<T>(boxedType as Class<T>, value)
    }
}
