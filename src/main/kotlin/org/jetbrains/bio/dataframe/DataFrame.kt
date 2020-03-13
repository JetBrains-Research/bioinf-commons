package org.jetbrains.bio.dataframe

import com.google.common.base.MoreObjects
import com.google.common.base.Preconditions.*
import com.google.common.collect.ImmutableList
import com.google.common.collect.ObjectArrays
import gnu.trove.map.hash.TObjectIntHashMap
import java.io.IOException
import java.nio.file.Path
import java.util.*
import java.util.stream.Stream

@Suppress("unused")
class DataFrame @JvmOverloads constructor(
        val rowsNumber: Int = 0,
        internal val columns: List<Column<*>> = emptyList()) {

    val labels = columns.map { it.label.intern() }.toTypedArray()

    val columnsNumber: Int get() = columns.size
    val iloc = Slicer()

    /**
     * Alters the size of each column to match [rowsNumber].
     *
     * @see Column.resize for precise semantics.
     */
    fun resize(rowsNumber: Int): DataFrame {
        return DataFrame(rowsNumber, columns.map { it.resize(rowsNumber) })
    }

    /**
     * Reorders each column based on the values in [on]
     *
     * @param reverse if `true`, descending order is used.
     */
    fun reorder(on: String, reverse: Boolean = false): DataFrame {
        val indices = this[on].sorted(reverse)
        return DataFrame(rowsNumber, columns.map { it.reorder(indices) })
    }

    fun test(pf: RowPredicateFactory, startRow: Int = 0, endRow: Int = rowsNumber): BitterSet {
        checkPositionIndexes(startRow, endRow, rowsNumber)
        val rowPredicate = pf(this)
        return BitterSet(endRow - startRow) { rowPredicate.test(startRow + it) }
    }

    fun filter(pf: RowPredicateFactory, startRow: Int = 0, endRow: Int = rowsNumber) : DataFrame {
        val mask = test(pf, startRow, endRow)

        val fullMask = if (startRow == 0 && endRow == rowsNumber) {
            mask
        } else {
            val newMask = BitterSet(rowsNumber)

            var ptr = mask.nextSetBit(0)
            while (ptr >= 0) {
                newMask.set(ptr + startRow)
                ptr = mask.nextSetBit(ptr + 1)
            }
            newMask
        }
        return filter(fullMask)
    }

    fun filter(mask: BitSet): DataFrame {
        return DataFrame(mask.cardinality(), columns.map { it.filter(mask) })
    }

    /**
     * Returns a data frame with a given subset of columns.
     */
    fun only(vararg labels: String) = DataFrame(rowsNumber, labels.map { get(it) })

    /**
     * Returns a data frame without a given subset of columns.
     */
    fun omit(label: String, vararg rest: String): DataFrame {
        val labels = ObjectArrays.concat(label, rest)
        return DataFrame(rowsNumber, columns.filter { it.label !in labels })
    }

    fun with(rowsNumber: Int, column: Column<*>): DataFrame {
        require(columns.isEmpty() || rowsNumber == this.rowsNumber) {
            "#columns = ${columns.size}, #rows = ${this.rowsNumber}, #new column rows = $rowsNumber"
        }

        val columnBuilder = ImmutableList.builder<Column<*>>()
        val idx = getLabelIndexUnsafe(column.label)
        if (idx >= 0) {
            val columns = ArrayList(columns)
            columns[idx] = column
            columnBuilder.addAll(columns)
        } else {
            columnBuilder.addAll(columns)
            columnBuilder.add(column)
        }

        return DataFrame(rowsNumber, columnBuilder.build())
    }

    fun with(label: String, data: ByteArray) = with(data.size, ByteColumn(label, data))

    fun with(label: String, data: ShortArray) = with(data.size, ShortColumn(label, data))

    fun with(label: String, data: IntArray) = with(data.size, IntColumn(label, data))

    fun with(label: String, data: LongArray) = with(data.size, LongColumn(label, data))

    fun with(label: String, data: FloatArray) = with(data.size, FloatColumn(label, data))

    fun with(label: String, data: DoubleArray) = with(data.size, DoubleColumn(label, data))

    // TODO: could be optimized for factors.
    fun with(label: String, data: Array<String>) = with(data.size, StringColumn(label, data))

    fun with(label: String, data: BitterSet): DataFrame {
        return with(data.size(), BooleanColumn(label, data))
    }

    fun <T : Enum<T>> with(label: String, enumType: Class<T>, data: Array<T>): DataFrame {
        return with(data.size, EnumColumn(label, enumType, data))
    }

    fun getAsByte(r: Int, label: String) = sliceAsByte(label)[r]

    fun getAsShort(r: Int, label: String) = sliceAsShort(label)[r]

    fun getAsInt(r: Int, label: String) = sliceAsInt(label)[r]

    fun getAsFloat(r: Int, label: String) = sliceAsFloat(label)[r]

    fun getAsDouble(r: Int, label: String) = sliceAsDouble(label)[r]

    fun <T> getAsObj(r: Int, label: String): T = sliceAsObj<T>(label)[r]

    fun getAsBool(r: Int, label: String) = sliceAsBool(label)[r]

    fun sliceAsByte(label: String) = columns[getLabelIndex(label)].data as ByteArray

    fun sliceAsShort(label: String) = columns[getLabelIndex(label)].data as ShortArray

    fun sliceAsInt(label: String) = columns[getLabelIndex(label)].data as IntArray

    fun sliceAsLong(label: String) = columns[getLabelIndex(label)].data as LongArray

    fun sliceAsFloat(label: String) = columns[getLabelIndex(label)].data as FloatArray

    fun sliceAsDouble(label: String) = columns[getLabelIndex(label)].data as DoubleArray

    @Suppress("unchecked_cast")
    fun <T> sliceAsObj(label: String): Array<T> = columns[getLabelIndex(label)].data as Array<T>

    fun sliceAsBool(label: String) = columns[getLabelIndex(label)].data as BitterSet

    // XXX remove me if you can.
    fun rowAsDouble(r: Int): DoubleArray {
        checkElementIndex(r, rowsNumber)
        val result = DoubleArray(columnsNumber)
        for (c in 0 until columnsNumber) {
            result[c] = columns[c].getAsDouble(r)
        }
        return result
    }

    operator fun get(label: String) = columns[getLabelIndex(label)]

    private fun getLabelIndexUnsafe(label: String): Int {
        val strings = labels
        val size = strings.size
        for (i in 0 until size) {
            if (strings[i] === label) {
                return i
            }
        }
        return -1
    }

    private fun getLabelIndex(label: String): Int {
        val idx = getLabelIndexUnsafe(label)
        require(idx >= 0) {
            listOf("Unknown label '$label'.",
                    "Make sure that you are using interned string as a label!",
                    "Reference equality is used for lookup.",
                    "Known labels ${labels.contentToString()}").joinToString(" ")
        }

        return idx
    }

    @Throws(IOException::class)
    fun save(path: Path) = DataFrameMappers.forPath(path).save(path, this)

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is DataFrame -> false
        else -> rowsNumber == other.rowsNumber && columns == other.columns
    }

    override fun hashCode() = Objects.hash(rowsNumber, columns)

    override fun toString() = MoreObjects.toStringHelper(this)
            .add("rowsNumber", rowsNumber)
            .add("columns", '[' + labels.joinToString(", ") + ']')
            .toString()

    inner class Slicer {
        operator fun get(rowsRange: IntProgression) = this@DataFrame.apply {
            val mask = BitterSet(rowsNumber)
            if (rowsRange is IntRange) {
                val startRow = rowsRange.first
                val endRow = rowsRange.last + 1
                checkPositionIndexes(startRow, endRow, rowsNumber)
                mask.set(startRow, endRow)
            } else {
                for (i in rowsRange) {
                    checkPositionIndex(i, rowsNumber)
                    mask.set(i)
                }
            }
            return filter(mask)
        }
    }

    companion object {
        @Throws(IOException::class)
        @JvmStatic fun load(path: Path) = DataFrameMappers.forPath(path).load(path)

        /**
         * Performs an inner join of a list of data frames.
         * (~intersect by all common 'on' values)
         *
         * @param on column to join on, should be present in all data
         *           frames; duplicate values aren't allowed.
         * @param dfs data frames.
         * @return new data frame with join result sorted wrt to the
         *         join column.
         */
        @Suppress("unchecked_cast")
        @JvmStatic fun mergeInner(on: String, vararg dfs: DataFrame): DataFrame {
            require(dfs.size >= 2) { "expected at least two data frames" }

            var predicate: ObjIntPredicate<*> = ObjIntPredicate<Any> { _, _ -> true }
            for (i in 1 until dfs.size) {
                val l = dfs[i - 1][on]
                val r = dfs[i][on]
                // The '#and' call should be safe, because if 'l' is of different
                // type than 'r' we'll get an exception during the intersection.
                val copy = predicate as ObjIntPredicate<Any>
                predicate = (l.intersect(r) as ObjIntPredicate<Any>) and
                        ObjIntPredicate { data, j -> copy.test(data, j) }
            }

            val combinedPredicate = predicate
            val filtered = dfs.map { df ->
                val mask = df[on].test(combinedPredicate)
                df.filter(mask).reorder(on)
            }.toTypedArray()

            val first = filtered.first()
            return columnBind(setOf(on), *filtered)
                    .with(first.rowsNumber, first[on])
        }


        /**
         * Performs an outer join of a list of data frames
         * (~union by all possible 'on' values and fill undefined cells with default values)
         *
         * @param on column to join on, should be present in all data
         *           frames; duplicate values aren't allowed.
         * @param dfs data frames.
         * @return new data frame with join result sorted wrt to the
         *         join column.
         */
        @JvmStatic fun mergeOuter(on: String, vararg dfs: DataFrame): DataFrame {
            require(dfs.size >= 2) { "expected at least two data frames" }

            val combinedColumn = dfs.asSequence().map { it[on] }
                    .reduce(Column<*>::merge)
                    .let { it.reorder(it.sorted()) }

            val rowsNumber = combinedColumn.size
            val resized = dfs.map { df ->
                df.omit(on).resize(rowsNumber)
                        .with(rowsNumber, df[on].merge(combinedColumn))
                        .reorder(on)
            }.toTypedArray()

            return columnBind(setOf(on), *resized)
                    .with(rowsNumber, combinedColumn)
        }

        /**
         * Combines a list of data frames.
         *
         * The columns with same labels are suffixed with a number, e.g.
         * the column `"n"` will be renamed to `"n1"`. Exceptions are
         * excluded columns.
         *
         * @param exclude a set of column labels to drop from the combined
         *                data frame.
         * @param dfs a list of data frames to combine, each having the
         *            same number of rows and (possibly) different number
         *            of columns.
         * @return a new data frame.
         */
        @JvmStatic fun columnBind(exclude: Set<String>,
                                  vararg dfs: DataFrame): DataFrame {
            require(dfs.isNotEmpty()) { "no data" }

            val summary = Stream.of(*dfs).mapToInt { it.rowsNumber }
                    .summaryStatistics()
            val rowsNumber = summary.max
            require(rowsNumber == summary.min) { "different number of rows" }

            val common = TObjectIntHashMap<String>()
            for (label in dfs.flatMap { it.labels.asList() }) {
                common.adjustOrPutValue(label, 1, 1)
            }

            val columns = ArrayList<Column<*>>()
            val counts = TObjectIntHashMap<String>(common.size())
            for (column in dfs.flatMap { it.columns }) {
                val label = column.label
                if (label in exclude) {
                    continue
                }

                if (common[label] > 1) {
                    val suffix = counts.adjustOrPutValue(label, 1, 1).toString()
                    columns.add(column.rename(label + suffix))
                } else {
                    columns.add(column)
                }
            }

            return DataFrame(rowsNumber, columns)
        }

        @JvmStatic fun columnBind(vararg dfs: DataFrame): DataFrame {
            return columnBind(emptySet(), *dfs)
        }

        /**
         * Combines data frames with the same set of columns by rows.
         */
        @JvmStatic fun rowBind(df1: DataFrame, df2: DataFrame): DataFrame {
            val labels = df1.labels
            if (!labels.contentEquals(df2.labels)) {
                val chunks = arrayOf("columns do not match: ",
                        labels.contentToString(), " ",
                        df2.labels.contentToString())
                throw IllegalArgumentException(chunks.joinToString("\n"))
            }

            return when {
                df1.rowsNumber == 0 -> df2
                df2.rowsNumber == 0 -> df1
                else -> {
                    val columns = labels.map { df1[it] + df2[it] }
                    DataFrame(df1.rowsNumber + df2.rowsNumber, columns)
                }
            }
        }

        @JvmStatic fun rowBind(dfs: Array<DataFrame>): DataFrame {
            require(dfs.isNotEmpty()) { "expected at least one dataframe" }
            return dfs.reduce { a, b -> rowBind(a, b) }
        }
    }
}
