package org.jetbrains.bio.dataframe

import java.util.*

/**
 * A row-by-row builder for [DataFrame].
 *
 * Note, that constructing the data frame via [DataFrame.with]
 * calls might be more efficient.
 *
 * @author Evgeny Kurbatsky
 */
class DataFrameBuilder(private val spec: DataFrameSpec) {
    private val accumulator
            = if (spec.synchronized) Collections.synchronizedList(ArrayList<Array<Any>>()) else ArrayList<Array<Any>>()

    fun add(vararg row: Any) {
        accumulator.add(arrayOf(*row))
    }

    fun build(): DataFrame {
        val rowsNumber = accumulator.size
        val columns = spec.columns.mapIndexed { col, column ->
            val newColumn = column.resize(rowsNumber)
            for (row in 0 until rowsNumber) {
                val value = accumulator[row][col]
                try {
                    check(newColumn.boxedType.isInstance(value))
                    newColumn.load(row, value.toString())
                } catch (e: Exception) {
                    throw IllegalArgumentException("Wrong type: $row-th arg (column: ${newColumn.label}) value " +
                            "$value is expected to be of type ${newColumn.typeName()}")
                }
            }
            newColumn
        }.toList()

        return DataFrame(rowsNumber, columns)
    }
}
