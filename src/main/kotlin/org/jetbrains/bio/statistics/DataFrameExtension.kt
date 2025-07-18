package org.jetbrains.bio.statistics

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.DoubleColumn
import org.jetbrains.bio.dataframe.FloatColumn
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array

// Support save states to file as [Float], see #1163
fun F64Array.toFloatArray(): FloatArray {
    return FloatArray(this.length) { this[it].toFloat() }
}

// Support loading states to file as [Float], see #1163
fun DataFrame.f64Array(column: String): F64Array =
    when {
        this[column] is DoubleColumn -> this.sliceAsDouble(column).asF64Array()
        this[column] is FloatColumn -> {
            val floats = this.sliceAsFloat(column)
            DoubleArray(floats.size) { floats[it].toDouble() }.asF64Array()
        }

        else -> throw IllegalStateException("Expected Doubles or Floats: $column")
    }
