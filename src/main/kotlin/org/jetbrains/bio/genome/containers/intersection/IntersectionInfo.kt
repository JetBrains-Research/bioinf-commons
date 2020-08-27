package org.jetbrains.bio.genome.containers.intersection

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode
import org.jetbrains.bio.dataframe.Column
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.DoubleColumn
import org.jetbrains.bio.dataframe.StringColumn
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path

/**
 * @author Roman.Chernyatchik
 * @date 2018-11-06
 */
data class IntersectionInfo(
        val data: F64Array,
        val mergedLociNumber: IntArray,
        val rowNames: List<String>,
        val colNames: List<String>) {


    fun saveToCSV(path: Path) {
        val (nrows, ncols) = data.shape

        CSVFormat.DEFAULT
                .withQuoteMode(QuoteMode.MINIMAL)
                .withCommentMarker('#')!!
                .print(path.bufferedWriter()).use { p ->

                    // header
                    p.print("")
                    p.print("#")
                    p.printRecord(colNames)

                    // data
                    for (i in 0 until nrows) {
                        p.print(rowNames[i])
                        p.print(mergedLociNumber[i])
                        for (cIdx in 0 until ncols) {
                            p.print(data[i, cIdx])
                        }
                        p.println()
                    }
                }
    }

    companion object {
        /**
         * Supports same format, as our WASHU tables
         */
        @Suppress("unused")
        fun loadFromCSV(path: Path) = CSVFormat.DEFAULT
                .withQuoteMode(QuoteMode.MINIMAL)
                .withCommentMarker('#')!!
                .parse(path.bufferedReader()).use {
                    var colNames: List<String>? = null

                    val data = ArrayList<DoubleArray>()
                    val rowNames = ArrayList<String>()
                    it.forEachIndexed { i, csvRecord ->
                        if (i == 0) {
                            colNames = csvRecord!!.toList().subList(1, csvRecord.size())
                        } else {
                            val records = csvRecord.toList()
                            val values = DoubleArray(colNames!!.size)
                            records.forEachIndexed { i, value ->
                                if (i == 0) {
                                    rowNames.add(value)
                                } else {
                                    values[i - 1] = value.toDouble()
                                }
                            }
                            data.add(values)
                        }
                    }

                    val cols = ArrayList<Column<*>>()
                    cols.add(StringColumn("rows", rowNames.toTypedArray()))
                    rowNames.zip(data).forEach { (colName, values) ->
                        cols.add(DoubleColumn(colName, values))
                    }

                    val mergedLociNumber: IntArray
                    val df = DataFrame(rowNames.size, cols).let { df ->
                        mergedLociNumber = when {
                            "#" in df.labels -> df.sliceAsInt("#")
                            else -> IntArray(rowNames.size) { 100 }    // legacy format w/o '#' column
                        }
                        df.omit("#")

                    }

                    val loadedRowNames = df.sliceAsObj<String>("rows").toList()
                    val loadedColNames = df.labels.toList().subList(1, df.columnsNumber)
                    val dataRowsN = loadedRowNames.size
                    val dataColsN = loadedColNames.size

                    val dataColumns = loadedColNames.map { col ->
                        df.sliceAsDouble(col)
                    }
                    val dataTable = F64Array(dataRowsN, dataColsN) { i, j ->
                        dataColumns[i][j]
                    }
                    IntersectionInfo(dataTable, mergedLociNumber, loadedRowNames, loadedColNames)
                }
    }


    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IntersectionInfo) return false

        if (data != other.data) return false
        if (!mergedLociNumber.contentEquals(other.mergedLociNumber)) return false
        if (rowNames != other.rowNames) return false
        if (colNames != other.colNames) return false

        return true
    }

    override fun hashCode(): Int {
        var result = data.hashCode()
        result = 31 * result + mergedLociNumber.contentHashCode()
        result = 31 * result + rowNames.hashCode()
        result = 31 * result + colNames.hashCode()
        return result
    }
}