package org.jetbrains.bio.genome.containers.intersection

import gnu.trove.list.array.TIntArrayList
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.viktor.F64Array
import java.io.StringWriter
import java.io.Writer
import java.nio.file.Path

/**
 * @author Roman.Chernyatchik
 * @date 2018-11-06
 */
data class IntersectionInfo(
    val data: F64Array,
    val rowsMergedLociNumber: IntArray,
    val colsMergedLociNumber: IntArray,
    val rowNames: List<String>,
    val colNames: List<String>
) {

    fun saveToCSV(writer: Writer) {
        val (nrows, ncols) = data.shape

        CSVFormat.DEFAULT
            .withQuoteMode(QuoteMode.MINIMAL)
            .withCommentMarker('#')!!
            .print(writer).use { p ->
                // header
                p.print("a_names")
                p.print("a_merged_loci_count")
                p.printRecord(colNames)

                // cols sizes
                p.print("b_merged_loci_count")
                p.print(".")
                for (i in 0 until ncols) {
                    p.print(colsMergedLociNumber[i])
                }
                p.println()

                // data
                for (i in 0 until nrows) {
                    p.print(rowNames[i])
                    p.print(rowsMergedLociNumber[i])
                    for (cIdx in 0 until ncols) {
                        p.print(data[i, cIdx])
                    }
                    p.println()
                }
            }
    }

    fun saveToCSV(path: Path) {
        saveToCSV(path.bufferedWriter())
    }

    override fun toString() = StringWriter()
        .also { saveToCSV(it) }
        .toString()
        .replace("\r", "")

    companion object {
        /**
         * Supports same format, as our WASHU tables
         */
        @Suppress("unused")
        fun loadFromCSV(path: Path) = CSVFormat.DEFAULT
            .withQuoteMode(QuoteMode.MINIMAL)
            .withCommentMarker('#')!!
            .parse(path.bufferedReader()).use { parser ->
                var colNames: List<String>? = null
                var colsMergedLociNumber: IntArray? = null
                val rowsMergedLociNumber = TIntArrayList()

                val data = ArrayList<DoubleArray>()
                val rowNames = ArrayList<String>()
                parser.forEachIndexed { recIdx, csvRecord ->
                    when (recIdx) {
                        0 -> {
                            colNames = csvRecord!!.toList().subList(2, csvRecord.size())
                        }
                        1 -> {
                            requireNotNull(colNames)
                            colsMergedLociNumber = IntArray(colNames!!.size) { i ->
                                csvRecord.get(i + 2).toInt()
                            }
                        }
                        else -> {
                            val records = csvRecord.toList()
                            val rowValues = DoubleArray(colNames!!.size)
                            records.forEachIndexed { i, value ->
                                when (i) {
                                    0 -> rowNames.add(value)
                                    1 -> rowsMergedLociNumber.add(value.toInt())
                                    else -> rowValues[i - 2] = (value.toDouble())
                                }
                            }
                            data.add(rowValues)
                        }
                    }
                }

                if (colNames == null) {
                    IntersectionInfo(
                        F64Array(0, 0),
                        IntArray(0), IntArray(0),
                        emptyList(), emptyList()
                    )
                } else {
                    requireNotNull(colNames)
                    requireNotNull(colsMergedLociNumber)

                    IntersectionInfo(
                        F64Array(rowNames.size, colNames!!.size) { i, j ->
                            data[i][j]
                        },
                        rowsMergedLociNumber.toArray(), colsMergedLociNumber!!,
                        rowNames, colNames!!
                    )
                }
            }
    }


    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IntersectionInfo) return false

        if (data != other.data) return false
        if (!rowsMergedLociNumber.contentEquals(other.rowsMergedLociNumber)) return false
        if (!colsMergedLociNumber.contentEquals(other.colsMergedLociNumber)) return false
        if (rowNames != other.rowNames) return false
        if (colNames != other.colNames) return false

        return true
    }

    override fun hashCode(): Int {
        var result = data.hashCode()
        result = 31 * result + rowsMergedLociNumber.contentHashCode()
        result = 31 * result + colsMergedLociNumber.contentHashCode()
        result = 31 * result + rowNames.hashCode()
        result = 31 * result + colNames.hashCode()
        return result
    }
}