package org.jetbrains.bio.dataframe

import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.size
import java.nio.file.Path

class CSVLike(val format: CSVFormat) : DataFrameMapper() {
    private val LOG = Logger.getLogger(CSVLike::class.java)

    companion object {
        const val DELIMITER = "; "
    }

    /**
     * Tries to guess data frame spec for a given [path].
     */
    override fun guess(path: Path): DataFrameSpec {
        val row = format.parse(path.bufferedReader()).use {
            checkNotNull(it.firstOrNull()) {
                "${path.toAbsolutePath()} is empty, no header given."
            }
        }


        val names = row.toList()
        val types = checkNotNull(row.comment) {
            "${path.toAbsolutePath()} doesn't contain typed header"
        }.trim().split("\\s*;\\s*".toRegex())
        // ^^^ the format is # [Type; ]+

        return DataFrameSpec.fromNamesAndTypes(names, types, path)
    }

    override fun load(path: Path, spec: DataFrameSpec) = load(path, spec, true)

    fun load(path: Path, spec: DataFrameSpec, header: Boolean): DataFrame {
        LOG.info("Loading data frame $path ${path.size}")
        val linesNumber = path.bufferedReader().use {
            it.lines().mapToInt { line -> if (line[0] != format.commentMarker) 1 else 0 }.sum()
        }
        val rowsNumber = if (linesNumber == 0) 0 else linesNumber - (if (header) 1 else 0)
        LOG.info("Columns: ${spec.columns.size} Rows: $rowsNumber")
        val progress = Progress { title = "Reading data frame $path" }.bounded(rowsNumber.toLong())
        val df = DataFrame(rowsNumber, spec.columns.map { it.resize(rowsNumber) })
        path.bufferedReader().use {
            val format = format.withHeader(*df.labels).withSkipHeaderRecord(header)
            for ((i, row) in format.parse(it).withIndex()) {
                // XXX we allow row to contain more columns because it's often
                // the case for UCSC annotations :(
                check(row.size() >= df.columnsNumber) { "inconsistent record $row" }
                for (col in 0 until df.columnsNumber) {
                    val column = df.columns[col]
                    val value = row[col]
                    try {
                        column.load(i, value)
                    } catch (e: Exception) {
                        throw IllegalArgumentException("Wrong type: $i-th arg (column: ${column.label}) value " +
                                "$value is expected to be of type ${column.typeName()}")
                    }
                }
                progress.report()
            }
        }
        progress.done()
        return df
    }

    override fun save(path: Path, df: DataFrame) = save(path.bufferedWriter(), df, true, true)

    fun save(out: Appendable,
             df: DataFrame,
             header: Boolean = true,
             typed: Boolean = true,
             comments: Array<String>? = null) {
        format.print(out).use { csvPrinter ->
            comments?.forEach { csvPrinter.printComment(it) }

            if (typed) {
                csvPrinter.printComment(df.columns.joinToString(DELIMITER) { it.typeName() })
            }

            if (header) {
                csvPrinter.printRecord(*df.labels)
            }

            for (r in 0 until df.rowsNumber) {
                for (column in df.columns) {
                    csvPrinter.print(column.dump(r))
                }
                csvPrinter.println()
            }
        }
    }
}