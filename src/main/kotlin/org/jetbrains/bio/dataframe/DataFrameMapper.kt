package org.jetbrains.bio.dataframe

import com.google.common.primitives.Primitives
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode
import org.jetbrains.bio.npy.NpzFile
import org.jetbrains.bio.util.extension
import java.io.IOException
import java.nio.file.Path

/**
 * A mapper implements data frame loading and saving logic.
 *
 * In most of the cases you should be good with [DataFrame.load] and
 * [DataFrame.save]. An appropriate mapper would be detected from file
 * extension.
 */
abstract class DataFrameMapper {
    /**
     * Tries to guess data frame spec for a given [path].
     */
    abstract fun guess(path: Path): DataFrameSpec

    /** Loads a data frame from a given [path]. */
    @Throws(IOException::class)
    fun load(path: Path): DataFrame {
        val spec = guess(path)
        return load(path, spec)
    }

    /**
     * Loads a data frame from a given [path] using the [spec] as a guide.
     */
    @Throws(IOException::class)
    abstract fun load(path: Path, spec: DataFrameSpec): DataFrame

    /**
     * Saves a data frame to a given [path].
     */
    @Throws(IOException::class)
    abstract fun save(path: Path, df: DataFrame)

}

fun DataFrame.dumpHead(rowsCount: Int): String {
    require(rowsCount <= rowsNumber) {
        "Cannot dump $rowsCount from df size $rowsNumber"
    }

    val buff = StringBuilder()
    DataFrameMappers.TSV.save(buff, resize(rowsCount), true, true);
    return buff.toString()
}

object DataFrameMappers {
    val TSV = CSVLike(
        CSVFormat.TDF.withQuoteMode(QuoteMode.MINIMAL).withCommentMarker('#')
            .withRecordSeparator("\n")!!
    )

    val CSV = CSVLike(
        CSVFormat.DEFAULT.withQuoteMode(QuoteMode.MINIMAL).withCommentMarker('#')
            .withRecordSeparator("\n")!!
    )

    val NPZ = object : DataFrameMapper() {
        override fun guess(path: Path) = NpzFile.read(path).use { reader ->
            val meta = reader.introspect()
            val names = meta.map { it.name }
            val types = meta.map { Primitives.wrap(it.type).simpleName }
            DataFrameSpec.fromNamesAndTypes(names, types, path)
        }

        override fun load(path: Path, spec: DataFrameSpec): DataFrame {
            return NpzFile.read(path).use { reader ->
                val columns = spec.columns.map { column ->
                    val values = reader[column.label]
                    when (column) {
                        is ByteColumn -> column.wrap(values.asByteArray())
                        is ShortColumn -> column.wrap(values.asShortArray())
                        is IntColumn -> column.wrap(values.asIntArray())
                        is LongColumn -> column.wrap(values.asLongArray())
                        is FloatColumn -> column.wrap(values.asFloatArray())
                        is DoubleColumn -> column.wrap(values.asDoubleArray())
                        is BooleanColumn -> {
                            val witness = values.asBooleanArray()
                            column.wrap(BitterSet(witness.size) { witness[it] })
                        }
                        is StringColumn -> column.wrap(values.asStringArray())
                        else -> error("unsupported column type: ${column.javaClass.canonicalName}")
                    }
                }

                val rowsNumber = columns.map { it.size }.distinct().single()
                DataFrame(rowsNumber, columns)
            }
        }

        override fun save(path: Path, df: DataFrame) {
            NpzFile.write(path).use { writer ->
                for (column in df.columns) {
                    when (column) {
                        is ByteColumn -> writer.write(column.label, column.data)
                        is ShortColumn -> writer.write(column.label, column.data)
                        is IntColumn -> writer.write(column.label, column.data)
                        is LongColumn -> writer.write(column.label, column.data)
                        is FloatColumn -> writer.write(column.label, column.data)
                        is DoubleColumn -> writer.write(column.label, column.data)
                        is BooleanColumn -> writer.write(
                            column.label, with(column.data) { BooleanArray(size()) { get(it) } })
                        is StringColumn -> writer.write(column.label, column.data)
                        else -> error("unsupported column type: ${column.javaClass.canonicalName}")
                    }
                }
            }
        }
    }

    /** Determines the appropriate mapper from file extension. */
    fun forPath(path: Path) = when (path.extension.lowercase()) {
        "npz" -> NPZ
        "csv" -> CSV
        else -> TSV
    }
}
