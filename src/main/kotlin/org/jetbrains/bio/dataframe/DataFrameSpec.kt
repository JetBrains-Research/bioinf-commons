package org.jetbrains.bio.dataframe

import com.google.common.collect.ObjectArrays
import java.nio.file.Path

/**
 * A spec holds data frame column names and types.
 *
 * @author Sergei Lebedev
 */
data class DataFrameSpec(
    val columns: MutableList<Column<*>> = ArrayList(),
    val synchronized: Boolean = false
) {

    fun ints(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(IntColumn(name, IntArray(0)))
        }
        return this
    }

    fun bytes(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(ByteColumn(name, ByteArray(0)))
        }
        return this
    }

    fun booleans(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(BooleanColumn(name, BitList(0)))
        }
        return this
    }

    fun longs(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(LongColumn(name, LongArray(0)))
        }
        return this
    }

    fun floats(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(FloatColumn(name, FloatArray(0)))
        }
        return this
    }

    fun doubles(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(DoubleColumn(name, DoubleArray(0)))
        }
        return this
    }

    fun <T : Enum<T>> enums(
        name: String,
        valueType: Class<T>
    ): DataFrameSpec {
        columns.add(
            EnumColumn(
                name, valueType, ObjectArrays.newArray(valueType, 0)
            )
        )
        return this
    }

    fun strings(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(StringColumn(name, emptyArray<String>()))
        }
        return this
    }

    fun shorts(vararg names: String): DataFrameSpec {
        for (name in names) {
            columns.add(ShortColumn(name, ShortArray(0)))
        }
        return this
    }

    fun builder() = DataFrameBuilder(this)

    companion object {
        /** A hack to make Kotlin accept the loaded type. */
        private enum class Traitor

        @Suppress("unchecked_cast")
        internal fun fromNamesAndTypes(
            names: List<String>,
            types: List<String>,
            source: Path
        ): DataFrameSpec {
            require(types.size == names.size) {
                "expected types for ${names.size} columns ($names), " +
                        "but got ${types.size}: $types"
            }

            val header = DataFrameSpec()
            for (i in types.indices) {
                val type = types[i]
                val name = names[i]

                when (type) {
                    "Byte" -> header.bytes(name)
                    "Integer" -> header.ints(name)
                    "Short" -> header.shorts(name)
                    "Long" -> header.longs(name)
                    "Float" -> header.floats(name)
                    "Double" -> header.doubles(name)
                    "Boolean" -> header.booleans(name)
                    "String" -> header.strings(name)
                    else -> {  // enum fqn given.
                        val valueType = try {
                            Class.forName(type) as Class<Traitor>
                        } catch (e: ClassNotFoundException) {
                            error("column type not found: $type in ${source.toAbsolutePath()}")
                        }

                        header.enums(name, valueType)
                    }
                }
            }

            return header
        }
    }
}
