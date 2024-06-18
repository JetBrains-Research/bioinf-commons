package org.jetbrains.bio.genome

import org.jetbrains.bio.util.FileSize
import org.jetbrains.bio.util.formatLongNumber
import org.jetbrains.bio.util.presentablePath
import java.net.URI
import java.nio.file.Path

data class TrackAboutMetricValue<T>(
    val type: TrackAboutColumnType<T>,
    val value: T
)

abstract class TrackAboutColumnType<T> {
    abstract val name: String

    open fun comparator(): Comparator<T>? = null
    abstract fun valueClass(): Class<*>?
    open fun render(value: Any?) = value?.toString() ?: "-"
    infix fun to(value: T): TrackAboutMetricValue<T> = TrackAboutMetricValue(this, value)

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is TrackAboutColumnType<*>) return false
        if (name != other.name) return false
        if (valueClass() != other.valueClass()) return false
        return true
    }

    override fun hashCode(): Int {
        return name.hashCode() * 31 + valueClass().hashCode()
    }

}

object TrackAboutColumnTypes {
    val CT_FILE = TrackAboutStringColumnType("File")
    val CT_FILE_SIZE = TrackAboutFileSizeColumnType("File size")
}

class TrackAboutStringColumnType(
    override val name: String
) : TrackAboutColumnType<String>() {
    override fun valueClass() = String::class.java

    override fun comparator() = Comparator<String> { o1, o2 ->
        o1.lowercase().compareTo(o2.lowercase())
    }

    infix fun to(path: Path?) = TrackAboutMetricValue(this, path?.toString() ?: "n/a")
    infix fun to(uri: URI) = TrackAboutMetricValue(this, uri.presentablePath())
}

class TrackAboutLongColumnType(
    override val name: String
) : TrackAboutColumnType<Long>() {
    override fun valueClass() = Long::class.java

    override fun render(value: Any?) = when (value) {
        is Long -> value.toLong().formatLongNumber()
        else -> super.render(value)
    }

    override fun comparator() = Comparator<Long> { o1, o2 ->
        o1.compareTo(o2)
    }

    infix fun to(value: Int) = TrackAboutMetricValue(this, value.toLong())
}

class TrackAboutIntegerColumnType(
    override val name: String
) : TrackAboutColumnType<Int>() {
    override fun valueClass() = Int::class.java

    override fun render(value: Any?) = when (value) {
        is Int -> value.toLong().formatLongNumber()
        else -> super.render(value)
    }

    override fun comparator() = Comparator<Int> { o1, o2 ->
        o1.compareTo(o2)
    }
}

open class TrackAboutDoubleColumnType(
    override val name: String
) : TrackAboutColumnType<Double>() {

    override fun render(value: Any?): String {
        return if (value != null)"%.3f".format(value as Double).trimEnd('0') else "-"
    }

    override fun valueClass() = Double::class.java

    override fun comparator() = Comparator<Double> { o1, o2 ->
        o1.compareTo(o2)
    }
}

open class TrackAboutBooleanColumnType(
    override val name: String
) : TrackAboutColumnType<Boolean>() {
    override fun valueClass() = Boolean::class.java

    override fun comparator() = Comparator<Boolean> { o1, o2 ->
        o1.compareTo(o2)
    }
}


class TrackAboutPercentageColumnType(
    name: String
) : TrackAboutDoubleColumnType(name) {

    override fun render(value: Any?) = when (value) {
        is Double -> String.format("%.2f%%", value)
        else -> super.render(value)
    }
}

class TrackAboutFileSizeColumnType(
    override val name: String
) : TrackAboutColumnType<FileSize>() {

    override fun comparator() = Comparator<FileSize> { o1, o2 ->
        o1.bytes.compareTo(o2.bytes)
    }

    override fun valueClass() = FileSize::class.java
}