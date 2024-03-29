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

interface TrackAboutColumnType<T> {
    val name: String

    fun comparator(): Comparator<T>? = null
    fun valueClass(): Class<*>?
    fun render(value: Any?) = value?.toString() ?: "-"
    infix fun to(value: T): TrackAboutMetricValue<T> = TrackAboutMetricValue(this, value)
}

object TrackAboutColumnTypes {
    val CT_FILE = TrackAboutStringColumnType("File")
    val CT_FILE_SIZE = TrackAboutFileSizeColumnType("File size")
}

class TrackAboutStringColumnType(
    override val name: String
) : TrackAboutColumnType<String> {
    override fun valueClass() = String::class.java

    override fun comparator() = Comparator<String> { o1, o2 ->
        o1.lowercase().compareTo(o2.lowercase())
    }

    infix fun to(path: Path?) = TrackAboutMetricValue(this, path?.toString() ?: "n/a")
    infix fun to(uri: URI) = TrackAboutMetricValue(this, uri.presentablePath())
}

class TrackAboutLongColumnType(
    override val name: String
) : TrackAboutColumnType<Long> {
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
) : TrackAboutColumnType<Int> {
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
) : TrackAboutColumnType<Double> {
    override fun valueClass() = Double::class.java

    override fun comparator() = Comparator<Double> { o1, o2 ->
        o1.compareTo(o2)
    }
}

open class TrackAboutBooleanColumnType(
    override val name: String
) : TrackAboutColumnType<Boolean> {
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
) : TrackAboutColumnType<FileSize> {

    override fun comparator() = Comparator<FileSize> { o1, o2 ->
        o1.bytes.compareTo(o2.bytes)
    }

    override fun valueClass() = FileSize::class.java
}