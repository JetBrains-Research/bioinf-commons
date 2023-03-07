package org.jetbrains.bio.util

import com.google.common.collect.ImmutableList
import com.google.common.hash.Hashing
import com.google.common.primitives.Longs
import kotlinx.support.jdk7.use
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.io.BufferedWriter
import java.io.File
import java.io.InputStream
import java.io.OutputStream
import java.net.URI
import java.nio.file.*
import java.nio.file.StandardCopyOption.ATOMIC_MOVE
import java.nio.file.StandardCopyOption.REPLACE_EXISTING
import java.nio.file.attribute.BasicFileAttributes
import java.text.DecimalFormat
import java.text.DecimalFormatSymbols
import java.util.zip.*
import kotlin.math.ln
import kotlin.math.pow

private val LOG = LoggerFactory.getLogger("org.jetbrains.bio.util.PathExtensions")

fun Array<String>.toPath(): Path {
    check(isNotEmpty())
    return Paths.get(this[0], *copyOfRange(1, size))
}

const val HASH_PREFIX = '#'

/** Creates # 5 digit-letter hash for given string */
val String.sha: String get() = HASH_PREFIX + Hashing.sha1().hashString(this, Charsets.UTF_8).toString().substring(0, 5)

fun String.toPath(): Path =
    // Windows case URI.getPath() returns path in such a format /C:/path/to/file
    if (isWindows() && this.startsWith(":/", 2)) {
        Paths.get(this.substring(1, this.length))
    } else
        Paths.get(this)

fun Path.resolve(vararg chunks: String): Path {
    return arrayOf(toString(), *chunks).toPath()
}

operator fun Path.div(other: String): Path = div(other.toPath())
operator fun Path.div(other: Path): Path = resolve(other)
operator fun String.div(other: String): Path = div(other.toPath())
operator fun String.div(other: Path): Path = toPath() / other

/** Check exists and readable. */
fun Path.checkAccessible() {
    if (notExists) {
        error("Track file doesn't exist: $this")
    } else if (!isReadable) {
        error("Cannot read file: $this, size $size")
    }
}

fun Path.isAccessible() = exists && isReadable

/** Returns `false` if a file exists and `true` otherwise. */
val Path.notExists: Boolean get() = Files.notExists(this)

/** Returns `true` if a file exists and `false` otherwise. */
val Path.exists: Boolean get() = Files.exists(this)

/** Tests whether a file can be read. */
val Path.isReadable: Boolean get() = Files.isReadable(this)

/** Tests whether a path points to a directory. */
val Path.isDirectory: Boolean get() = Files.isDirectory(this)

/** Tests whether a path is a regular file. */
val Path.isRegularFile: Boolean get() = Files.isRegularFile(this)

/** Returns the name of the file or directory denoted by this path. */
val Path.name: String get() = fileName.toString()

/** Returns the name of the file or directory without extension. */
val Path.stem: String
    get() {
        return name.substringBeforeLast(".$extension")
    }

/** Returns the name of the file or directory without extension.gz. */
val Path.stemGz: String
    get() {
        return when {
            Regex(".*\\.[a-z0-9_]+\\.gz", RegexOption.IGNORE_CASE).matches(name) ->
                name.substringBeforeLast('.').substringBeforeLast('.')

            else -> stem
        }
    }

/**
 * Returns the extension of this path (not including the dot), or
 * an empty string if it doesn't have one.
 */
val Path.extension: String get() = toFile().extension


/** Prepare string to form correct path. */
fun String.escape() = replace("[^a-zA-Z0-9_\\.]".toRegex(), "_")
    .replace("_{2,}".toRegex(), "_")
    .replace("^_|_$".toRegex(), "")
    .replace("_\\.".toRegex(), ".")
    .lowercase()

/** Creates a 5 digit-letter hash string for file full path */
val Path.sha: String get() = this.toAbsolutePath().toString().sha

/** File size pretty-printer. */
data class FileSize(
    /** Size in bytes as returned by [Files.size]. */
    val bytes: Long
) : Comparable<FileSize> {

    init {
        require(bytes == -1L || bytes >= 0) { "size must be >=0 or -1, but was: $bytes" }
    }

    override fun compareTo(other: FileSize) = Longs.compare(bytes, other.bytes)

    fun isAccessible() = bytes != -1L

    fun isEmpty() = bytes == 0L
    fun isNotEmpty() = !isEmpty()

    override fun toString(): String {
        if (bytes == -1L) {
            return "<not accessible>"
        }

        if (bytes == 0L) {
            return "0 b"
        }

        val units = (ln(bytes.toDouble()) / ln(1024.0)).toInt()
        val value = bytes / 1024.0.pow(units.toDouble())
        return "${FORMAT.format(value)} ${UNITS[units]}"
    }

    companion object {
        private val FORMAT = DecimalFormat("#,##0.#").apply {
            isGroupingUsed = false
            decimalFormatSymbols = DecimalFormatSymbols.getInstance().apply {
                decimalSeparator = ','
            }
        }

        /** Add terabytes when we get enough disk space. */
        private val UNITS = arrayOf("b", "kb", "mb", "gb")

        const val KB = 1024L
        const val MB = KB * 1024
        const val GB = MB * 1024
    }
}

/**
 * Returns file size suitable for pretty-printing.
 *
 * Use [FileSize.bytes] to get the result of [Files.size].
 */
val Path.size: FileSize get() = FileSize(if (isAccessible()) Files.size(this) else -1)

/**
 * Returns file size suitable for pretty-printing.
 *
 * Use [FileSize.bytes] to get the result of [Files.size].
 */
val URI.size: FileSize get() = FileSize(if (isAccessible()) length() else -1)

/** Returns a new path with the [name] changed. */
fun Path.withName(newName: String): Path {
    check(name.isNotEmpty())
    return parent / newName
}

/** Returns a new path with the [extension] changed. */
fun Path.withExtension(newExtension: String): Path {
    require(newExtension.isNotEmpty())
    return parent / "$stem.$newExtension"
}

/** Returns a new path with the [stem] changed. */
fun Path.withStem(newStem: String): Path {
    return if (extension.isEmpty()) {
        parent / newStem
    } else {
        parent / "$newStem.$extension"
    }
}

fun Path.walkFileTree(walker: FileVisitor<Path>) {
    Files.walkFileTree(this, walker)
}

fun Path.list(): List<Path> {
    return Files.list(this).use { s ->
        return@use ImmutableList.copyOf(s.iterator())
    }
}

fun String.toPathMatcher(): PathMatcher {
    // We enforce glob syntax, because it's human-readable.
    return FileSystems.getDefault().getPathMatcher("glob:$this")
}

/**
 * Recursively glob a directory.
 * Note that the current implementation only supports glob for *files*.
 */
fun Path.glob(pattern: String): List<Path> {
    check(exists) { "$this doesn't exist" }
    check(isDirectory) { "$this is not a directory" }

    val matcher = "${toAbsolutePath()}${File.separatorChar}$pattern"
        .let { if (isWindows()) it.replace("\\", "\\\\") else it }
        .toPathMatcher()
    val matched = ArrayList<Path>()
    walkFileTree(object : SimpleFileVisitor<Path>() {
        override fun visitFile(file: Path, attrs: BasicFileAttributes): FileVisitResult {
            if (matcher.matches(file)) {
                matched.add(file)
            }

            return FileVisitResult.CONTINUE
        }
    })

    return matched
}

fun isWindows() = System.getProperty("os.name")?.lowercase()?.contains("windows") ?: false

fun Path.createDirectories(): Path {
    if (notExists) {
        Files.createDirectories(this)
    }
    return this
}

fun Path.delete() = Files.delete(this)

fun Path.deleteIfExists() {
    if (isDirectory) {
        if (exists) deleteDirectory()
        return
    }

    Files.deleteIfExists(this)
}

fun Path.deleteDirectory() {
    toFile().deleteRecursively()
}

fun InputStream.copy(path: Path, vararg options: StandardCopyOption) {
    Files.copy(this, path, *options)
}

fun Path.copy(outputStream: OutputStream) {
    Files.copy(this, outputStream)
}

fun Path.move(target: Path, vararg options: StandardCopyOption) {
    if (target.exists && target.isDirectory) {
        Files.move(this, target.resolve(this.fileName), *options)
    } else {
        Files.move(this, target, *options)
    }
}

fun Path.copy(target: Path, vararg options: StandardCopyOption) {
    if (target.exists && target.isDirectory) {
        Files.copy(this, target.resolve(this.fileName), *options)
    } else {
        Files.copy(this, target, *options)
    }
}

/**
 * Checks if [chunk] is in path
 */
fun Path.hasChunk(chunk: String): Boolean =
    chunk.lowercase() in fileName.toString().lowercase().split('_', '.')

fun <T> Path.useLines(consumer: (Sequence<String>) -> T): T = Files.newBufferedReader(this).useLines(consumer)

fun Path.read() = String(Files.readAllBytes(this))

fun Path.write(buf: String, vararg options: StandardOpenOption): Path {
    return write(buf.toByteArray(), *options)
}

fun Path.write(buf: ByteArray, vararg options: StandardOpenOption): Path {
    parent.createDirectories()
    return Files.write(this, buf, *options)
}

/**
 * Checks if file exists and launches computation in case of missing file, thread-safe.
 *
 * @param label a human-readable description of the computation.
 * @param isDirectory If true path to recalculate is a directory, instead of file
 * @param ignoreEmptyFile   Do not throw exception on produced empty file
 * @param timestamp Recalculate file if it exists, but older than timestamp
 * @param recalculate a function for populating the path. To ensure
 *                    atomicity the result of [recalculate] must be
 *                    first written to a temporary path and then moved
 *                    to this path.
 *

 */
@JvmOverloads
fun Path.checkOrRecalculate(
    label: String = "",
    isDirectory: Boolean = false,
    ignoreEmptyFile: Boolean = false,
    timestamp: Long = 0,
    recalculate: (PathWrapper) -> Unit
): Path {
    // IMPORTANT: synchronize on interned path string!
    val target = toAbsolutePath().normalize().toString().intern()
    synchronized(target) {
        val prefix = if (label.isNotEmpty()) "$label: " else ""
        val ts = if (exists) toFile().lastModified() else -1
        if (exists && (size.isNotEmpty() || ignoreEmptyFile || isDirectory) && (ts > timestamp)) {
            LOG.debug("${prefix}exists ($size): $target")
        } else {
            val targetType = if (isDirectory) "Directory" else "File"
            if (!exists) {
                LOG.debug("$prefix$targetType is missing: $target")
            } else if (ts <= timestamp) {
                LOG.debug("$prefix$targetType is outdated (file $ts <= $timestamp): $target")
            } else {
                LOG.debug("$prefix$targetType is zero-sized: $target")
            }

            val outdated = if (exists && toFile().lastModified() < timestamp) "outdated " else ""
            LOG.time(level = Level.INFO, message = "${prefix}processing $outdated$target") {
                parent.createDirectories()
                val block: (Path) -> Unit = { tmpPath ->
                    recalculate(PathWrapper(tmpPath))
                    if (ignoreEmptyFile || (isDirectory || tmpPath.size.isNotEmpty())) {
                        tmpPath.move(this, ATOMIC_MOVE, REPLACE_EXISTING)
                        LOG.debug("${prefix}processed $target ($size written)")
                    } else {
                        tmpPath.deleteIfExists()
                        throw IllegalStateException("${prefix}produced a zero-sized file: $tmpPath")
                    }
                }
                if (isDirectory) {
                    withTempDirectory("${stem}_tmp", dir = parent, block = block)
                } else {
                    withTempFile("${stem}_tmp", ".$extension", dir = parent, block = block)
                }
            }
        }
        return this@checkOrRecalculate
    }
}

/**
 * This is to make sure the caller of [checkOrRecalculate] didn't
 * ignore the output path, which is easy to do in Kotlin.
 */
data class PathWrapper(val path: Path)

/**
 * Creates temp file, passes it to a given closure and deletes after the closure has been invoked.
 *
 * If you open related temp file inside the closure and leave it opened, this method will fail to
 * delete temp file on Windows with error like:
 *      java.nio.file.FileSystemException: ...tmp\track305430928875073443.bw: The process cannot access the
 *      file because it is being used by another process.
 *
 * For some reason you'll get no error on linux, but error is quite expected and reasonable.
 */
inline fun <T> withTempFile(
    prefix: String, suffix: String, dir: Path? = null,
    block: (Path) -> T
): T {
    val tempFile = if (dir == null) {
        Files.createTempFile(prefix, suffix)
    } else {
        dir.createDirectories()
        Files.createTempFile(dir, prefix, suffix)
    }
    return try {
        block(tempFile)
    } finally {
        tempFile.deleteIfExists()
    }
}

inline fun <T> withTempDirectory(
    prefix: String,
    dir: Path? = null,
    block: (Path) -> T
): T {
    val tempDir = if (dir == null) {
        Files.createTempDirectory(prefix)
    } else {
        dir.createDirectories()
        Files.createTempDirectory(dir, prefix)
    }

    return try {
        block(tempDir)
    } finally {
        tempDir.deleteIfExists()
    }
}

/**
 * Returns a buffered input stream for this path.
 */
fun Path.inputStream(vararg options: OpenOption) =
    Files.newInputStream(this, *options).buffered().streamFor(name)

fun InputStream.streamFor(path: String) = path.let {
    val lcPath = it.lowercase()
    val parentStream = when {
        IOMonitor.debugMode && this !is IOMonitorInputStream -> IOMonitorInputStream(this)
        else -> this
    }
    when {
        lcPath.endsWith(".gz") -> GZIPInputStream(parentStream)
        lcPath.endsWith(".zip") ->
            // This only works for single-entry ZIP files.
            ZipInputStream(parentStream).apply { nextEntry }

        else -> parentStream
    }
}

fun Path.bufferedReader(vararg options: OpenOption) = inputStream(*options).bufferedReader()

/**
 * Returns a buffered output stream for this path.
 */
fun Path.outputStream(vararg options: OpenOption): OutputStream {
    val outputStream = Files.newOutputStream(this, *options).buffered()
    return when (extension.lowercase()) {
        "gz" -> GZIPOutputStream(outputStream)
        "zip" -> ZipOutputStream(outputStream).apply {
            // Without this line ZIP file will be corrupted.
            putNextEntry(ZipEntry("jbr_rulezzz.txt"))
        }

        else -> outputStream
    }
}

fun Path.bufferedWriter(vararg options: OpenOption): BufferedWriter {
    // Use Path#toAbsolutePath(), otherwise it may have no parent
    toAbsolutePath().parent.createDirectories()
    return outputStream(*options).bufferedWriter()
}

fun Path.touch() = apply {
    // touch is also supposed to update file modification date
    com.google.common.io.Files.touch(toFile())
}
