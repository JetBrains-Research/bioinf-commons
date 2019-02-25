package org.jetbrains.bio.util

import com.google.common.io.Files
import com.google.common.io.Resources
import htsjdk.samtools.seekablestream.SeekableStreamFactory
import htsjdk.tribble.util.RemoteURLHelper
import java.io.InputStreamReader
import java.net.URI
import java.net.URLDecoder
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

/**
 * @author Roman.Chernyatchik
 */

fun URI.isFile() = this.scheme == null || this.scheme == "file"

private val SUPPORTED_URL_PROTOCOLS = listOf("http", "https", "ftp")
fun URI.isSupportedURL() = this.scheme in SUPPORTED_URL_PROTOCOLS

fun URI.toPath(): Path {
    require(isFile()) { "Cannot convert URL to path: $this"}
    return (if (scheme == null) Paths.get(toString()) else Paths.get(this)).normalize()
}

fun URI.presentablePath() = when {
    isFile() -> toPath().toString()
    else -> URLDecoder.decode(toString(), "UTF-8")
}

fun URI.shortName() = when {
    isFile() -> toPath().name
    !path.isNullOrEmpty() && query == null -> {
        URLDecoder.decode(path.substringAfterLast('/'), "UTF-8")
    }
    else -> URLDecoder.decode(toString(), "UTF-8")
}

private val URL_REGEXP = Pattern.compile("^[a-z]{3,10}:/(/){0,2}[a-z0-9].*$")
fun String.toUri(): URI {
    // 'file:/' urls:
    val lc = toLowerCase()
    if (lc.startsWith("file:/")) {
        // e.g. user defined yaml config, let's replace spaces at least
        return URI.create(this.replace(" ", "%20"))
    }
    
    val isURL = URL_REGEXP.matcher(lc).matches()
    return if (!isURL) {
        this.toPath().toUri()
    } else {
        // throw error if ur is malformed
        URI.create(this)
    }
}

fun URI.checkAccessible() {
    if (isFile()) {
        toPath().checkAccessible()
    } else {
        val url = this.toString()
        check(this.isSupportedURL()) { "URL not supported: $url" }

        // If ftp url isn't accessible socket timeout while reading takes long time
        // in order to do fast check do:
        if (!RemoteURLHelper(this.toURL()).exists()) {
            error("URL is not accessible $url")
        }

        try {
            // try to read one byre from url to be 100% sure:
            SeekableStreamFactory.getInstance().getStreamFor(this.toString()).use { it.read() }
        } catch (e: Exception) {
            throw IllegalStateException("Cannot read URL $url", e)
        }
    }
}

fun URI.isAccessible() = if (isFile()) {
        toPath().isAccessible()
    } else {
        // If ftp url isn't accessible socket timeout while reading takes long time
        // in order to do fast check do:
        var accessible = isSupportedURL() && RemoteURLHelper(this.toURL()).exists()
        if (accessible) {
            try {
                // try to read one byte from url
                SeekableStreamFactory.getInstance().getStreamFor(this.toString()).use { it.read() }
            } catch (e: Exception) {
                accessible = false
            }
        }
        accessible
    }

fun URI.asByteSource() = if (isFile()) {
    Files.asByteSource(toPath().toFile())!!
} else {
    Resources.asByteSource(toURL())!!
}

fun URI.hasExt(vararg extensions: String): Boolean {
    val suffixes = extensions.map { ".${it.toLowerCase()}" }

    // try path part, e.g.:
    //    file:///mnt/stripe/foo.bed
    //    https://github.com/PetrTsurinov/BigBedTest/blob/master/ENCFF575VMI.bigBed?raw=true
    //    https://github.com/PetrTsurinov/BigBedTest/blob/master/ENCFF575VMI.bigBed
    val path = path.toLowerCase().replace("%2e", ".")
    val pathHasExt = suffixes.any { suffix -> path.endsWith(suffix) }
    if (pathHasExt || isFile()) {
        return pathHasExt
    }

    // try whole query, e.g:
    //   https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63874&format=file&file=GSE63874%5FNa%5FK2%5Fall%5Fminus%5Fcov%2Etdf
    // Just decode '.' or use full solution: URLDecoder.decode(urlValue, "UTF-8")
    val url = toString().toLowerCase().replace("%2e", ".")
    return suffixes.any { suffix -> url.endsWith(suffix) }
}

fun URI.stream() = when {
    isFile() -> toPath().inputStream()
    else -> asByteSource().openBufferedStream().streamFor(this.path)
}

fun URI.reader()  = when {
    isFile() -> toPath().bufferedReader()
    else -> InputStreamReader(asByteSource().openBufferedStream().streamFor(this.path), Charsets.UTF_8).buffered()
}

/**
 * Source length in bytes
 */
fun URI.length() = if (isFile()) {
    toPath().toFile().length()
} else {
    SeekableStreamFactory.getInstance().getStreamFor(toString()).length()
}
