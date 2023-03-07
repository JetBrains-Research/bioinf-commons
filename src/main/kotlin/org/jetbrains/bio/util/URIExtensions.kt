package org.jetbrains.bio.util

import com.google.common.io.Files
import com.google.common.io.Resources
import htsjdk.samtools.seekablestream.SeekableStreamFactory
import htsjdk.tribble.util.RemoteURLHelper
import org.slf4j.LoggerFactory
import java.io.IOException
import java.io.InputStreamReader
import java.net.HttpURLConnection
import java.net.HttpURLConnection.*
import java.net.URI
import java.net.URLDecoder
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

/**
 * @author Roman.Chernyatchik
 */
private val LOG = LoggerFactory.getLogger("org.jetbrains.bio.util.URIExtensions")

fun URI.isFile() = this.scheme == null || this.scheme == "file"

private val SUPPORTED_URL_PROTOCOLS = listOf("http", "https", "ftp")
fun URI.isSupportedURL() = this.scheme in SUPPORTED_URL_PROTOCOLS

fun URI.toPath(): Path {
    require(isFile()) { "Cannot convert URL to path: $this" }
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
    val lc = lowercase()
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
        val urlStr = this.toString()
        check(this.isSupportedURL()) { "URL not supported: $urlStr" }

        // If ftp url isn't accessible socket timeout while reading takes long time
        // in order to do fast check do:
        val url = this.toURL()

        val (urlAccessible, details) = if (url.protocol.lowercase().startsWith("http")) {
            // http or https
            var conn: HttpURLConnection? = null
            try {
                conn = url.openConnection() as HttpURLConnection
                conn.requestMethod = "HEAD"

                val status = conn.responseCode
                when (status) {
                    HTTP_MOVED_TEMP, HTTP_MOVED_PERM, HTTP_SEE_OTHER -> {
                        // most likely HTTP -> HTTPS redirect, don't do it for
                        // security reasons (passwords / cookies insecure transfer)
                        val redirectUrl = conn.getHeaderField("Location")
                        false to (" (http response $status, use $redirectUrl instead)")
                    }

                    HttpURLConnection.HTTP_OK -> true to ""
                    else -> false to " (http response $status)"
                }
            } catch (e: IOException) {
                LOG.debug("URL not accessible: $urlStr (${e.message})", e)
                false to " (${e.message})"
            } finally {
                conn?.disconnect()
            }
        } else {
            RemoteURLHelper(url).exists() to ""
        }
        if (!urlAccessible) {
            error("URL is not accessible: $urlStr$details")
        }

        try {
            // try to read one byre from url to be 100% sure:
            SeekableStreamFactory.getInstance().getStreamFor(this.toString()).use { it.read() }
        } catch (e: Exception) {
            throw IllegalStateException("Cannot read URL $urlStr", e)
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

fun URI.extension(): String = path.substringAfterLast(".")

fun URI.hasExt(vararg extensions: String): Boolean {
    val suffixes = extensions.map { ".${it.lowercase()}" }

    // try path part, e.g.:
    //    file:///mnt/stripe/foo.bed
    //    https://github.com/PetrTsurinov/BigBedTest/blob/master/ENCFF575VMI.bigBed?raw=true
    //    https://github.com/PetrTsurinov/BigBedTest/blob/master/ENCFF575VMI.bigBed
    val path = path.lowercase().replace("%2e", ".")
    val pathHasExt = suffixes.any { suffix -> path.endsWith(suffix) }
    if (pathHasExt || isFile()) {
        return pathHasExt
    }

    // try whole query, e.g:
    //   https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63874&format=file&file=GSE63874%5FNa%5FK2%5Fall%5Fminus%5Fcov%2Etdf
    // Just decode '.' or use full solution: URLDecoder.decode(urlValue, "UTF-8")
    val url = toString().lowercase().replace("%2e", ".")
    return suffixes.any { suffix -> url.endsWith(suffix) }
}

fun URI.stream() = when {
    isFile() -> toPath().inputStream()
    else -> asByteSource().openBufferedStream().streamFor(this.path)
}

fun URI.reader() = when {
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
