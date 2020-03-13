package org.jetbrains.bio.genome

import org.apache.commons.net.ftp.FTPClient
import org.apache.http.HttpEntity
import org.apache.http.client.HttpClient
import org.apache.http.client.HttpResponseException
import org.apache.http.client.config.RequestConfig
import org.apache.http.client.methods.HttpGet
import org.apache.http.client.utils.HttpClientUtils
import org.apache.http.config.SocketConfig
import org.apache.http.conn.ConnectTimeoutException
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.impl.client.LaxRedirectStrategy
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.FileOutputStream
import java.io.IOException
import java.net.SocketTimeoutException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

private fun HttpClient.tryGet(url: String): HttpEntity? {
    val response = execute(HttpGet(url))
    val statusLine = response.statusLine
    if (statusLine.statusCode / 100 != 2) {
        throw HttpResponseException(statusLine.statusCode, "${statusLine.reasonPhrase}: $url")
    }

    return response.entity
}

/**
 * Downloads a URL to a given path.
 *
 * @param outputPath path to write to.
 * @param timeout HTTP connection timeout (in seconds).
 * @param maxRetries maximum number of download attempts.
 */
fun String.downloadTo(outputPath: Path,
                      timeout: Int = 30,
                      maxRetries: Int = 10) {
    outputPath.parent.createDirectories()

    try {
        UCSC.LOG.info("Download [${this}] to $outputPath")
        val timeoutMs = timeout * 1000
        if (this.startsWith("ftp://")) {
            downloadFtp(outputPath)
        } else {
            downloadHttp(timeoutMs, maxRetries, outputPath)
        }
    } catch (e: Exception) {
        throw RuntimeException("Failed to download $this to $outputPath", e)
    }
}

private fun String.downloadFtp(outputPath: Path) {
    val prefix = "ftp://anonymous@"
    check(this.startsWith(prefix)) {
        "Only anonymous FTP access is supported, but url is: ${this}"
    }
    val user = "anonymous"
    val pass = "anonymous"

    val hostAndPath = this.substring(prefix.length)
    val idx = hostAndPath.indexOf('/')
    check(idx > 0) {
        "Unsupported url: ${this}"
    }
    val ftpHost = hostAndPath.substring(0, idx)
    val ftpPath = hostAndPath.substring(idx)

    val ftpClient = FTPClient()

    ftpClient.connect(ftpHost, 21)
    ftpClient.login(user, pass)
    ftpClient.enterLocalPassiveMode()
    ftpClient.setFileType(FTPClient.BINARY_FILE_TYPE)

    val output = FileOutputStream(outputPath.toFile())

    ftpClient.retrieveFile(ftpPath, output)

    output.close()

    ftpClient.logout()
    ftpClient.disconnect()
}

private fun String.downloadHttp(timeoutMs: Int, maxRetries: Int, outputPath: Path): Boolean {
    val httpClient = HttpClientBuilder.create()
            .setDefaultRequestConfig(RequestConfig.custom().setConnectTimeout(timeoutMs).build())
            .setDefaultSocketConfig(SocketConfig.custom().setSoTimeout(timeoutMs).build())
            .setRedirectStrategy(LaxRedirectStrategy())
            .build()

    for (trial in 1..maxRetries) {
        try {
            val entity = httpClient.tryGet(this)
            entity?.content?.use {
                it.copy(outputPath, StandardCopyOption.REPLACE_EXISTING)
            }

            return true   // VICTORY!
        } catch (e: Exception) {
            if (trial == maxRetries ||
                    e !is SocketTimeoutException || e !is ConnectTimeoutException) {
                throw e
            }

            UCSC.LOG.warn("Connection timeout, retry ($trial/$maxRetries) ...")
            try {
                Thread.sleep(timeoutMs.toLong())
            } catch (ignore: InterruptedException) {
            }
        } finally {
            HttpClientUtils.closeQuietly(httpClient)
        }
    }
    return false
}

/**
 * A proxy for UCSC data.
 *
 * @author Sergei Lebedev
 */
object UCSC {
    internal val LOG = LoggerFactory.getLogger(UCSC::class.java)

    /**
     * Downloads multiple gzipped files from UCSC into a single local file.
     *
     * @param outputPath path to gzipped file.
     * @param build genome version.
     * @param chunks an array of URI components with the last component being
     *               a string template accepting chromosome name, e.g.
     *               `"%s_rmsk.txt.gz"`.
     * @throws IOException if any of the I/O operations do so.
     */
    fun downloadBatchTo(outputPath: Path, genome: Genome, rootUrl: String, template: String) {
        val tmpDir = Files.createTempDirectory("batch")
        val targetPath = tmpDir / "target.gz"

        try {
            outputPath.outputStream().use { merged ->
                for (chromosome in genome.chromosomes) {
                    "$rootUrl${template.format(chromosome.name)}".downloadTo(targetPath)
                    targetPath.inputStream().use { it.copyTo(merged) }
                }
            }
        } catch (e: IOException) {
            outputPath.deleteIfExists()  // no semi-merged files.
            throw e
        } finally {
            tmpDir.deleteDirectory()     // cleanup.
        }
    }
}
