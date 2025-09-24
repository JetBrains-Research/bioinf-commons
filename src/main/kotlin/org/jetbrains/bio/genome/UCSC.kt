package org.jetbrains.bio.genome

import org.apache.commons.net.ftp.FTPClient
import org.apache.hc.client5.http.ConnectTimeoutException
import org.apache.hc.client5.http.HttpResponseException
import org.apache.hc.client5.http.classic.methods.HttpGet
import org.apache.hc.client5.http.config.RequestConfig
import org.apache.hc.client5.http.impl.DefaultRedirectStrategy
import org.apache.hc.client5.http.impl.classic.CloseableHttpClient
import org.apache.hc.client5.http.impl.classic.HttpClientBuilder
import org.apache.hc.core5.http.HttpEntity
import org.apache.hc.core5.io.CloseMode
import org.apache.hc.core5.util.Timeout
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.FileOutputStream
import java.io.IOException
import java.net.SocketTimeoutException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.concurrent.TimeUnit

private fun CloseableHttpClient.tryGet(url: String): HttpEntity? {
    val response = execute(HttpGet(url))
    val code = response.code
    if (code / 100 != 2) {
        throw HttpResponseException(code, "${response.reasonPhrase}: $url")
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
fun String.downloadTo(
    outputPath: Path,
    timeout: Int = 30,
    maxRetries: Int = 10
) {
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

    FileOutputStream(outputPath.toFile()).use { output ->
        ftpClient.retrieveFile(ftpPath, output)
    }

    ftpClient.logout()
    ftpClient.disconnect()
}

private fun String.downloadHttp(timeoutMs: Int, maxRetries: Int, outputPath: Path): Boolean {
    val httpClient = HttpClientBuilder.create()
        .setDefaultRequestConfig(RequestConfig.custom().setConnectTimeout(Timeout.of(timeoutMs.toLong(), TimeUnit.MILLISECONDS)).build())
        .setRedirectStrategy(DefaultRedirectStrategy())
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
                e !is SocketTimeoutException || e !is ConnectTimeoutException
            ) {
                throw e
            }

            UCSC.LOG.warn("Connection timeout, retry ($trial/$maxRetries) ...")
            try {
                Thread.sleep(timeoutMs.toLong())
            } catch (ignore: InterruptedException) {
            }
        } finally {
            httpClient.close(CloseMode.GRACEFUL)
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
