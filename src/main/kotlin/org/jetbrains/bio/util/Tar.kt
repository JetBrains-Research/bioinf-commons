package org.jetbrains.bio.util

import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.utils.IOUtils
import org.slf4j.LoggerFactory
import java.io.File
import java.io.FileInputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.nio.file.Path

object Tar {

    private val LOG = LoggerFactory.getLogger(Tar::class.java)

    fun compress(path: Path, vararg files: File) {
        getTarArchiveOutputStream(path).use { out ->
            for (file in files) {
                addToArchiveCompression(out, file, ".")
            }
        }
    }
    
    fun decompress(input: InputStream, folder: File, skipExisting: Boolean = false) {
        folder.mkdirs()
        TarArchiveInputStream(input).use { fin ->
            while (true) {
                val entry = fin.nextTarEntry ?: break
                val outputFile = File(folder, entry.name)
                if (entry.isDirectory) {
                    if (!outputFile.exists()) {
                        val res = outputFile.mkdirs()
                        check(res) {
                            "Cannot untar file: Failed to created directory '${outputFile.path}'"
                        }
                    }
                } else {
                    if (outputFile.exists()) {
                        if (skipExisting) {
                            continue
                        } else {
                            LOG.info("File updated: $outputFile")
                            outputFile.delete()
                        }
                    }
                    FileOutputStream(outputFile).use { output ->
                        IOUtils.copy(fin, output)
                    }
                }
            }
        }
    }


    fun decompress(path: Path, folder: File) {
        decompress(FileInputStream(path.toFile()), folder)
    }

    private fun getTarArchiveOutputStream(path: Path): TarArchiveOutputStream {
        val taos = TarArchiveOutputStream(FileOutputStream(path.toFile()))
        // TAR has an 8 gig file limit by default, this gets around that
        taos.setBigNumberMode(TarArchiveOutputStream.BIGNUMBER_STAR)
        // TAR originally didn't support long file names, so enable the support for it
        taos.setLongFileMode(TarArchiveOutputStream.LONGFILE_GNU)
        taos.setAddPaxHeadersForNonAsciiNames(true)
        return taos
    }

    private fun addToArchiveCompression(out: TarArchiveOutputStream, file: File, dir: String) {
        val entry = dir + File.separator + file.name
        if (file.isFile) {
            out.putArchiveEntry(TarArchiveEntry(file, entry))
            FileInputStream(file).use { IOUtils.copy(it, out) }
            out.closeArchiveEntry()
        } else if (file.isDirectory) {
            val children = file.listFiles()
            if (children != null) {
                for (child in children) {
                    addToArchiveCompression(out, child, entry)
                }
            }
        } else {
            LOG.error("${file.name} is not supported")
        }
    }
}