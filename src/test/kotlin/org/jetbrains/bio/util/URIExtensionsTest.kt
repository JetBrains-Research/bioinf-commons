package org.jetbrains.bio.util

import org.jetbrains.bio.Tests
import org.junit.Assume
import org.junit.Test
import java.io.File
import java.io.InputStreamReader
import java.net.URI
import java.nio.file.Paths
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * @author Roman.Chernyatchik
 */
class URIExtensionsTest {

    @Test
    fun isFile() {
        assertTrue(URI.create("file:///mnt/stripe/foo.bw").isFile())
        assertTrue(URI.create("file:/C:/mnt/stripe/foo.bw").isFile())
        assertTrue(URI.create("file:///C:/mnt/stripe/foo.bw").isFile())

        assertFalse(URI.create("ftp://hgdownload.soe.ucsc.edu/goldenPath/foo.bw").isFile())
        assertFalse(URI.create("https://www.encodeproject.org/files/@@download/foo.bigWig").isFile())
    }

    @Test
    fun extension() {
        assertEquals("omni", "file:///mnt/stripe/foo.omni".toUri().extension())
        assertEquals("span", "file:///mnt/stripe/foo.span".toUri().extension())
        assertEquals("narrowPeak", "file:///mnt/stripe/f.o.o.narrowPeak".toUri().extension())
        assertEquals("bigWig", ("https://www.encodeproject.org/files/@@download/foo.bigWig").toUri().extension())
        assertEquals("bw", "ftp://hgdownload.soe.ucsc.edu/goldenPath/foo.bw".toUri().extension())
        assertEquals("bw", "file:///mnt/stripe/boo%20doo%20foo.bw".toUri().extension())
        assertEquals("bigBed", "file:/C:/mnt/stripe/foo.bigBed".toUri().extension())
        assertEquals("gz", "file:/C:/mnt/stripe/foo.tsv.gz".toUri().extension())
        assertEquals("tsv", "file:///C:/mnt/stripe/foo.tsv".toUri().extension())
        assertEquals("peak", "file:/C:/mnt/stripe/boo%20foo.peak".toUri().extension())

        // Decode url:
        assertEquals("bw", "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498%5FCC%5FGVO%5FRNAseq%5Frep1%2Ebw".toUri().extension())

        // Not extension:
        assertEquals("org?request_id=foo", "https://www.doo.org?request_id=foo".toUri().extension())
        assertEquals("", "file:///mnt/stripe/foo".toUri().extension())
        assertEquals("", "file:/C:/mnt/stripe/bigBed".toUri().extension())
    }
    @Test
    fun isSupportedURL() {
        assertTrue(URI.create("ftp://hgdownload.soe.ucsc.edu/goldenPath/foo.bw").isSupportedURL())
        assertTrue(URI.create("https://www.encodeproject.org/files/@@download/foo.bigWig").isSupportedURL())
        assertTrue(URI.create("http://www.encodeproject.org/files/@@download/foo.bigWig").isSupportedURL())

        assertFalse(URI.create("file:///mnt/stripe/foo.bw").isSupportedURL())
        assertFalse(URI.create("file:/C:/mnt/stripe/foo.bw").isSupportedURL())
        assertFalse(URI.create("file:///C:/mnt/stripe/foo.bw").isSupportedURL())
        assertFalse(URI.create("mailto:java-net@java.sun.com").isSupportedURL())
        assertFalse(URI.create("telnet://user:password@host.com").isSupportedURL())
        assertFalse(URI.create("//localhost/index.html").isSupportedURL())
        assertFalse(URI.create("nntp://host.com:119/newsgroup").isSupportedURL())
    }

    @Test
    fun toPath() {
        assertEquals("foo.bw", URI.create("file:///mnt/stripe/foo.bw").toPath().name)
        assertEquals("boo doo foo.bw", URI.create("file:///mnt/stripe/boo%20doo%20foo.bw").toPath().name)
        assertEquals("boo doo foo.bw", "file:///mnt/stripe/boo doo foo.bw".toUri().toPath().name)

        if (isWindows()) {
            assertEquals(
                "C:\\mnt\\stripe\\foo.bw",
                URI.create("file:///C:/mnt/stripe/foo.bw").toPath().toString()
            )
            assertEquals(
                "C:\\mnt\\stripe\\foo.bw",
                URI.create("file:/C:/mnt/stripe/foo.bw").toPath().toString()
            )
            assertEquals(
                "C:\\mnt\\stripe\\boo foo.bw",
                URI.create("file:/C:/mnt/stripe/boo%20foo.bw").toPath().toString()
            )
        } else {
            assertEquals(
                "/mnt/stripe/foo.bw",
                URI.create("file:///mnt/stripe/foo.bw").toPath().toString()
            )
            assertEquals(
                "/mnt/stripe/foo.bw",
                URI.create("file:/mnt/stripe/foo.bw").toPath().toString()
            )
            assertEquals(
                "/mnt/stripe/boo foo.bw",
                URI.create("file:/mnt/stripe/boo%20foo.bw").toPath().toString()
            )
        }
    }

    @Test
    fun toPathForNonFileUrl() {
        val url = "http://www.encodeproject.org/files/@@download/foo.bigWig"

        Tests.assertThrowsWithMessage(
            "Cannot convert URL to path: $url",
            IllegalArgumentException::class.java,
        ) {
            URI.create(url).toPath()
        }
    }

    @Test
    fun toPathForUrlWithoutProtocol() {
        assertEquals("/foo/boo".replace("/", File.separator), URI.create("/foo/boo").toPath().toString())
    }

    @Test
    fun presentablePath() {
        if (isWindows()) {
            assertEquals(
                "C:\\mnt\\stripe\\foo.bw",
                URI.create("file:/C:/mnt/stripe/foo.bw").presentablePath()
            )
            assertEquals(
                "C:\\mnt\\stripe\\foo.bw",
                URI.create("file:///C:/mnt/stripe/foo.bw").presentablePath()
            )
        } else {
            assertEquals(
                "/mnt/stripe/foo.bw",
                URI.create("file:///mnt/stripe/foo.bw").presentablePath()
            )
        }

        assertEquals(
            "https://www.encodeproject.org/files/@@download/foo.bigWig",
            URI.create("https://www.encodeproject.org/files/@@download/foo.bigWig").presentablePath()
        )
        assertEquals(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498_CC_GVO_RNAseq_rep1.bw",
            URI.create("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498%5FCC%5FGVO%5FRNAseq%5Frep1%2Ebw")
                .presentablePath()
        )
    }

    @Test
    fun shortNameFile() {
        withTempFile("foo", "boo") { path ->
            assertEquals(path.name, path.toUri().shortName())
        }
    }

    @Test
    fun shortNameURL() {
        assertEquals(
            "foo.bigWig",
            "https://www.doo.org/foo.bigWig".toUri().shortName()
        )
        assertEquals(
            "foo.bigWig",
            "https://www.doo.org/aaa/foo.bigWig".toUri().shortName()
        )

        assertEquals(
            "https://www.doo.org",
            "https://www.doo.org".toUri().shortName()
        )
        assertEquals(
            "https://www.doo.org/foo.bigWig?request_id=foo",
            "https://www.doo.org/foo.bigWig?request_id=foo".toUri().shortName()
        )
        assertEquals(
            "https://www.doo.org?request_id=foo",
            "https://www.doo.org?request_id=foo".toUri().shortName()
        )
        assertEquals(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498_CC_GVO_RNAseq_rep1.bw",
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498%5FCC%5FGVO%5FRNAseq%5Frep1%2Ebw".toUri()
                .shortName()
        )
    }

    @Test
    fun toURI() {
        if (isWindows()) {
            assertEquals(URI.create("file:///C:/mnt/stripe/foo%20foo.bw"), "C:/mnt/stripe/foo foo.bw".toUri())
            assertEquals(URI.create("file:///C:/mnt/stripe/foo.bw"), "C:/mnt/stripe/foo.bw".toUri())
            assertEquals(URI.create("file:///C:/mnt/stripe/foo.bw"), "C:\\mnt\\stripe\\foo.bw".toUri())
            assertEquals(
                URI.create(
                    "file:///${
                        Paths.get(".").toAbsolutePath().normalize().toString().replace(File.separatorChar, '/')
                    }/stripe/foo.bw"
                ),
                "stripe\\foo.bw".toUri()
            )
            assertEquals(URI.create("file:///C:/mnt/stripe/foo.bw"), "file:///C:/mnt/stripe/foo.bw".toUri())
            assertEquals(URI.create("file:///C:/mnt/stripe/foo.bw"), "file:/C:/mnt/stripe/foo.bw".toUri())
            //assertEquals(URI.create("file:///C:/mnt/stripe/foo.bw"), "file://C:/mnt/stripe/foo.bw".toUri())
        } else {
            assertEquals(URI.create("file:///mnt/stripe/foo.bw"), "/mnt/stripe/foo.bw".toUri())
            assertEquals(URI.create("file:///mnt/stripe/foo%20foo.bw"), "/mnt/stripe/foo foo.bw".toUri())
            assertEquals(
                URI.create("file://${Paths.get(".").toAbsolutePath().normalize()}/stripe/foo.bw"),
                "stripe/foo.bw".toUri()
            )
            assertEquals(URI.create("file:///mnt/stripe/boo%20foo.bw"), "file:/mnt/stripe/boo foo.bw".toUri())
        }
        assertEquals(URI.create("file:///mnt/stripe/foo.bw"), "file:///mnt/stripe/foo.bw".toUri())
        assertEquals(URI.create("file:///mnt/stripe/foo%20foo.bw"), "file:///mnt/stripe/foo foo.bw".toUri())
        assertEquals(URI.create("file:///mnt/stripe/boo%20foo%20foo.bw"), "file:///mnt/stripe/boo foo foo.bw".toUri())
        assertEquals(URI.create("https://www.doo.org/foo.bigWig"), "https://www.doo.org/foo.bigWig".toUri())
    }

    @Test
    fun toURIMalformedPath() {
        // if from path with %20 - path will be malformed
        if (isWindows()) {
            assertEquals("file:///C:/foo/boo%2520doo", "C:\\foo\\boo%20doo".toUri().toString())

            val driveLetter = Paths.get(".").toAbsolutePath().toString()[0]
            assertEquals("file:///$driveLetter:/foo/boo%2520doo", "/foo/boo%20doo".toUri().toString())
            assertEquals("$driveLetter:\\foo\\boo%20doo", "/foo/boo%20doo".toUri().toPath().toString())
        } else {
            assertEquals("file:///foo/boo%2520doo", "/foo/boo%20doo".toUri().toString())
            assertEquals("/foo/boo%20doo", "/foo/boo%20doo".toUri().toPath().toString())
        }
    }

    @Test(expected = IllegalArgumentException::class)
    fun toURIMalFormedURL() {
        "https://www.doo.org/boo foo.bigWig".toUri()
    }

    @Test
    fun toURIMalFormedPathWindows() {
        Assume.assumeTrue(isWindows())

        Tests.assertThrowsWithMessage(
            "Illegal character in path at index 10: file:///C:\\mnt\\stripe\\foo.bw",
            IllegalArgumentException::class.java,
        ) {
            "file:///C:\\mnt\\stripe\\foo.bw".toUri()
        }


    }

    @Test
    fun checkAccessibleFile() {
        withTempFile("foo", "boo") { path ->
            path.toUri().checkAccessible()
        }
    }

    @Test
    fun isAccessibleFile() {
        withTempFile("foo", "boo") { path ->
            assertTrue(path.toUri().isAccessible())
        }
    }

    @Test
    fun checkAccessibleFileNotExist() {
        withTempFile("foo", "boo") { path ->
            path.deleteIfExists()

            Tests.assertThrowsWithMessage(
                "Track file doesn't exist: $path",
                IllegalStateException::class.java,
            ) {
                path.toUri().checkAccessible()
            }
        }
    }

    @Test
    fun isAccessibleFileNotExist() {
        withTempFile("foo", "boo") { path ->
            path.deleteIfExists()
            assertFalse(path.toUri().isAccessible())
        }
    }

    @Test
    fun checkAccessibleURI() {
        Tests.assertThrowsWithMessage(
            "URL not supported: foo://boo",
            IllegalStateException::class.java,
        ) {
            "foo://boo".toUri().checkAccessible()
        }

    }

    @Test
    fun isAccessibleURI() {
        assertFalse("foo://boo".toUri().isAccessible())
    }

    @Suppress("SpellCheckingInspection")
    @Test
    fun hasExtURL() {
        assertTrue("https://github.com/blob/master/ENCFF575VMI.bigBed".toUri().hasExt("bigBed"))
        assertTrue("https://github.com/blob/master/ENCFF575VMI.bigBed".toUri().hasExt("bigbed"))
        assertTrue("https://github.com/blob/master/ENCFF575VMI.bigbed".toUri().hasExt("bigBed"))
        assertTrue(
            ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63874&format=file&file=" +
                    "GSE63874%5FNa%5FK2%5Fall%5Fminus%5Fcov.tdf").toUri().hasExt("tdf")
        )
        assertTrue(
            ("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/" +
                    "GSE63874%5FNa%5FK2%5Fall%5Fminus%5Fcov%2Etdf").toUri().hasExt("tdf")
        )

        assertTrue(
            ("https://github.com/PetrTsurinov/BigBedTest/blob/master/" +
                    "ENCFF575VMI.bigBed?raw=true").toUri().hasExt("bigBed")
        )

        assertFalse("https://github.com/ENCFF575VMI.bed".toUri().hasExt("bam", "foo"))
        assertTrue("https://github.com/ENCFF575VMI.bed".toUri().hasExt("bam", "bed", "foo"))

        assertTrue(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3074498&format=file&file=GSM3074498%5FCC%5FGVO%5FRNAseq%5Frep1%2Ebw".toUri()
                .hasExt("bw")
        )
        assertTrue(
            "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3074nnn/GSM3074498/suppl/GSM3074498%5FCC%5FGVO%5FRNAseq%5Frep1%2Ebw".toUri()
                .hasExt("bw")
        )

    }

    @Suppress("SpellCheckingInspection")
    @Test
    fun hasExtFile() {

        assertFalse("/foo/boo/boobigBed".toUri().hasExt("bigBed"))
        assertTrue("/foo/boo/doo.bigBed".toUri().hasExt("bigBed"))
        assertFalse("/foo/boo/foo".toUri().hasExt("bigBed"))

        assertTrue("/foo/boo/boo.bed.gz".toUri().hasExt("gz"))
        assertTrue("/foo/boo/boo.bed.gz".toUri().hasExt("bed.gz"))

        assertTrue("C:/foo/boo/doo.bigBed".toUri().hasExt("bigBed"))
        assertTrue("C:///foo/boo/doo.bigBed".toUri().hasExt("bigBed"))
    }

    @Test
    fun asByteSource() {
        withTempFile("foo", ".boo") { path ->
            val content = "line1\nline2\nline3"
            path.write(content)
            val reader = InputStreamReader(path.toUri().asByteSource().openBufferedStream(), Charsets.UTF_8)
            assertEquals(content, reader.useLines { it.joinToString("\n") })
        }
    }
}