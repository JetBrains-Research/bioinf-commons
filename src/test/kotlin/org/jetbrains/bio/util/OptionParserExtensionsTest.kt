package org.jetbrains.bio.util

import joptsimple.OptionParser
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.PrintStream
import java.nio.file.Path
import java.nio.file.Paths
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

class OptionParserExtensionsTest {
    private var prevSuppressExitValue: String? = null

    companion object {
        private val OUT = System.out
        private val ERR = System.err
    }

    private lateinit var stdErr: ByteArrayOutputStream
    private lateinit var stdOut: ByteArrayOutputStream

    @Before
    fun setUp() {
        stdErr = ByteArrayOutputStream()
        stdOut = ByteArrayOutputStream()
        System.setOut(PrintStream(stdOut))
        System.setErr(PrintStream(stdErr))

        prevSuppressExitValue = System.getProperty(JOPTSIMPLE_SUPPRESS_EXIT)
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true")
    }

    @After
    fun tearDown() {
        System.setOut(OUT)
        System.setErr(ERR)
        try {
            stdOut.close()
        } catch (e: Exception) {
            // Do nothing
        }
        try {
            stdErr.close()
        } catch (e: Exception) {
            // Do nothing
        }
        if (prevSuppressExitValue == null) {
            System.getProperties().remove(JOPTSIMPLE_SUPPRESS_EXIT)
        } else {
            System.getProperties().setProperty(JOPTSIMPLE_SUPPRESS_EXIT, prevSuppressExitValue)
        }
    }
    @Test
    fun unrecognized() {
        with(OptionParser()) {
            parse(arrayOf("foo")) { _ ->
            }
        }
        assert("Unrecognized options: [foo]" in stdErr.toString())
        assertEquals("", stdOut.toString())
    }

    @Test
    fun exception() {
        with(OptionParser()) {
            parse(arrayOf()) { _ ->
                throw RuntimeException("Bah!")
            }
        }

        assert("Bah!" in stdErr.toString())
        assertEquals("", stdOut.toString())
    }

    @Test
    fun description() {
            with(OptionParser()) {
                accepts("foo", "some option")

                parse(arrayOf("-h"), description = "My description") { _ ->
                }
            }
        val stdErrText = stdErr.toString().replace("\r", "")
        assertTrue("""
            |My description
            |
            |
            |Arguments: [-h]""".trimMargin().trim() in stdErrText)
        
        assertEquals("", stdOut.toString())
    }


    @Test
    fun pathConverterNoCheck() {
        with(OptionParser()) {
            acceptsAll(listOf("path"))
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.noCheck())
                    .required()
            parse(arrayOf("--path", "foo")) { options ->
                val path = options.valueOf("path") as Path
                require(path.notExists)
                assertEquals("foo", path.name)
            }
        }
        assertEquals("", stdErr.toString())
        assertEquals("", stdOut.toString())
    }

    @Test
    fun pathConverterExists_MissingFile() {
        with(OptionParser()) {
            acceptsAll(listOf("path"))
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())
                    .required()
            parse(arrayOf("--path", "foo")) { options ->
                require(Paths.get("foo").notExists)
                assertNull(options.valueOf("path"))
            }
        }

        assertMatches(stdErr, "ERROR: Path (.*)foo does not exist(.|\\n)*".toRegex())
        assertEquals("", stdOut.toString())
    }

    @Test
    fun pathConverterExists() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.exists())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun pathConverterExists_Ext() {
        withTempFile("optParser", suffix = "suffix.bed") { file ->
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.exists("bed"))
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }


    @Test
    fun pathConverterExists_ExtOtherCase() {
        withTempFile("optParser", suffix = "suffix.bEd") { file ->
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.exists("bed"))
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun pathConverterExists_BadExt() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.exists("bed"))
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("ERROR: Expected *.bed file, but was ${file.name}" in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_MissingFile() {
        with(OptionParser()) {
            acceptsAll(listOf("path"))
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                    .required()
            parse(arrayOf("--path", "foo")) { options ->
                require(Paths.get("foo").notExists)
                assertNull(options.valueOf("path"))
            }
        }

        assertMatches(stdErr, "ERROR: Path (.*)foo does not exist(.|\\n)*".toRegex())
        assertEquals("", stdOut.toString())
    }

    @Test
    fun bedtoolsValidFile_BadExt() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile("bed"))
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("ERROR: Expected *.bed file, but was ${file.name}" in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }


    @Test
    fun bedtoolsValidFile_NotTabbed() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1 2 3")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("ERROR: Expected TAB separated file, but separator is [ ]" in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_UnsupportedSeparator() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1#2#3")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("""
                |Unknown BED format:
                |chr1#2#3
                |Fields number in BED file is between 3 and 15, but was 1
                """.trimMargin().trim() in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_FieldsLess3() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1 2")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("""
                |Unknown BED format:
                |chr1 2
                |Fields number in BED file is between 3 and 15, but was 2
                """.trimMargin().trim() in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_BedFieldsLess3() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1 2 0.24")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("""
                           |Unknown BED format:
                           |chr1 2 0.24
                           |Fields number in BED file is between 3 and 15, but was 2
                           """.trimMargin().trim() in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_Bed4() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\t0.24")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_Bed6() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t24\t.")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_Bed6p1() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t24\t+\t0.33")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_Bed5() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t200")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertEquals(file, options.valueOf("path") as Path)
                }
            }

            assertEquals("", stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    @Test
    fun bedtoolsValidFile_BedFieldsUnknownStrand() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t0.24\t0.33")
            require(file.exists)

            with(OptionParser()) {
                acceptsAll(listOf("path"))
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                        .required()
                parse(arrayOf("--path", file.toString())) { options ->
                    assertNull(options.valueOf("path"))
                }
            }

            assertTrue("ERROR: Bedtools will fail recognize strand column, your format is [bed4+2]" in stdErr.toString())
            assertEquals("", stdOut.toString())
        }
    }

    private fun assertMatches(outputStream: ByteArrayOutputStream, regex: Regex) {
        outputStream.toString().let { content ->
            assertTrue(
                    regex.matches(content.replace("\r", "")),
                    "Doesn't match content:\n<$content>"
            )
        }
    }

}