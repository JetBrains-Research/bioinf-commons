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

    private var stdErr: ByteArrayOutputStream? = null
    private var stdOut: ByteArrayOutputStream? = null

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
        restoreCapturedStreams()
        if (prevSuppressExitValue == null) {
            System.getProperties().remove(JOPTSIMPLE_SUPPRESS_EXIT)
        } else {
            System.getProperties().setProperty(JOPTSIMPLE_SUPPRESS_EXIT, prevSuppressExitValue)
        }
    }

    private fun restoreCapturedStreams(): Pair<String, String> {
        System.setOut(OUT)
        System.setErr(ERR)

        val out = stdOut.toString().replace("\r", "")
        if (stdOut != null) {
            try {
                stdOut!!.close()
                stdOut = null
            } catch (e: Exception) {
                // Do nothing
            }
        }

        val err = stdErr.toString().replace("\r", "")
        if (stdErr != null) {
            try {
                stdErr!!.close()
                stdErr = null
            } catch (e: Exception) {
                // Do nothing
            }
        }

        println("STDERR: <$out>")
        println("STDOUT: <$err>")
        return out to err
    }

    @Test
    fun unrecognized() {
        with(OptionParser()) {
            parse(arrayOf("foo")) { _ ->
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertTrue("Unrecognized options: [foo]" in stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun notDumbTerminalWarning() {
        with(OptionParser()) {
            accepts("foo", "some option")

            parse(arrayOf("--foo")) { _ ->
                print("Done")
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertEquals("", stdErr)
        assertEquals("Done", stdOut)
    }

    @Test
    fun dumbTerminalWarningsOnWrongOption() {
        // Fixes current behaviour: show warning only if need to show help msg

        with(OptionParser()) {
            accepts("foo", "some option")

            parse(arrayOf("--boo")) { _ ->
                print("Done")
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertTrue("ERROR: boo is not a recognized option" in stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun dumbTerminalWarningsOnHelp() {
        // Fixes current behaviour: show warning only if need to show help msg

        with(OptionParser()) {
            accepts("foo", "some option")

            parse(arrayOf("-h")) { _ ->
                print("Done")
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertTrue("WARNING: Unable to create a system terminal, creating a dumb terminal" in stdErr)
        assertTrue("-?, -h, --help  Show help" in stdErr, stdErr)
        assertEquals("Done", stdOut, stdOut)
    }

    @Test
    fun exception() {
        with(OptionParser()) {
            parse(arrayOf()) { _ ->
                throw RuntimeException("Bah!")
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertTrue("Bah!" in stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun description() {
        with(OptionParser()) {
            accepts("foo", "some option")

            parse(arrayOf("-h"), description = "My description") { _ ->
            }
        }

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertTrue("""
            |My description
            |
            |
            |Arguments: [-h]""".trimMargin().trim() in stdErr)

        assertEquals("", stdOut)
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

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertEquals("", stdErr)
        assertEquals("", stdOut)
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

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertMatches(stdErr, "(.|\\n)*ERROR: Path (.*)foo does not exist(.|\\n)*".toRegex())
        assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("ERROR: Expected *.bed file, but was ${file.name}" in stdErr)
            assertEquals("", stdOut)
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

        val (stdOut, stdErr) = restoreCapturedStreams()
        assertMatches(stdErr, "(.|\\n)*ERROR: Path (.*)foo does not exist(.|\\n)*".toRegex())
        assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("ERROR: Expected *.bed file, but was ${file.name}" in stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("ERROR: Expected TAB separated file, but separator is [ ]" in stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("""
                |Unknown BED format:
                |chr1#2#3
                |Fields number in BED file is between 3 and 15, but was 1
                """.trimMargin().trim() in stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("""
                |Unknown BED format:
                |chr1 2
                |Fields number in BED file is between 3 and 15, but was 2
                """.trimMargin().trim() in stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("""
                           |Unknown BED format:
                           |chr1 2 0.24
                           |Fields number in BED file is between 3 and 15, but was 2
                           """.trimMargin().trim() in stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertEquals("", stdErr)
            assertEquals("", stdOut)
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

            val (stdOut, stdErr) = restoreCapturedStreams()
            assertTrue("ERROR: Bedtools will fail recognize strand column, your format is [bed4+2]" in stdErr)
            assertEquals("", stdOut)
        }
    }

    private fun assertMatches(output: String, regex: Regex) {
        assertTrue(regex.matches(output), "Doesn't match content:\n<$output>")
    }
}