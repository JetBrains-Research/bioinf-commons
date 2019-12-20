package org.jetbrains.bio.util

import joptsimple.OptionParser
import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertMatches
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.nio.file.Path
import java.nio.file.Paths
import kotlin.test.assertEquals
import kotlin.test.assertNull

class OptionParserExtensionsTest {
    private var prevSuppressExitValue: String? = null

    @Before
    fun setUp() {
        prevSuppressExitValue = System.getProperty(JOPTSIMPLE_SUPPRESS_EXIT)
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true")
    }

    @After
    fun tearDown() {
        if (prevSuppressExitValue == null) {
            System.getProperties().remove(JOPTSIMPLE_SUPPRESS_EXIT)
        } else {
            System.getProperties().setProperty(JOPTSIMPLE_SUPPRESS_EXIT, prevSuppressExitValue)
        }
    }

    @Test
    fun unrecognized() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
            with(OptionParser()) {
                parse(arrayOf("foo")) { _ ->
                }
            }
        }
        assertIn("Unrecognized options: [foo]", stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun exception() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
            with(OptionParser()) {
                parse(arrayOf()) { _ ->
                    throw RuntimeException("Bah!")
                }
            }
        }
        assertIn("Bah!", stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun description() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
            with(OptionParser()) {
                accepts("foo", "some option")

                parse(arrayOf("--foo", "-h"), description = "My description") { _ ->
                }
            }
        }

        assertIn("My description", stdErr)
        assertIn("Arguments:\n    --foo\n    -h", stdErr)

        assertEquals("", stdOut)
    }


    @Test
    fun pathConverterNoCheck() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
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
        }

        assertEquals("", stdErr)
        assertEquals("", stdOut)
    }

    @Test
    fun pathConverterExists_MissingFile() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
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
        }
        assertMatches(stdErr, "(.|\\n|\\r\\n)*ERROR: Path (.*)foo does not exist(.|\\n|\\r\\n)*".toRegex())
        assertEquals("", stdOut)
    }

    @Test
    fun pathConverterExists() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.exists())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun pathConverterExists_Ext() {
        withTempFile("optParser", suffix = "suffix.bed") { file ->
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.exists("bed"))
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }


    @Test
    fun pathConverterExists_ExtOtherCase() {
        withTempFile("optParser", suffix = "suffix.bEd") { file ->
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.exists("bed"))
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun pathConverterExists_BadExt() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.exists("bed"))
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("ERROR: Expected *.bed file, but was ${file.name}", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_MissingFile() {
        val (stdOut, stdErr) = Logs.captureLoggingOutput {
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
        }

        assertMatches(stdErr, "(.|\\n|\\r\\n)*ERROR: Path (.*)foo does not exist(.|\\n|\\r\\n)*".toRegex())
        assertEquals("", stdOut)
    }

    @Test
    fun bedtoolsValidFile_BadExt() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile("bed"))
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("ERROR: Expected *.bed file, but was ${file.name}", stdErr)
            assertEquals("", stdOut)
        }
    }


    @Test
    fun bedtoolsValidFile_NotTabbed() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1 2 3")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("ERROR: Expected TAB separated file, but separator is [ ]", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_UnsupportedSeparator() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1#2#3\nchr1#3#4")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("Unknown BED format:", stdErr)
            assertIn("chr1#3#4", stdErr)
            assertIn("Fields number in BED file is between 3 and 15, but was 1", stdErr)

            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_FieldsLess3() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\nchr1\t3")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("Unknown BED format:", stdErr)
            assertIn("chr1\t3", stdErr)
            assertIn("Fields number in BED file is between 3 and 15, but was 2", stdErr)

            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_BedFieldsLess3() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t0.24\nchr1\t3\t0.42")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("Unknown BED format:", stdErr)
            assertIn("chr1\t3\t0.42", stdErr)
            assertIn("Fields number in BED file is between 3 and 15, but was 2", stdErr)

            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_Bed4() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\t0.24")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_Bed6() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t24\t.")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_Bed6p1() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t24\t+\t0.33")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_Bed5() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t200")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertEquals(file, options.valueOf("path") as Path)
                    }
                }
            }

            assertEquals("", stdErr)
            assertEquals("", stdOut)
        }
    }

    @Test
    fun bedtoolsValidFile_BedFieldsUnknownStrand() {
        withTempFile("optParser", suffix = "suffix.txt") { file ->
            file.write("chr1\t2\t3\tfoo\t0.24\t0.33")
            require(file.exists)

            val (stdOut, stdErr) = Logs.captureLoggingOutput {
                with(OptionParser()) {
                    acceptsAll(listOf("path"))
                            .withRequiredArg()
                            .withValuesConvertedBy(PathConverter.bedtoolsValidFile())
                            .required()
                    parse(arrayOf("--path", file.toString())) { options ->
                        assertNull(options.valueOf("path"))
                    }
                }
            }

            assertIn("ERROR: Bedtools will fail to recognize the strand column, your format is [bed4+]", stdErr)
            assertEquals("", stdOut)
        }
    }
}