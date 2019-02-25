package org.jetbrains.bio.dataframe

import org.jetbrains.bio.util.withTempFile
import org.jetbrains.bio.util.write
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

private val testDataFrameSpec = DataFrameSpec()
        .ints("i").shorts("s").bytes("b").longs("l").bools("f").strings("str")
        .enums("e", TestEnum::class.java)

private val testDataFrame: DataFrame
    get() {
    val bitSet = BitterSet(5)
    bitSet.set(2)
    bitSet.set(3)
    return DataFrame().with("i", intArrayOf(5, 1, 3, 4, 2))
            .with("s", shortArrayOf(5, 1, 3, 4, 2))
            .with("b", byteArrayOf(5, 1, 3, 4, 2))
            .with("l", longArrayOf(5, 1, 3, 4, 2))
            .with("f", bitSet)
            .with("str", arrayOf("str5", "str1", "str3", "str4", "str2"))
            .with("e", TestEnum::class.java, arrayOf(TestEnum.VAL5, TestEnum.VAL1,
                    TestEnum.VAL3, TestEnum.VAL4,
                    TestEnum.VAL2))
}

class CsvMapperTest {
    // TODO: csv

    @Test fun testGuess() {
        val blob = arrayOf("# Integer; Short; Byte; Long; Boolean; String; org.jetbrains.bio.dataframe.TestEnum",
                           "i\ts\tb\tl\tf\tstr\te",
                           "5\t5\t5\t5\t0\tstr5\tVAL5",
                           "1\t1\t1\t1\t0\tstr1\tVAL1",
                           "3\t3\t3\t3\t1\tstr3\tVAL3",
                           "4\t4\t4\t4\t1\tstr4\tVAL4",
                           "2\t2\t2\t2\t0\tstr2\tVAL2").joinToString("\n")

        withTempFile("df", ".csv") { path ->
            path.write(blob)

            val spec = DataFrameMappers.TSV.guess(path)
            assertTrue(spec.columns[0] is IntColumn)
            assertTrue(spec.columns[1] is ShortColumn)
            assertTrue(spec.columns[2] is ByteColumn)
            assertTrue(spec.columns[3] is LongColumn)
            assertTrue(spec.columns[4] is BooleanColumn)
            assertTrue(spec.columns[5] is StringColumn)
            assertTrue(spec.columns[6] is EnumColumn<*>)
            assertEquals(TestEnum::class.java, spec.columns[6].boxedType)
        }
    }

    @Test fun testSaveLoad() {
        val suffixes = arrayOf(".csv.gz", ".csv")
        for (suffix in suffixes) {
            withTempFile("df", suffix) { path ->
                val df = testDataFrame
                df.save(path)
                val loaded = DataFrameMappers.TSV.load(path, testDataFrameSpec)

                assertEquals(df, loaded, message = suffix)
            }
        }
    }

    @Test fun testSaveHeader() {
        val df = testDataFrame
        withTempFile("df", ".csv") { path ->
            DataFrameMappers.TSV.save(path, df)

            val blob = listOf("# Integer; Short; Byte; Long; Boolean; String; org.jetbrains.bio.dataframe.TestEnum",
                              "i\ts\tb\tl\tf\tstr\te",
                              "5\t5\t5\t5\t0\tstr5\tVAL5",
                              "1\t1\t1\t1\t0\tstr1\tVAL1",
                              "3\t3\t3\t3\t1\tstr3\tVAL3",
                              "4\t4\t4\t4\t1\tstr4\tVAL4",
                              "2\t2\t2\t2\t0\tstr2\tVAL2")
            assertEquals(blob, path.toFile().readLines())
        }
    }

    @Test fun testLoadNoHeader() {
        val blob = arrayOf("5\t5\t5\t5\t0\tstr5",
                           "1\t1\t1\t1\t0\tstr1",
                           "3\t3\t3\t3\t1\tstr3",
                           "4\t4\t4\t4\t1\tstr4",
                           "2\t2\t2\t2\t0\tstr2").joinToString("\n")

        withTempFile("df", ".csv") { path ->
            path.write(blob)

            val spec = DataFrameSpec().ints("i").shorts("s").bytes("b").longs("l").bools("f").strings("str")
            assertEquals(testDataFrame.omit("e"),
                    DataFrameMappers.TSV.load(path, spec, header = false))
        }
    }

    @Test fun testLoadNonNumericQuoting() {
        // backward compatibility with previous CSV format: quote non-numeric
        val blob = arrayOf("\"i\"\t\"s\"\t\"b\"\t\"l\"\t\"f\"\t\"str\"",
                           "5\t5\t5\t5\t0\t\"str5\"",
                           "1\t1\t1\t1\t0\t\"str1\"",
                           "3\t3\t3\t3\t1\t\"str3\"",
                           "4\t4\t4\t4\t1\t\"str4\"",
                           "2\t2\t2\t2\t0\t\"str2\"").joinToString("\n")

        withTempFile("df", ".csv") { path ->
            path.write(blob)

            val spec = DataFrameSpec().ints("i").shorts("s").bytes("b").longs("l").bools("f").strings("str")
            assertEquals(testDataFrame.omit("e"), DataFrameMappers.TSV.load(path, spec))
        }
    }

    @Test fun testDump() {
        val expected = arrayOf("# Integer; Short; Byte; Long; Boolean; String; org.jetbrains.bio.dataframe.TestEnum",
                               "i\ts\tb\tl\tf\tstr\te",
                               "5\t5\t5\t5\t0\tstr5\tVAL5",
                               "1\t1\t1\t1\t0\tstr1\tVAL1",
                               "3\t3\t3\t3\t1\tstr3\tVAL3",
                               "4\t4\t4\t4\t1\tstr4\tVAL4",
                               "2\t2\t2\t2\t0\tstr2\tVAL2",
                               "").joinToString("\r\n")

        val buff = StringBuilder()
        DataFrameMappers.TSV.save(buff, testDataFrame, true, true);
        assertEquals(expected, buff.toString())
    }

    @Test fun testDumpHead() {
        val dump = testDataFrame.dumpHead(3)
        val expected = arrayOf("# Integer; Short; Byte; Long; Boolean; String; org.jetbrains.bio.dataframe.TestEnum",
                               "i\ts\tb\tl\tf\tstr\te",
                               "5\t5\t5\t5\t0\tstr5\tVAL5",
                               "1\t1\t1\t1\t0\tstr1\tVAL1",
                               "3\t3\t3\t3\t1\tstr3\tVAL3",
                               "").joinToString("\r\n")
        assertEquals(expected, dump)
    }
}

class NpzMapperTest {
    private val specNoEnum = DataFrameSpec()
            .ints("i").shorts("s").bytes("b").longs("l").bools("f").strings("str")

    @Test fun testGuess() {
        withTempFile("df", ".npz") { path ->
            val df = testDataFrame.omit("e")
            DataFrameMappers.NPZ.save(path, df)
            assertEquals(specNoEnum, DataFrameMappers.NPZ.guess(path))
        }
    }

    @Test fun testSaveLoad() {
        withTempFile("df", ".npz") { path ->
            val df = testDataFrame.omit("e")
            DataFrameMappers.NPZ.save(path, df)
            assertEquals(df, DataFrameMappers.NPZ.load(path, specNoEnum))
        }
    }
}
