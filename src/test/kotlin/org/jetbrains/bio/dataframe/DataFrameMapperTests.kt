package org.jetbrains.bio.dataframe

import org.jetbrains.bio.util.withLocaleEU
import org.jetbrains.bio.util.withLocaleUS
import org.jetbrains.bio.util.withTempFile
import org.jetbrains.bio.util.write
import org.junit.Assert
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

private val testDataFrameSpec = DataFrameSpec()
    .ints("i").shorts("s").bytes("b").longs("l").booleans("f").strings("str")
    .enums("e", TestEnum::class.java)

private val testDataFrame: DataFrame
    get() {
        val bitList = BitList(5)
        bitList.set(2)
        bitList.set(3)
        return DataFrame().with("i", intArrayOf(5000000, 1000000, 3000000, 4000000, 2000000))
            .with("s", shortArrayOf(5000, 1000, 3000, 4000, 2000))
            .with("b", byteArrayOf(5, 1, 3, 4, 2))
            .with("l", longArrayOf(5000000, 1000000, 3000000, 4000000, 2000000))
            .with("f", bitList)
            .with("str", arrayOf("str5", "str1", "str3", "str4", "str2"))
            .with(
                "e", TestEnum::class.java, arrayOf(
                    TestEnum.VAL5, TestEnum.VAL1,
                    TestEnum.VAL3, TestEnum.VAL4,
                    TestEnum.VAL2
                )
            )
    }

class CsvMapperTest {
    @Test
    fun testGuess() {
        val blob = arrayOf(
            "# Integer;\tShort;\tByte;\tLong;\tBoolean;\tString;\torg.jetbrains.bio.dataframe.TestEnum",
            "i\ts\tb\tl\tf\tstr\te",
            "5000000\t5000\t5\t5000000\t0\tstr5\tVAL5",
            "1000000\t1000\t1\t1000000\t0\tstr1\tVAL1",
            "3\t3\t3\t3\t1\tstr3\tVAL3",
            "4000000\t4000\t4\t4000000\t1\tstr4\tVAL4",
            "2000000\t2000\t2\t2000000\t0\tstr2\tVAL2"
        ).joinToString("\n")

        withTempFile("df", ".tsv") { path ->
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

    @Test
    fun testSaveLoad() {
        val suffixes = arrayOf(".tsv.gz", ".tsv")
        for (suffix in suffixes) {
            withTempFile("df", suffix) { path ->
                val df = testDataFrame
                df.save(path)
                val loaded = DataFrameMappers.TSV.load(path, testDataFrameSpec)

                assertEquals(df, loaded, message = suffix)
            }
        }
    }

    @Test
    fun testSaveHeader() {
        val df = testDataFrame
        withTempFile("df", ".tsv") { path ->
            DataFrameMappers.TSV.save(path, df)

            val blob = listOf(
                "# Integer;\tShort;\tByte;\tLong;\tBoolean;\tString;\torg.jetbrains.bio.dataframe.TestEnum",
                "i\ts\tb\tl\tf\tstr\te",
                "5000000\t5000\t5\t5000000\t0\tstr5\tVAL5",
                "1000000\t1000\t1\t1000000\t0\tstr1\tVAL1",
                "3000000\t3000\t3\t3000000\t1\tstr3\tVAL3",
                "4000000\t4000\t4\t4000000\t1\tstr4\tVAL4",
                "2000000\t2000\t2\t2000000\t0\tstr2\tVAL2"
            )
            assertEquals(blob, path.toFile().readLines())
        }
    }

    @Test
    fun testLoadNoHeader() {
        val blob = arrayOf(
            "5000000\t5000\t5\t5000000\t0\tstr5",
            "1000000\t1000\t1\t1000000\t0\tstr1",
            "3000000\t3000\t3\t3000000\t1\tstr3",
            "4000000\t4000\t4\t4000000\t1\tstr4",
            "2000000\t2000\t2\t2000000\t0\tstr2"
        ).joinToString("\n")

        withTempFile("df", ".tsv") { path ->
            path.write(blob)

            val spec = DataFrameSpec().ints("i").shorts("s").bytes("b").longs("l").booleans("f").strings("str")
            assertEquals(
                testDataFrame.omit("e"),
                DataFrameMappers.TSV.load(path, spec, header = false)
            )
        }
    }

    @Test
    fun testLoadNonNumericQuoting() {
        // backward compatibility with previous TSV format: quote non-numeric
        val blob = arrayOf(
            "\"i\"\t\"s\"\t\"b\"\t\"l\"\t\"f\"\t\"str\"",
            "5000000\t5000\t5\t5000000\t0\t\"str5\"",
            "1000000\t1000\t1\t1000000\t0\t\"str1\"",
            "3000000\t3000\t3\t3000000\t1\t\"str3\"",
            "4000000\t4000\t4\t4000000\t1\t\"str4\"",
            "2000000\t2000\t2\t2000000\t0\t\"str2\""
        ).joinToString("\n")

        withTempFile("df", ".tsv") { path ->
            path.write(blob)

            val spec = DataFrameSpec().ints("i").shorts("s").bytes("b").longs("l").booleans("f").strings("str")
            assertEquals(testDataFrame.omit("e"), DataFrameMappers.TSV.load(path, spec))
        }
    }

    @Test
    fun testDump() {
        val expected = arrayOf(
            "# Integer;\tShort;\tByte;\tLong;\tBoolean;\tString;\torg.jetbrains.bio.dataframe.TestEnum",
            "i\ts\tb\tl\tf\tstr\te",
            "5000000\t5000\t5\t5000000\t0\tstr5\tVAL5",
            "1000000\t1000\t1\t1000000\t0\tstr1\tVAL1",
            "3000000\t3000\t3\t3000000\t1\tstr3\tVAL3",
            "4000000\t4000\t4\t4000000\t1\tstr4\tVAL4",
            "2000000\t2000\t2\t2000000\t0\tstr2\tVAL2",
            ""
        ).joinToString("\n")

        val buff = StringBuilder()
        DataFrameMappers.TSV.save(buff, testDataFrame, true, true)
        assertEquals(expected, buff.toString())
    }

    @Test
    fun testDumpHead_US() {
        withLocaleUS {
            val dump = testDataFrame.dumpHead(3)
            val expected = arrayOf(
                "# Integer;\tShort;\tByte;\tLong;\tBoolean;\tString;\torg.jetbrains.bio.dataframe.TestEnum",
                "i\ts\tb\tl\tf\tstr\te",
                "5000000\t5000\t5\t5000000\t0\tstr5\tVAL5",
                "1000000\t1000\t1\t1000000\t0\tstr1\tVAL1",
                "3000000\t3000\t3\t3000000\t1\tstr3\tVAL3",
                ""
            ).joinToString("\n")
            assertEquals(expected, dump)
        }
    }
    @Test
    fun testDumpHead_EU() {
        withLocaleEU {
            val dump = testDataFrame.dumpHead(3)
            val expected = arrayOf(
                "# Integer;\tShort;\tByte;\tLong;\tBoolean;\tString;\torg.jetbrains.bio.dataframe.TestEnum",
                "i\ts\tb\tl\tf\tstr\te",
                "5000000\t5000\t5\t5000000\t0\tstr5\tVAL5",
                "1000000\t1000\t1\t1000000\t0\tstr1\tVAL1",
                "3000000\t3000\t3\t3000000\t1\tstr3\tVAL3",
                ""
            ).joinToString("\n")
            assertEquals(expected, dump)
        }
    }

    @Test
    fun testDumpChangedFormat() {
        val expected = """
            |# Integer;	Short;	Byte;	Long;	Boolean;	String;	org.jetbrains.bio.dataframe.TestEnum;	Double;	Float
            |i	s	b	l	f	str	e	d	fl
            |i:5000000	s:5000	b:5	l:5000000	f:false	str5	VAL5	1666666.6667	1666666.63
            |i:1000000	s:1000	b:1	l:1000000	f:false	str1	VAL1	333333.3333	333333.34
            |i:3000000	s:3000	b:3	l:3000000	f:true	str3	VAL3	1000000.0000	1000000.00
            |i:4000000	s:4000	b:4	l:4000000	f:true	str4	VAL4	1333333.3333	1333333.38
            |i:2000000	s:2000	b:2	l:2000000	f:false	str2	VAL2	666666.6667	666666.69
            |
        """.trimMargin()

        val buff = StringBuilder()
        val tdf = testDataFrame
        val newDf = tdf
            .with("i", tdf.sliceAsInt("i")) { "i:${it}" }
            .with("s", tdf.sliceAsShort("s")) { "s:${it}" }
            .with("b", tdf.sliceAsByte("b")) { "b:${it}" }
            .with("l", tdf.sliceAsLong("l")) { "l:${it}" }
            .with("f", tdf.sliceAsBool("f")) { "f:${it}" }
            .with("str", tdf.sliceAsObj("str")) { "str:${it}" }
            .with("e", TestEnum::class.java, tdf.sliceAsObj("e")) { "e:${it}" }
            .with("d", tdf.sliceAsInt("i").map { it / 3.0 }.toDoubleArray()) { Formatter(Locale.US).format("%.4f", it).toString() }
            .with("fl", tdf.sliceAsInt("i").map { it / 3.0f }.toFloatArray()) { Formatter(Locale.US).format("%.2f", it).toString() }
        DataFrameMappers.TSV.save(buff, newDf, true, true)
        Assert.assertEquals(
            expected,
            buff.toString()
        )
    }

    @Test
    fun testDumpChangedFormatAndRenameCols() {
        val expected = """
            |# Integer;	Short;	Byte;	Long;	Boolean;	String;	org.jetbrains.bio.dataframe.TestEnum;	Double;	Float
            |[i]	[s]	[b]	[l]	[f]	[str]	[e]	[d]	[fl]
            |i:5000000	s:5000	b:5	l:5000000	f:false	str5	VAL5	1666666.6667	1666666.63
            |i:1000000	s:1000	b:1	l:1000000	f:false	str1	VAL1	333333.3333	333333.34
            |i:3000000	s:3000	b:3	l:3000000	f:true	str3	VAL3	1000000.0000	1000000.00
            |i:4000000	s:4000	b:4	l:4000000	f:true	str4	VAL4	1333333.3333	1333333.38
            |i:2000000	s:2000	b:2	l:2000000	f:false	str2	VAL2	666666.6667	666666.69
            |
        """.trimMargin()

        val buff = StringBuilder()
        val tdf = testDataFrame
        val newDf = tdf
            .with("i", tdf.sliceAsInt("i")) { "i:${it}" }
            .with("s", tdf.sliceAsShort("s")) { "s:${it}" }
            .with("b", tdf.sliceAsByte("b")) { "b:${it}" }
            .with("l", tdf.sliceAsLong("l")) { "l:${it}" }
            .with("f", tdf.sliceAsBool("f")) { "f:${it}" }
            .with("str", tdf.sliceAsObj("str")) { "str:${it}" }
            .with("e", TestEnum::class.java, tdf.sliceAsObj("e")) { "e:${it}" }
            .with("d", tdf.sliceAsInt("i").map { it / 3.0 }.toDoubleArray()) { Formatter(Locale.US).format("%.4f", it).toString() }
            .with("fl", tdf.sliceAsInt("i").map { it / 3.0f }.toFloatArray()) { Formatter(Locale.US).format("%.2f", it).toString() }

        val renamedDf = DataFrame(
            rowsNumber = newDf.rowsNumber,
            columns = newDf.columns.map { it.rename("[${it.label}]") }
        )
        DataFrameMappers.TSV.save(buff, renamedDf, true, true);
        Assert.assertEquals(
            expected,
            buff.toString()
        )
    }

    @Test
    fun testDumpChangedFormatAndResizeDown() {
        val expected = """
            |# Integer;	Short;	Byte;	Long;	Boolean;	String;	org.jetbrains.bio.dataframe.TestEnum;	Double;	Float
            |i	s	b	l	f	str	e	d	fl
            |i:5000000	s:5000	b:5	l:5000000	f:false	str5	VAL5	1666666.6667	1666666.63
            |i:1000000	s:1000	b:1	l:1000000	f:false	str1	VAL1	333333.3333	333333.34
            |i:3000000	s:3000	b:3	l:3000000	f:true	str3	VAL3	1000000.0000	1000000.00
            |
        """.trimMargin()

        val buff = StringBuilder()
        val tdf = testDataFrame
        val newDf = tdf
            .with("i", tdf.sliceAsInt("i")) { "i:${it}" }
            .with("s", tdf.sliceAsShort("s")) { "s:${it}" }
            .with("b", tdf.sliceAsByte("b")) { "b:${it}" }
            .with("l", tdf.sliceAsLong("l")) { "l:${it}" }
            .with("f", tdf.sliceAsBool("f")) { "f:${it}" }
            .with("str", tdf.sliceAsObj("str")) { "str:${it}" }
            .with("e", TestEnum::class.java, tdf.sliceAsObj("e")) { "e:${it}" }
            .with("d", tdf.sliceAsInt("i").map { it / 3.0 }.toDoubleArray()) { Formatter(Locale.US).format("%.4f", it).toString() }
            .with("fl", tdf.sliceAsInt("i").map { it / 3.0f }.toFloatArray()) { Formatter(Locale.US).format("%.2f", it).toString() }

        val resizedDf = newDf.resize(newDf.rowsNumber - 2)
        DataFrameMappers.TSV.save(buff, resizedDf, true, true)
        Assert.assertEquals(
            expected,
            buff.toString()
        )
    }


}

class NpzMapperTest {
    private val specNoEnum = DataFrameSpec()
        .ints("i").shorts("s").bytes("b").longs("l").booleans("f").strings("str")

    @Test
    fun testGuess() {
        withTempFile("df", ".npz") { path ->
            val df = testDataFrame.omit("e")
            DataFrameMappers.NPZ.save(path, df)
            assertEquals(specNoEnum, DataFrameMappers.NPZ.guess(path))
        }
    }

    @Test
    fun testSaveLoad() {
        withTempFile("df", ".npz") { path ->
            val df = testDataFrame.omit("e")
            DataFrameMappers.NPZ.save(path, df)
            assertEquals(df, DataFrameMappers.NPZ.load(path, specNoEnum))
        }
    }
}
