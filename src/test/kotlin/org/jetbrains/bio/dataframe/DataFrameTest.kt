package org.jetbrains.bio.dataframe

import org.apache.commons.math3.util.Precision
import org.junit.Assert.assertArrayEquals
import org.junit.Test
import java.util.function.IntPredicate
import kotlin.test.assertEquals

class DataFrameTest {
    val TRUE = RowPredicateFactory { IntPredicate { true } }

    @Test fun testWithAdd() {
        val x = doubleArrayOf(42.0)
        val y = doubleArrayOf(24.0)
        val df = DataFrame().with("x", x).with("y", y)
        assertEquals(2, df.columnsNumber)
        assertArrayEquals(y, df.sliceAsDouble("y"), Precision.EPSILON)
    }

    @Test fun testWithReplace() {
        val x = doubleArrayOf(42.0)
        val y = doubleArrayOf(24.0)
        val df = DataFrame().with("x", x).with("x", y)
        assertEquals(1, df.columnsNumber)
        assertArrayEquals(y, df.sliceAsDouble("x"), Precision.EPSILON)
    }

    @Test fun testSlice() {
        val df = DataFrame()
            .with("x", doubleArrayOf(1.0, 2.0, 3.0, 4.0))
            .with("y", doubleArrayOf(5.0, 6.0, 7.0, 8.0))

        assertEquals(4, df.sliceAsDouble("x").size)
        assertEquals(1.0, df.sliceAsDouble("x")[0])
    }

    @Test fun testGet() {
        val df = testDataFrame
        assertEquals(1, df.getAsInt(0, "x"))
        assertEquals(1.0, df.getAsDouble(0, "y"))
        assertEquals("a", df.getAsObj<Any>(0, "z"))
    }

    @Test fun testRowAsDouble() {
        val df = DataFrame()
                .with("x", doubleArrayOf(1.0, 2.0, 3.0, 4.0))
                .with("y", doubleArrayOf(5.0, 6.0, 7.0, 8.0))

        assertEquals(2, df.rowAsDouble(0).size)
        assertEquals(1.0, df.rowAsDouble(0)[0])
        assertEquals(2.0, df.rowAsDouble(1)[0])
        assertEquals(8.0, df.rowAsDouble(3)[1])
    }

    @Test fun testFilter() {
        val filtered = testDataFrame.filter(byInt("x") { it % 2 == 0 })

        assertArrayEquals(testDataFrame.labels, filtered.labels)
        assertEquals(1, filtered.rowsNumber)
        assertArrayEquals(intArrayOf(2), filtered.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(2.0), filtered.sliceAsDouble("y"), Precision.EPSILON)
        assertArrayEquals(arrayOf("b"), filtered.sliceAsObj("z"))
    }

    @Test fun ilocRangeStep1() {
        assertEquals(testDataFrame.filter(TRUE, 0, 2),
                     testDataFrame.iloc[0 until 2])

        assertEquals(testDataFrame.filter(TRUE, 0, 2),
                     testDataFrame.iloc[0 .. 1])

        assertEquals(testDataFrame.filter(TRUE, 1, 3),
                     testDataFrame.iloc[1 until 3])
    }

    @Test fun ilocRangeStep2() {
        assertEquals(testDataFrame.filter(RowPredicateFactory { IntPredicate { row -> row % 2 == 0 } }, 0, 3),
                     testDataFrame.iloc[0..2 step 2])

        assertEquals(testDataFrame.filter(RowPredicateFactory { IntPredicate { row -> row % 2 == 0 } }, 0, 2),
                     testDataFrame.iloc[0..1 step 2])

        assertEquals(testDataFrame.filter(TRUE, 0, 3),
                     testDataFrame.iloc[2 downTo 0])

        assertEquals(testDataFrame.filter(RowPredicateFactory { IntPredicate { row -> row % 2 == 0 } }, 0, 3),
                     testDataFrame.iloc[2 downTo 0 step 2])

        assertEquals(testDataFrame.filter(RowPredicateFactory { IntPredicate { row -> row % 2 == 0 } }, 0, 3),
                     testDataFrame.iloc[0..3 step 2])

    }

    @Test(expected = IndexOutOfBoundsException::class)
    fun ilocOOE1() {
        testDataFrame.iloc[0..5]
    }

    @Test(expected = IndexOutOfBoundsException::class)
    fun ilocOOE2() {
        testDataFrame.iloc[2 downTo -5]
    }

    @Test fun testAsMask() {
        val df = testDataFrame
        assertEquals("3@{1}", df.test(byInt("x") { it % 2 == 0 }).toString())
        assertEquals("3@{0, 2}", df.test(byInt("x") { it % 2 == 1 }).toString())
    }

    @Test fun testColumnBindNoExclude() {
        val df1 = testDataFrame
        val df2 = testDataFrame

        val df = DataFrame.columnBind(df1, df2)
        assertEquals(6, df.columnsNumber)
        assertEquals(3, df.rowsNumber)
        val labels = arrayOf("x1", "y1", "z1", "x2", "y2", "z2")
        assertArrayEquals(labels, df.labels)
    }

    @Test fun testColumnBindWithExclude() {
        val df1 = testDataFrame
        val df2 = testDataFrame

        val df = DataFrame.columnBind(setOf("x"), df1, df2)
        assertEquals(4, df.columnsNumber)
        assertEquals(3, df.rowsNumber)
        val labels = arrayOf("y1", "z1", "y2", "z2")
        assertArrayEquals(labels, df.labels)
    }

    @Test fun testColumnBindThreesome() {
        val df1 = testDataFrame
        val df2 = DataFrame()
            .with("x", df1.sliceAsInt("x"))
            .with("y", doubleArrayOf(2.0, 3.0, 4.0))
        val df3 = DataFrame()
            .with("x", df1.sliceAsInt("x"))
            .with("z", doubleArrayOf(3.0, 4.0, 5.0))
            .with("u", doubleArrayOf(3.0, 4.0, 5.0))

        val df = DataFrame.columnBind(setOf("x"), df1, df2, df3)
        assertEquals(5, df.columnsNumber)
        assertEquals(3, df.rowsNumber)
        val labels = arrayOf("y1", "z1", "y2", "z2", "u")
        assertArrayEquals(labels, df.labels)
    }

    @Test fun testRowBind() {
        val df1 = testDataFrame
        val df2 = DataFrame()
            .with("x", intArrayOf(4))
            .with("y", doubleArrayOf(4.0))
            .with("z", arrayOf("d"))

        val df = DataFrame.rowBind(df1, df2)
        assertEquals(df1.columnsNumber, df.columnsNumber)
        assertEquals((df1.rowsNumber + df2.rowsNumber), df.rowsNumber)
        assertArrayEquals(intArrayOf(1, 2, 3, 4), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 4.0),
                          df.sliceAsDouble("y"), Precision.EPSILON)
        assertArrayEquals(arrayOf("a", "b", "c", "d"),
                          df.sliceAsObj<String>("z"))

    }

    @Test fun testResizeSmaller() {
        val df = testDataFrame
        val dfSmaller = df.resize(2)
        assertEquals(df.columnsNumber, dfSmaller.columnsNumber)
        assertEquals(2, dfSmaller.rowsNumber)
        assertArrayEquals(intArrayOf(1, 2), dfSmaller.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0),
                          dfSmaller.sliceAsDouble("y"), Precision.EPSILON)
        assertArrayEquals(arrayOf("a", "b"), dfSmaller.sliceAsObj<String>("z"))
    }

    @Test fun testResizeLarger() {
        val df = testDataFrame
        val dfLarger = df.resize(4)
        assertEquals(df.columnsNumber, dfLarger.columnsNumber)
        assertEquals(4, dfLarger.rowsNumber)
        assertArrayEquals(intArrayOf(1, 2, 3, 0), dfLarger.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 0.0),
                          dfLarger.sliceAsDouble("y"), Precision.EPSILON)
        assertArrayEquals(arrayOf("a", "b", "c", null), dfLarger.sliceAsObj<String>("z"))
    }

    @Test fun testReorderMovesLinkedValues() {
        val df = DataFrame()
            .with("x", intArrayOf(5, 4, 3))
            .with("y", doubleArrayOf(3.0, 4.0, 5.0))

        val orderedByX = df.reorder("x")
        assertEquals(df.columnsNumber.toLong(), orderedByX.columnsNumber.toLong())
        assertEquals(df.rowsNumber.toLong(), orderedByX.rowsNumber.toLong())
        assertArrayEquals(intArrayOf(3, 4, 5), orderedByX.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(5.0, 4.0, 3.0), orderedByX.sliceAsDouble("y"),
                          Precision.EPSILON)
    }

    @Test fun testReorder() {
        val bitSet = BitterSet(5)
        bitSet.set(2)
        bitSet.set(3)
        val df = DataFrame()
            .with("i", intArrayOf(5, 1, 3, 4, 2))
            .with("d", doubleArrayOf(5.0, 1.0, 3.0, 4.0, 2.0))
            .with("s", shortArrayOf(5, 1, 3, 4, 2))
            .with("b", byteArrayOf(5, 1, 3, 4, 2))
            .with("l", longArrayOf(5, 1, 3, 4, 2))
            .with("F34", bitSet)
            .with("e", TestEnum::class.java,
                  arrayOf(TestEnum.VAL5, TestEnum.VAL1, TestEnum.VAL3,
                          TestEnum.VAL4, TestEnum.VAL2))
                .with("t", arrayOf("boo", "baz", "bar", "yada", "do"))

        assertArrayEquals(intArrayOf(1, 2, 3, 4, 5), df.reorder("i").sliceAsInt("i"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 4.0, 5.0),
                          df.reorder("d").sliceAsDouble("d"), Precision.EPSILON)
        assertArrayEquals(shortArrayOf(1, 2, 3, 4, 5), df.reorder("s").sliceAsShort("s"))
        assertArrayEquals(byteArrayOf(1, 2, 3, 4, 5), df.reorder("b").sliceAsByte("b"))
        assertArrayEquals(longArrayOf(1, 2, 3, 4, 5), df.reorder("l").sliceAsLong("l"))
        assertArrayEquals(arrayOf(TestEnum.VAL1, TestEnum.VAL2, TestEnum.VAL3,
                TestEnum.VAL4, TestEnum.VAL5),
                          df.reorder("e").sliceAsObj<Any>("e"))
        assertEquals("5@{2, 3}", df.reorder("l").sliceAsBool("F34").toString())

        assertArrayEquals(intArrayOf(1, 2, 3, 4, 5),
                          df.reorder("i", reverse = false).sliceAsInt("i"))

        assertArrayEquals(intArrayOf(5, 4, 3, 2, 1),
                          df.reorder("i", reverse = true).sliceAsInt("i"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 4.0, 5.0),
                          df.reorder("d").sliceAsDouble("d"), Precision.EPSILON)
        assertArrayEquals(shortArrayOf(5, 4, 3, 2, 1),
                          df.reorder("s", reverse = true).sliceAsShort("s"))
        assertArrayEquals(byteArrayOf(5, 4, 3, 2, 1),
                          df.reorder("b", reverse = true).sliceAsByte("b"))
        assertArrayEquals(longArrayOf(5, 4, 3, 2, 1),
                          df.reorder("l", reverse = true).sliceAsLong("l"))
        assertArrayEquals(arrayOf(TestEnum.VAL5, TestEnum.VAL4, TestEnum.VAL3,
                TestEnum.VAL2, TestEnum.VAL1),
                          df.reorder("e", reverse = true).sliceAsObj<Any>("e"))
        assertEquals("5@{1, 2}", df.reorder("l", reverse = true).sliceAsBool("F34").toString())

        assertArrayEquals(arrayOf("bar", "baz", "boo", "do", "yada"),
                          df.reorder("t", reverse = false).sliceAsObj<Any>("t"))
        assertArrayEquals(arrayOf("yada", "do", "boo", "baz", "bar"),
                          df.reorder("t", reverse = true).sliceAsObj<Any>("t"))
    }

    @Test fun testReorderByFlag() {
        val bitSet = BitterSet(5)
        bitSet.set(2)
        bitSet.set(3)
        val df = DataFrame()
                .with("i", intArrayOf(5, 1, 3, 4, 2))
                .with("F34", bitSet)

        assertArrayEquals(intArrayOf(5, 1, 2, 3, 4),
                          df.reorder("F34").sliceAsInt("i"))
        assertEquals("5@{3, 4}",
                     df.reorder("F34").sliceAsBool("F34").toString())

        assertArrayEquals(intArrayOf(5, 1, 2, 3, 4),
                          df.reorder("F34", reverse = false).sliceAsInt("i"))
        assertEquals("5@{3, 4}", df.reorder("F34", reverse = false).sliceAsBool("F34").toString())

        assertArrayEquals(intArrayOf(3, 4, 5, 1, 2),
                          df.reorder("F34", reverse = true).sliceAsInt("i"))
        assertEquals("5@{0, 1}", df.reorder("F34", reverse = true).sliceAsBool("F34").toString())
    }

    @Test fun testMergeInner() {
        val df1 = testDataFrame
        val df2 = DataFrame()
                .with("x", intArrayOf(3, 4, 5))
                .with("y", doubleArrayOf(3.0, 4.0, 5.0))

        val df = DataFrame.mergeInner("x", df1, df2)
        assertEquals(4, df.columnsNumber)
        assertEquals(1, df.rowsNumber)

        assertArrayEquals(intArrayOf(3), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y1"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y2"), Precision.EPSILON)
        assertArrayEquals(arrayOf("c"), df.sliceAsObj<String>("z"))
    }

    @Test fun testMergeInnerTIntHasSetRetainSideEffect() {
        val df1 = testDataFrame
        val df2 = DataFrame()
                .with("x", intArrayOf(5, 4, 3))
                .with("y", doubleArrayOf(5.0, 4.0, 3.0))

        val df = DataFrame.mergeInner("x", df1, df2)
        assertEquals(4, df.columnsNumber)
        assertEquals(1, df.rowsNumber)

        assertArrayEquals(intArrayOf(3), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y1"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y2"), Precision.EPSILON)
        assertArrayEquals(arrayOf("c"), df.sliceAsObj<String>("z"))
    }

    @Test fun testMergeInnerThreesome() {
        val df1 = testDataFrame
        val df2 = DataFrame()
                .with("x", intArrayOf(2, 3, 4, 5))
                .with("y", doubleArrayOf(2.0, 3.0, 4.0, 5.0))
        val df3 = DataFrame()
                .with("x", intArrayOf(3, 4, 5))
                .with("y", doubleArrayOf(3.0, 4.0, 5.0))

        val df = DataFrame.mergeInner("x", df1, df2, df3)
        assertEquals(5, df.columnsNumber)
        assertEquals(1, df.rowsNumber)

        assertArrayEquals(intArrayOf(3), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y1"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y2"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(3.0), df.sliceAsDouble("y3"), Precision.EPSILON)
        assertArrayEquals(arrayOf("c"), df.sliceAsObj<String>("z"))
    }

    @Test fun testMergeOuter() {
        val df1 = testDataFrame
        val df2 = DataFrame()
                .with("x", intArrayOf(5, 4, 3))
                .with("y", doubleArrayOf(5.0, 4.0, -3.0))

        val df = DataFrame.mergeOuter("x", df1, df2)
        assertEquals(4, df.columnsNumber)
        assertEquals(5, df.rowsNumber)

        assertArrayEquals(intArrayOf(1, 2, 3, 4, 5), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 0.0, 0.0),
                          df.sliceAsDouble("y1"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(0.0, 0.0, -3.0, 4.0, 5.0),
                          df.sliceAsDouble("y2"), Precision.EPSILON)
        assertArrayEquals(arrayOf("a", "b", "c", null, null),
                          df.sliceAsObj<String>("z"))
    }

    @Test fun testMergeOuterThreesome() {
        val df1 = testDataFrame
        val df2 = DataFrame()
                .with("x", intArrayOf(2, 3, 4, 5))
                .with("y", doubleArrayOf(-2.0, -3.0, -4.0, -5.0))
        val df3 = DataFrame()
                .with("x", intArrayOf(3, 4, 5))
                .with("y", doubleArrayOf(0.3, 0.4, 0.5))

        val df = DataFrame.mergeOuter("x", df1, df2, df3)
        assertEquals(5, df.columnsNumber)
        assertEquals(5, df.rowsNumber)

        assertArrayEquals(intArrayOf(1, 2, 3, 4, 5), df.sliceAsInt("x"))
        assertArrayEquals(doubleArrayOf(1.0, 2.0, 3.0, 0.0, 0.0),
                          df.sliceAsDouble("y1"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(0.0, -2.0, -3.0, -4.0, -5.0),
                          df.sliceAsDouble("y2"), Precision.EPSILON)
        assertArrayEquals(doubleArrayOf(0.0, 0.0, 0.3, 0.4, 0.5),
                          df.sliceAsDouble("y3"), Precision.EPSILON)
        assertArrayEquals(arrayOf("a", "b", "c", null, null),
                          df.sliceAsObj<String>("z"))
    }

    private val testDataFrame: DataFrame
        get() = DataFrame()
            .with("x", intArrayOf(1, 2, 3))
            .with("y", doubleArrayOf(1.0, 2.0, 3.0))
            .with("z", arrayOf("a", "b", "c"))
}
