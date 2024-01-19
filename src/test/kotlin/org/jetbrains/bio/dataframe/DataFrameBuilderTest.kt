package org.jetbrains.bio.dataframe

import org.jetbrains.bio.Tests
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class DataFrameBuilderTest {
    @Test
    fun testBuilder() {
        val builder = DataFrameSpec()
            .ints("i")
            .bytes("b")
            .doubles("d")
            .strings("s")
            .booleans("f")
            .enums("e", TestEnum::class.java)
            .builder()
        builder.add(1000, 1.toByte(), 1.0, "a", true, TestEnum.VAL1)
        builder.add(2000, 2.toByte(), 2.0, "b", false, TestEnum.VAL2)
        builder.add(3000, 3.toByte(), 3.0, "c", true, TestEnum.VAL1)

        val bitList = BitList(3)
        bitList.set(0)
        bitList.set(2)

        val expectedDf = DataFrame()
            .with("i", intArrayOf(1000, 2000, 3000))
            .with("b", byteArrayOf(1, 2, 3))
            .with("d", doubleArrayOf(1.0, 2.0, 3.0))
            .with("s", arrayOf("a", "b", "c"))
            .with("f", bitList)
            .with(
                "e", TestEnum::class.java,
                arrayOf(TestEnum.VAL1, TestEnum.VAL2, TestEnum.VAL1)
            )


        val df = builder.build()
        assertTrue(df["i"] is IntColumn)
        assertTrue(df["b"] is ByteColumn)
        assertTrue(df["d"] is DoubleColumn)
        assertTrue(df["s"] is StringColumn)
        assertTrue(df["f"] is BooleanColumn)
        assertTrue(df["e"] is EnumColumn<*>)

        assertEquals(expectedDf, df)
    }

    @Test
    fun testCheckWrongIntRow() {
        val builder = DataFrameSpec().ints("v").builder()

        Tests.assertThrowsWithMessage(
            "Wrong type: 0-th arg (column: v) value 1 is expected to be of type Integer",
            IllegalArgumentException::class.java,
        ) {
            builder.add(1.toLong())
            builder.build()
        }
    }

    @Test
    fun testCheckWrongLongRow() {
        val builder = DataFrameSpec().longs("v").builder()

        Tests.assertThrowsWithMessage(
            "Wrong type: 0-th arg (column: v) value 1 is expected to be of type Long",
            IllegalArgumentException::class.java,
        ) {
            builder.add(1)
            builder.build()
        }
    }

    @Test
    fun testCheckWrongByteRow() {
        val builder = DataFrameSpec().bytes("v").builder()

        Tests.assertThrowsWithMessage(
            "Wrong type: 0-th arg (column: v) value 1 is expected to be of type Byte",
            IllegalArgumentException::class.java,
        ) {
            builder.add(1)
            builder.build()
        }
    }

    @Test
    fun testCheckWrongDoubleRow() {
        val builder = DataFrameSpec().doubles("v").builder()

        Tests.assertThrowsWithMessage(
            "is expected to be of type Double",
            IllegalArgumentException::class.java, partialMessageMatch = true
        ) {
            builder.add(1)
            builder.build()
        }

    }

    @Test
    fun testCheckWrongStringRow() {
        val builder = DataFrameSpec().strings("v").builder()

        Tests.assertThrowsWithMessage(
            "Wrong type: 0-th arg (column: v) value 1 is expected to be of type String",
            IllegalArgumentException::class.java,
        ) {
            builder.add(1)
            builder.build()
        }

    }
}