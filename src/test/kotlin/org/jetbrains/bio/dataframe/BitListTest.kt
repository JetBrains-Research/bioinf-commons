package org.jetbrains.bio.dataframe

import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.Range
import org.junit.Test
import org.junit.runner.RunWith
import org.junit.runners.Parameterized
import org.junit.runners.Parameterized.Parameters
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class BitListTest {
    @Test
    fun setIOOB() {
        Tests.assertThrowsWithMessage(
            IndexOutOfBoundsException::class.java,
            "bitIndex >= size: 10 >= 10",
        ) {
            BitList(10).set(10)
        }
    }

    @Test
    fun setIOOB2() {
        Tests.assertThrowsWithMessage(
            IndexOutOfBoundsException::class.java,
            "toIndex > size: 0 > 10",
        ) {
            BitList(10).set(0, 11, true)
        }
    }

    @Test
    fun clearIOOB() {
        Tests.assertThrowsWithMessage(
            IndexOutOfBoundsException::class.java,
            "bitIndex >= size: 10 >= 10",
        ) {
            BitList(10).clear(10)
        }

    }

    @Test
    fun clearIOOB2() {
        Tests.assertThrowsWithMessage(
            IndexOutOfBoundsException::class.java,
            "toIndex > size: 0 > 10",
        ) {
            BitList(10).clear(0, 11)
        }
    }


    @Test
    fun copy() {
        val b = BitList(10)
        val copy = b.copy()
        assertEquals(b, copy)
        assertEquals(b.hashCode(), copy.hashCode())
        copy.set(4)
        assertFalse(b[4])
    }

    @Test
    fun constructor() {
        val b = BitList(10) { it < 4 }
        for (i in 0 until 4) {
            assertTrue(b[i])
        }
    }

    @Test
    fun empty() {
        assertEquals(0, BitList(0).size())
    }

    @Test
    fun equalsEmpty() {
        assertEquals(BitList(0), BitList(0))
    }

    @Test
    fun hashCodeConsistency() {
        val b1 = BitList(6) { it in setOf(2, 4, 5) }
        val b2 = BitList(6) { it in setOf(2, 4, 5) }
        assertEquals(b1, b2)
        assertEquals(b1.hashCode(), b2.hashCode())
    }

    @Test
    fun plusEmpty() {
        val bs = BitList(10)
        assertEquals(bs, bs + BitList(0))
    }

    @Test
    fun plusNonEmpty() {
        val bs1 = BitList(10) { it in intArrayOf(1, 2) }
        val bs2 = BitList(5) { it == 1 }

        val bs = bs1 + bs2
        assertEquals(2, bs1.cardinality())
        assertEquals(10, bs1.size())
        assertEquals(1, bs2.cardinality())
        assertEquals(5, bs2.size())

        assertEquals(bs1.cardinality() + bs2.cardinality(), bs.cardinality())
        assertEquals(bs1.size() + bs2.size(), bs.size())
        assertTrue(bs[1] && bs[2])
        assertTrue(bs[11])
    }

    @Test
    fun iteratorEmpty() {
        val bs = BitList(10)
        var count = 0
        bs.iterator().forEach { count++ }
        assertEquals(0, count)
    }

    @Test
    fun iteratorNotEmpty() {
        val bs = BitList(6) { it in setOf(2, 4, 5) }
        val setBits = mutableListOf<Int>()
        bs.iterator().forEach { setBits.add(it) }
        assertEquals(3, setBits.size)
        assertEquals(listOf(2, 4, 5), setBits)
    }

    @Test
    fun toBitterSet() {
        assertEquals(booleanArrayOf(false, false, true, false, true, true).toBitList(),
            BitList(6) { it in setOf(2, 4, 5) })
        assertEquals(
            booleanArrayOf(false, false, false, false).toBitList(), BitList(4)
        )
        assertEquals(booleanArrayOf(true, true, true, true, true).toBitList(), BitList(5) { true })
    }

    @Test
    fun getOutOfUniverse() {
        val bs = BitList(6) { it in setOf(2, 4, 5) }
        require(bs.get(5))

        assertFalse(bs.get(6))
        assertFalse(bs.get(7))
        assertFalse(bs.get(1000))
    }
}

@RunWith(Parameterized::class)
class AggregateTest(private val expected: List<Range>, private val bits: IntArray) {
    @Test
    fun aggregate() {
        val actual = bits.toBitList()
        assertEquals(expected, actual.aggregate())
    }

    companion object {
        @JvmStatic
        @Parameters
        fun `data`(): Collection<Array<out Any>> = listOf(
            arrayOf(emptyList<Range>(), intArrayOf()),
            arrayOf(listOf(Range(1, 4)), intArrayOf(1, 2, 3)),
            arrayOf(listOf(Range(1, 4), Range(5, 7)), intArrayOf(1, 2, 3, 5, 6)),
            arrayOf(
                listOf(Range(1, 4), Range(5, 7), Range(10, 11), Range(15, 16)), intArrayOf(1, 2, 3, 5, 6, 10, 15)
            ),
            arrayOf(emptyList<Range>(), intArrayOf()),
            arrayOf(listOf(Range(1, 3), Range(4, 6)), intArrayOf(1, 2, 4, 5))
        )
    }
}