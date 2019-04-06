package org.jetbrains.bio.dataframe

import org.junit.Test
import org.junit.runner.RunWith
import org.junit.runners.Parameterized
import org.junit.runners.Parameterized.Parameters
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class BitterSetTest {
    @Test fun copy() {
        val b = BitterSet(10)
        val copy = b.copy()
        assertEquals(b, copy)
        assertEquals(b.hashCode(), copy.hashCode())
        copy.set(4)
        assertFalse(b[4])
    }

    @Test fun constructor() {
        val b = BitterSet(10) { it < 4 }
        for (i in 0 until 4) {
            assertTrue(b[i])
        }
    }

    @Test fun empty() {
        assertEquals(0, BitterSet(0).size())
    }

    @Test fun equalsEmpty() {
        assertEquals(BitterSet(0), BitterSet(0))
    }

    @Test fun hashCodeConsistency() {
        val b1 = BitterSet(6) { it in setOf(2, 4, 5) }
        val b2 = BitterSet(6) { it in setOf(2, 4, 5) }
        assertEquals(b1, b2)
        assertEquals(b1.hashCode(), b2.hashCode())
    }

    @Test fun plusEmpty() {
        val bs = BitterSet(10)
        assertEquals(bs, bs + BitterSet(0))
    }

    @Test fun plusNonEmpty() {
        val bs1 = BitterSet(10) { it in intArrayOf(1, 2) }
        val bs2 = BitterSet(5) { it == 1 }

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

    @Test fun iteratorEmpty() {
        val bs = BitterSet(10)
        var count = 0
        bs.iterator().forEach { count++ }
        assertEquals(0, count)
    }

    @Test fun iteratorNotEmpty() {
        val bs = BitterSet(6) { it in setOf(2, 4, 5) }
        val setBits = mutableListOf<Int>()
        bs.iterator().forEach { setBits.add(it) }
        assertEquals(3, setBits.size)
        assertEquals(listOf(2,4, 5), setBits)
    }

    @Test fun toBitterSet() {
        assertEquals(
                booleanArrayOf(false, false, true, false, true, true).toBitterSet(),
                BitterSet(6) { it in setOf(2, 4, 5) }
        )
        assertEquals(
                booleanArrayOf(false, false, false, false).toBitterSet(),
                BitterSet(4)
        )
        assertEquals(
                booleanArrayOf(true, true, true, true, true).toBitterSet(),
                BitterSet(5) { true }
        )
    }

    @Test fun getOutOfUniverse() {
        val bs = BitterSet(6) { it in setOf(2, 4, 5) }
        require(bs.get(5))
        
        assertFalse(bs.get(6))
        assertFalse(bs.get(7))
        assertFalse(bs.get(1000))
    }
}

@RunWith(Parameterized::class)
class AggregateTest(private val expected: List<BitRange>, private val bits: IntArray, private val gap: Int) {
    @Test fun aggregate() {
        val actual = BitterSet((bits.max() ?: 0) + 1) { it in bits }
        assertEquals(expected, actual.aggregate(gap))
    }

    companion object {
        @JvmStatic
        @Parameters
        fun `data`(): Collection<Array<out Any>> = listOf(
                arrayOf(emptyList<BitRange>(), intArrayOf(), 0),
                arrayOf(listOf(BitRange(1, 4)), intArrayOf(1, 2, 3), 0),
                arrayOf(listOf(BitRange(1, 4), BitRange(5, 7)), intArrayOf(1, 2, 3, 5, 6), 0),
                arrayOf(listOf(BitRange(1, 4), BitRange(5, 7), BitRange(10, 11), BitRange(15, 16)),
                        intArrayOf(1, 2, 3, 5, 6, 10, 15), 0),
                arrayOf(emptyList<BitRange>(), intArrayOf(), 10),
                arrayOf(listOf(BitRange(1, 18), BitRange(30, 41)), intArrayOf(1, 2, 3, 4, 10, 15, 17, 30, 35, 40), 5),
                arrayOf(listOf(BitRange(1, 3), BitRange(4, 6)), intArrayOf(1, 2, 4, 5), 0),
                arrayOf(listOf(BitRange(1, 6)), intArrayOf(1, 2, 4, 5), 1))
    }
}