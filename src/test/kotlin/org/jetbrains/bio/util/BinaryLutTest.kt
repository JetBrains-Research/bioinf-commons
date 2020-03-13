package org.jetbrains.bio.util

import org.junit.Test
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class BinaryLutTest {
    @Test fun random() {
        val data = RANDOM.ints(512, 0, 1024).toArray()
        data.sort()

        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = RANDOM.nextInt(2048)
            assertSameIndex(data, key, lut)
        }
    }

    @Test fun randomNegative() {
        val data = RANDOM.ints(512, 0, 1024).toArray()
        data.sort()

        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = -RANDOM.nextInt(1024) - 1
            assertSameIndex(data, key, lut)
        }
    }

    @Test fun randomNearestElementDist() {
        val data = RANDOM.ints(512, 0, 256).toArray()
        data.sort()

        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = RANDOM.nextInt(256)
            val (low, high) = lut.nearestElemDist(data, key)
            val insertPoint = data.binarySearch(key).let { if (it < 0) it.inv() else it }
            assertTrue(insertPoint in low-1..high+1,
                    "Expected insert point $insertPoint to be inside range [${low-1},${high+1}]")
        }
    }

    @Test fun randomNearestElementLR() {
        val data = RANDOM.ints(512, 0, 256).toArray()
        data.sort()

        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = RANDOM.nextInt(256)
            val (low, high) = lut.nearestElemLR(data, key)
            val insertPoint = data.binarySearch(key).let { if (it < 0) it.inv() else it }
            assertTrue(insertPoint in low-1..high+1,
                    "Expected insert point $insertPoint to be inside range [${low-1},${high+1}]")
        }
    }

    @Test fun randomElementsWithinDist() {
        val data = RANDOM.ints(512, 0, 256).toArray()
        data.sort()

        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = RANDOM.nextInt(256)
            val left = RANDOM.nextInt(128)
            val right = RANDOM.nextInt(128)
            val (low, high) = lut.elemWithinDist(data, key, left, right)
            assertTrue((low == -1 && high == -1) || (low..high).map { data[it] in key - left..key + right }.all { it },
                    "Not all elements in range are within given distance to the key")
            assertTrue((low == -1 && high == -1) || ((0..low - 1).map { data[it] !in key - left..key + right }.all { it }
                    && (high + 1..data.size - 1).map { data[it] !in key - left..key + right }.all { it }),
                    "Some elements outside of range are within given distance to the key")
            if (low == -1 || high == -1) {
                assertTrue(low == -1 && high == -1, "Only one end of range is -1")
                assertTrue((0..data.size - 1).map { data[it] !in key - left..key + right }.all { it },
                        "Some elements outside of range are within given distance to the key")
            }
        }
    }

    @Test fun powersOfTwo() {
        val data = intArrayOf(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
        val lut = BinaryLut.of(data, 24)
        for (iter in 0..99) {
            val key = RANDOM.nextInt(data.max()!!)
            assertSameIndex(data, key, lut)
        }
    }

    @Test fun specialCase() {
        val data = intArrayOf(4076, 9675)
        val lut = BinaryLut.of(data, 24)
        val actual = lut.binarySearch(data, 368)
        assertEquals(-1, actual, "Expected -1, got $actual, see issue #1033")
    }

    @Test fun nearestElemIdxEmptyArray() {
        val data = intArrayOf()
        val lut = BinaryLut.of(data, 24)
        assertEquals(-1, lut.nearestElemIdx(data, 0))
    }

    @Test fun nearestElemDistEmptyArray() {
        val data = intArrayOf()
        val lut = BinaryLut.of(data, 24)
        assertEquals(-1 to -1, lut.nearestElemDist(data, 0))
    }

    @Test fun nearestElemLREmptyArray() {
        val data = intArrayOf()
        val lut = BinaryLut.of(data, 24)
        assertEquals(-1 to -1, lut.nearestElemLR(data, 0))
    }

    @Test fun elemWithinDistEmptyArray() {
        val data = intArrayOf()
        val lut = BinaryLut.of(data, 24)
        assertEquals(-1 to -1, lut.elemWithinDist(data, 0, 1))
    }

    @Test fun nearestElemIdx() {
        val data = intArrayOf(4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
        val lut = BinaryLut.of(data, 8)
        assertEquals(4, data[lut.nearestElemIdx(data, 0)])
        assertEquals(4, data[lut.nearestElemIdx(data, 3)])
        assertEquals(4, data[lut.nearestElemIdx(data, 4)])
        assertEquals(4, data[lut.nearestElemIdx(data, 5)])
        assertEquals(4, data[lut.nearestElemIdx(data, 6)])
        assertEquals(8, data[lut.nearestElemIdx(data, 7)])
        assertEquals(512, data[lut.nearestElemIdx(data, 511)])
        assertEquals(2048, data[lut.nearestElemIdx(data, 2047)])
        assertEquals(2048, data[lut.nearestElemIdx(data, 2048)])
        assertEquals(2048, data[lut.nearestElemIdx(data, 2049)])
        assertEquals(2048, data[lut.nearestElemIdx(data, 5000)])
    }

    @Test fun nearestElemIdxWithDuplicates() {
        val data = intArrayOf(4, 8, 32, 32, 64, 128, 256, 512, 1024, 2048)
        val lut = BinaryLut.of(data, 8)
        assertEquals(2 to 32, lut.nearestElemIdx(data, 31).let { it to data[it]})
        assertEquals(2 to 32, lut.nearestElemIdx(data, 32).let { it to data[it]})
        assertEquals(3 to 32, lut.nearestElemIdx(data, 33).let { it to data[it]})
    }

    @Test fun nearestElemDistWithDuplicates() {
        val data = intArrayOf(4, 8, 32, 32, 64, 128, 256, 512, 1024, 2048)
        val lut = BinaryLut.of(data, 8)
        assertEquals(0 to 0, lut.nearestElemDist(data, 2))
        assertEquals(0 to 1, lut.nearestElemDist(data, 6))
        assertEquals(1 to 1, lut.nearestElemDist(data, 7))
        assertEquals(1 to 1, lut.nearestElemDist(data, 8))
        assertEquals(1 to 3, lut.nearestElemDist(data, 20))
        assertEquals(2 to 3, lut.nearestElemDist(data, 31))
        assertEquals(2 to 3, lut.nearestElemDist(data, 32))
        assertEquals(2 to 3, lut.nearestElemDist(data, 33))
    }

    @Test fun nearestElemLRWithDuplicates() {
        val data = intArrayOf(4, 8, 32, 32, 64, 128, 256, 512, 1024, 2048)
        val lut = BinaryLut.of(data, 8)
        assertEquals(0 to 0, lut.nearestElemLR(data, 2))
        assertEquals(0 to 1, lut.nearestElemLR(data, 6))
        assertEquals(0 to 1, lut.nearestElemLR(data, 7))
        assertEquals(0 to 3, lut.nearestElemLR(data, 8))
        assertEquals(1 to 3, lut.nearestElemLR(data, 20))
        assertEquals(1 to 3, lut.nearestElemLR(data, 31))
        assertEquals(1 to 4, lut.nearestElemLR(data, 32))
        assertEquals(2 to 4, lut.nearestElemLR(data, 33))
    }

    @Test fun elemWithinDistWithDuplicates() {
        val data = intArrayOf(4, 8, 32, 32, 64, 128, 256, 512, 1024, 2048)
        val lut = BinaryLut.of(data, 8)
        assertEquals(0 to 0, lut.elemWithinDist(data, 2, 2))
        assertEquals(0 to 1, lut.elemWithinDist(data, 6, 2))
        assertEquals(0 to 1, lut.elemWithinDist(data, 7, 8))
        assertEquals(0 to 3, lut.elemWithinDist(data, 8, 4, 24))
        assertEquals(1 to 3, lut.elemWithinDist(data, 20, 12))
        assertEquals(2 to 3, lut.elemWithinDist(data, 31, 22, 32))
        assertEquals(1 to 4, lut.elemWithinDist(data, 32, 26, 40))
        assertEquals(2 to 4, lut.elemWithinDist(data, 33, 1, 90))
    }

    private fun assertSameIndex(data: IntArray, key: Int, lut: BinaryLut) {
        val i = lut.binarySearch(data, key)
        val j = data.binarySearch(key)
        assertTrue(i == j
                   || (i >= 0 && j >= 0 && data[i] == data[j])
                   || (i < 0 && j < 0 && data[-(i + 1)] == data[-(j + 1)]))
    }

    companion object {
        private val RANDOM = Random()
    }
}