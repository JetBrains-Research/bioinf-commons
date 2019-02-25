package org.jetbrains.bio.dataframe

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class BooleanColumnTest {
    @Test fun testResizeSmaller() {
        val c = BooleanColumn("mask", BitterSet(10))
        c.data.set(0, 10)

        val resized = c.resize(5)
        assertEquals(5, resized.size)

        for (i in 0 until resized.size) {
            assertTrue(resized.data[i])
        }
        for (i in resized.size until c.size) {
            assertFalse(resized.data[i])
        }
    }

    @Test fun testResizeLarger() {
        val c = BooleanColumn("mask", BitterSet(10))
        c.data.set(0, 10)

        val resized = c.resize(20)
        assertEquals(20, resized.size)

        for (i in 0 until c.size) {
            assertTrue(resized.data[i])
        }
        for (i in c.size until resized.size) {
            assertFalse(resized.data[i])
        }
    }
}