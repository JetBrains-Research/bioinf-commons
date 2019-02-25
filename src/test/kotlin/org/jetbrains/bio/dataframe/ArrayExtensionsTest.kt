package org.jetbrains.bio.dataframe

import org.junit.Assert
import org.junit.Test

class ArrayExtensionsTest {
    @Test fun testSortedReorder() {
        val values = intArrayOf(42, 2, -1, 0, 4, 2);
        val indices = values.argSort();
        Assert.assertArrayEquals(intArrayOf(2, 3, 1, 5, 4, 0), indices);
        values.reorder(indices)
        Assert.assertArrayEquals(intArrayOf(-1, 0, 2, 2, 4, 42), values);
    }
}
