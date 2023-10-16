package org.jetbrains.bio.util

import org.junit.Test
import java.text.DecimalFormatSymbols
import java.util.*
import kotlin.test.assertEquals

class StringUtilTest {

    @Test
    fun testSeparatorToSystem() {
        assertEquals("1\\1\\1\\1\\1\\1\\", "1/1/1/1/1/1/".separatorsToSystem('\\'))
        assertEquals("1\\2\\1", "1/2\\1".separatorsToSystem('\\'))


        assertEquals("1/1/1/1/1/1/", "1\\1\\1\\1\\1\\1\\".separatorsToSystem('/'))
        assertEquals("1/2/1", "1/2\\1".separatorsToSystem('/'))


        assertEquals("100", "100".separatorsToSystem('/'))
        assertEquals("100", "100".separatorsToSystem('\\'))
    }

    @Test
    fun testIntFormatIntLongNumber() {
        val gs = DecimalFormatSymbols(Locale.getDefault()).groupingSeparator
        assertEquals("0", 0.formatLongNumber())
        assertEquals("10", 10.formatLongNumber())
        assertEquals("1${gs}000", 1000.formatLongNumber())
        assertEquals("10${gs}000", 10000.formatLongNumber())

        assertEquals("0", 0.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("10", 10.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("1 000", 1000.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("10 000", 10000.formatLongNumber(universalThousandsSeparator = true))
    }

    @Test
    fun testLongFormatLongNumber() {
        val gs = DecimalFormatSymbols(Locale.getDefault()).groupingSeparator
        assertEquals("0", 0L.formatLongNumber())
        assertEquals("10", 10L.formatLongNumber())
        assertEquals("1${gs}000", 1000L.formatLongNumber())
        assertEquals("10${gs}000", 10000L.formatLongNumber())
        assertEquals("1${gs}000${gs}000", 1000000L.formatLongNumber())

        assertEquals("0", 0L.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("10", 10L.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("1 000", 1000L.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("10 000", 10000L.formatLongNumber(universalThousandsSeparator = true))
        assertEquals("1 000 000", 1000000L.formatLongNumber(universalThousandsSeparator = true))
    }
}