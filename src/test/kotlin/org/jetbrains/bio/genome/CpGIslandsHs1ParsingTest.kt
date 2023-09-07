package org.jetbrains.bio.genome

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.assertNull

class CpGIslandsHs1ParsingTest {

    private val chr1 = Chromosome(Genome["to1"], "chr1")

    @Test
    fun shouldParseAnIslandWithAnAdditionalField() {
        val result = CpGIslandsHs1.parseIsland(chr1, 2, 2705, "CpG: 26\t373\t26\t205\t13.9\t55.0\t0.93\taddedLater")

        val expectedIsland = CpGIsland(26, 205, 0.93, Location(2, 2705, chr1, Strand.PLUS))
        assertNotNull(result)
        assertEquals(expectedIsland, result)
    }

    @Test
    fun shouldSkipInvalidLine() {
        val result = CpGIslandsHs1.parseIsland(chr1, 2, 2705, "invalid;line")

        assertNull(result)
    }
}