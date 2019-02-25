package org.jetbrains.bio.genome

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertSame

class ChromosomeTest {
    @Test fun reuseInvoke() {
        assertSame(Chromosome("to1", "chr1"), Chromosome("to1", "chr1"))
        assertSame(Chromosome(Genome["to1"], "chr1"), Chromosome("to1", "chr1"))
        assertSame(Chromosome(Genome["to1"], "1"), Chromosome("to1", "chr1"))
    }

    @Test fun reuseGson() {
        val chromosome = Chromosome("to1", "chr1")
        assertSame(with(Chromosome.ADAPTER) { fromJson(toJson(chromosome)) }, chromosome)
    }

    @Test fun canonicalName() {
        assertEquals("chr1", Chromosome("to1", "chr1").name)
        assertEquals("chr1", Chromosome("to1", "1").name)
    }
}