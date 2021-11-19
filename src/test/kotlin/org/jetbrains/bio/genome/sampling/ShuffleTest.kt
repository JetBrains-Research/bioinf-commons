package org.jetbrains.bio.genome.sampling

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail


class ShuffleTest {
    val genomeQuery = GenomeQuery(Genome["to1"])
    val chromosomes = genomeQuery.get()

    private fun checkShuffled(regions: List<ChromosomeRange>, shuffled: List<ChromosomeRange>) {
        val l1 = regions.map { it.length() }.sorted().toList()
        val l2 = shuffled.map { it.length() }.sorted().toList()
        assertEquals(l1, l2, "Shuffled regions must have same length set")
        for (i1 in shuffled.indices) {
            for (i2 in shuffled.indices) {
                if (i1 == i2) {
                    continue
                }
                if (shuffled[i1].chromosome == shuffled[i2].chromosome &&
                    shuffled[i1].toRange() intersects shuffled[i2].toRange()
                ) {
                    fail("Regions must not intersect")
                }
            }
        }
    }

    @Test
    fun shuffleRegionsWithBackgroud() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val background = listOf(
            ChromosomeRange(100, 300, chromosomes[0]),
            ChromosomeRange(1000, 2000, chromosomes[1]),
            ChromosomeRange(2000, 3000, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(genomeQuery, regions, background)

        checkShuffled(regions, shuffled)
        for (r in shuffled) {
            var inBackground = false
            for (b in background) {
                if (r.chromosome == b.chromosome && r.startOffset in b.toRange()) {
                    inBackground = true
                }
            }
            assertTrue(inBackground, "All region start must be in backgroud.")
        }
    }


    @Test
    fun shuffleRegions() {
        val regions = listOf(
            ChromosomeRange(0, 100, chromosomes[0]),
            ChromosomeRange(0, 200, chromosomes[1]),
            ChromosomeRange(0, 300, chromosomes[2])
        )

        val shuffled = shuffleChromosomeRanges(genomeQuery, regions)
        checkShuffled(regions, shuffled)
    }

}