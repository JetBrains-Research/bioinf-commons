package org.jetbrains.bio.genome.coverage

import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.processReads
import org.jetbrains.bio.genome.format.toBedEntry
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals

class FragmentSizeTest {
    @Test
    fun testDetectFragmentSize() {
        withTempFile("track", ".bed") { trackPath ->
            val bedFormat = BedFormat()
            bedFormat.print(trackPath).use {
                it.print(Location(0, 50, chromosome1, Strand.PLUS).toBedEntry())
                it.print(Location(0, 50, chromosome1, Strand.PLUS).toBedEntry())
                it.print(Location(0, 50, chromosome1, Strand.PLUS).toBedEntry())
                it.print(Location(10, 60, chromosome1, Strand.PLUS).toBedEntry())
                it.print(Location(20, 70, chromosome1, Strand.PLUS).toBedEntry())
                it.print(Location(51, 101, chromosome1, Strand.MINUS).toBedEntry())
                it.print(Location(51, 101, chromosome1, Strand.MINUS).toBedEntry())
                it.print(Location(51, 101, chromosome1, Strand.MINUS).toBedEntry())
                it.print(Location(111, 161, chromosome1, Strand.MINUS).toBedEntry())
                it.print(Location(121, 171, chromosome1, Strand.MINUS).toBedEntry())
            }

            val builder = SingleEndCoverage.builder(genomeQuery).apply {
                processReads(genomeQuery, trackPath) {
                    process(it)
                }
            }
            // this is required to sort the builder data
            builder.build(false)
            // the fragment detection method should ignore duplicate tags,
            // so 100 is a wrong answer!
            assertEquals(
                150,
                FragmentSize.detectFragmentSize(builder.data, 20.0)
            )
        }
    }

    companion object {
        internal var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
        internal var chromosome1: Chromosome = genomeQuery.get()[0]
    }
}