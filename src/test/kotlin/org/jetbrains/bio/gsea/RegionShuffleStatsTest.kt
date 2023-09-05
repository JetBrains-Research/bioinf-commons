package org.jetbrains.bio.gsea

import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.util.withTempBedFile
import org.junit.Test
import kotlin.test.assertEquals

class RegionShuffleStatsTest {
    @Test
    fun loadGenomeAllowedLocationsFilter() {
        val ll: LocationsMergingList = withTempBedFile(
            listOf(
                Location(10, 100, chr1, Strand.PLUS),
                Location(50, 130, chr1, Strand.PLUS),
            )
        ) { bedPath ->
            RegionShuffleStats.loadGenomeAllowedLocationsFilter(bedPath, gq)
        }
        assertEquals(listOf(Location(10,130, chr1)), ll.toList())
    }
    @Test
    fun loadComplementaryToMaskedGenomeRegionsFilterForPlusStrandLocations() {
        val ll: LocationsMergingList = withTempBedFile(
            listOf(
                Location(10, 100, chr1),
                Location(50, 130, chr1),
                Location(200, 230, chr1),
                Location(20, 50, chr2),
            )
        ) { bedPath ->
            RegionShuffleStats.loadComplementaryToMaskedGenomeRegionsFilter(bedPath, gq)
        }
        val expected = arrayListOf(
            Location(0, 10, chr1),
            Location(130, 200, chr1),
            Location(230, chr1.length, chr1),
            Location(0, 20, chr2),
            Location(50, chr2.length, chr2),
        )
        gq.get().drop(2).forEach { chr ->
            expected.add(Location(0, chr.length, chr))
        }
        gq.get().forEach { chr ->
            expected.add(Location(0, chr.length, chr, Strand.MINUS))
        }
        assertEquals(
            expected.sorted(),
            ll.toList().sorted()
        )
    }

    @Test
    fun loadComplementaryToMaskedGenomeRegionsFilterForBothStrandLocations() {
        val ll: LocationsMergingList = withTempBedFile(
            listOf(
                Location(10, 100, chr1, Strand.PLUS),
                Location(50, 130, chr1, Strand.MINUS),
                Location(200, 230, chr1, Strand.PLUS),
            )
        ) { bedPath ->
            RegionShuffleStats.loadComplementaryToMaskedGenomeRegionsFilter(bedPath, gq)
        }
        val expected = arrayListOf(
            Location(0, 10, chr1),
            // XXX: strand is ignored by #RegionShuffleStats.loadComplementaryToMaskedGenomeRegionsFilter()
            Location(130, 200, chr1),
            Location(230, chr1.length, chr1),
            Location(0, chr1.length, chr1, Strand.MINUS)
        )
        gq.get().drop(1).forEach { chr ->
            expected.add(Location(0, chr.length, chr, Strand.PLUS))
            expected.add(Location(0, chr.length, chr, Strand.MINUS))
        }
        assertEquals(
            expected.sorted(),
            ll.toList().sorted()
        )
    }

    @Test
    fun makeAllowedRegionsFilterOnlyAllowed() {
        val loci = listOf(
            Location(10, 100, chr1),
        )
        val llActual = withTempBedFile(
            loci
        ) { path ->
            RegionShuffleStats.makeAllowedRegionsFilter(null, path, gq)
        }!!

        val llExpected = withTempBedFile(
            loci
        ) { path ->
            RegionShuffleStats.loadGenomeAllowedLocationsFilter(path, gq)
        }

        assertEquals(
            llExpected.toList().sorted(),
            llActual.toList().sorted()
        )
    }

    @Test
    fun makeAllowedRegionsFilterOnlyMasked() {
        val loci = listOf(
            Location(10, 100, chr1),
        )
        val llActual = withTempBedFile(
            loci
        ) { path ->
            RegionShuffleStats.makeAllowedRegionsFilter(path, null, gq)
        }!!

        val llExpected = withTempBedFile(
            loci
        ) { path ->
            RegionShuffleStats.loadComplementaryToMaskedGenomeRegionsFilter(path, gq)
        }

        assertEquals(
            llExpected.toList().sorted(),
            llActual.toList().sorted()
        )
    }

    @Test
    fun makeAllowedRegionsFilter() {
        val ll = withTempBedFile(
            listOf(
                Location(10, 100, chr1),
            )
        ) { genomeAllowedMaskedLociPath ->
            withTempBedFile(
                listOf(
                    Location(50, 80, chr1, Strand.PLUS),
                )
            ) { genomeMaskedLociPath ->
                RegionShuffleStats.makeAllowedRegionsFilter(genomeMaskedLociPath, genomeAllowedMaskedLociPath, gq)
            }
        }!!

        val expected = arrayListOf(
            Location(10, 50, chr1),
            Location(80, 100, chr1),
        )
        assertEquals(
            expected.sorted(),
            ll.toList().sorted()
        )
    }

    companion object {
        internal val gq: GenomeQuery = GenomeQuery(Genome["to1"])
        internal val chr1: Chromosome = gq.get()[0]
        internal val chr2: Chromosome = gq.get()[1]
    }
}
