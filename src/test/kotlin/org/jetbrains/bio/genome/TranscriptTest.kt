package org.jetbrains.bio.genome

import com.google.common.collect.Ordering
import org.jetbrains.bio.util.name
import org.jetbrains.bio.util.toPath
import org.junit.Test
import java.util.*
import kotlin.test.*

class TranscriptTest {
    private val chromosome = Chromosome(Genome["to1"], "chr1")
    private val RANDOM = Random(1234)

    @Test
    fun exonsIntronsSorted() {
        val transcripts = chromosome.transcripts.filter { it.exons.size > 1 }
        val ordering = Ordering.natural<Location>()

        val transcriptPlus = transcripts.first { it.strand.isPlus() }
        assertTrue(ordering.isOrdered(transcriptPlus.exons), "Exons aren't ordered")
        assertTrue(ordering.isOrdered(transcriptPlus.introns), "Introns aren't ordered")

        val transcriptMinus = transcripts.first { it.strand.isMinus() }
        assertTrue(ordering.isOrdered(transcriptMinus.exons), "Exons aren't ordered")
        assertTrue(ordering.isOrdered(transcriptMinus.introns), "Introns aren't ordered")
    }

    @Test
    fun equalsHashCode() {
        val gene1 = Transcript(
            "foo", "foo_gene", "symbol1",
            Location(0, 100, chromosome, Strand.PLUS), null, -1, listOf()
        )
        val gene2 = Transcript(
            "foo", "foo_gene", "symbol2",
            Location(0, 200, chromosome, Strand.PLUS), null, -1, listOf()
        )
        val gene3 = Transcript(
            "foo2", "foo_gene", "symbol3",
            Location(0, 300, chromosome, Strand.PLUS), null, -1, listOf()
        )

        assertEquals(gene1.hashCode(), gene2.hashCode())
        assertNotEquals(gene1.hashCode(), gene3.hashCode())
        assertEquals(gene1, gene2)
        assertNotEquals(gene1, gene3)
        assertNotEquals(gene2, gene3)
    }


    @Test
    fun transcriptsSorted() {
        assert(Ordering.natural<Int>().isOrdered(chromosome.transcripts.map { it.location.get5Bound() }))
    }

    @Test
    fun transcripts5IndexSorted() {
        val bounds5 = Transcripts.bound5Index(chromosome.genome, false)[chromosome]!!.second
        assert(Ordering.natural<Int>().isOrdered(bounds5.toList()))
    }

    @Test
    fun transcripts5AllGenes() {
        val (transcripts, bounds5) = Transcripts.bound5Index(chromosome.genome, false)[chromosome]!!
        assertEquals(transcripts.size, bounds5.size)
        assertTrue(transcripts.any { it.isCoding }, "No coding genes found")
        assertTrue(transcripts.any { !it.isCoding }, "No non-coding genes found")
    }

    @Test
    fun transcripts5CodingGenes() {
        val (transcripts, bounds5) = Transcripts.bound5Index(chromosome.genome, true)[chromosome]!!
        assertEquals(transcripts.size, bounds5.size)
        assertTrue(transcripts.all { it.isCoding }, "Coding genes query returned some non-coding genes")
    }

    /**
     * The simulated transcripts don't intersect, so there's only one transcript closest to the start.
     */
    @Test
    fun associatedTranscriptsStart() {
        val result = Transcripts.associatedTranscriptsSingle(Location(0, 1, chromosome), limit = chromosome.length)
        assertEquals(
            result.size, 1, "Expected a single transcript to be associated with the chromosome start, " +
                    "but got ${result.size} transcripts"
        )
        assertEquals(chromosome.transcripts.first(), result.first())
    }

    /**
     * The simulated transcripts don't intersect, so there's only one transcript framing the start.
     */
    @Test
    fun associatedTranscripts2Start() {
        val result = Transcripts.associatedTranscriptsTwo(Location(0, 1, chromosome), limit = chromosome.length)
        assertEquals(
            result.size, 1, "Expected a single transcript to be associated with the chromosome start, " +
                    "but got ${result.size} transcripts"
        )
        assertEquals(chromosome.transcripts.first(), result.first())
    }

    /**
     * The simulated transcripts don't intersect, so there's only one transcript closest to the end.
     */
    @Test
    fun associatedTranscriptsEnd() {
        val result = Transcripts.associatedTranscriptsSingle(
            Location(chromosome.length - 1, chromosome.length, chromosome),
            limit = chromosome.length
        )
        assertEquals(
            result.size, 1, "Expected a single transcript to be associated with the chromosome end, " +
                    "but got ${result.size} transcripts"
        )
        assertEquals(chromosome.transcripts.last(), result.first())
    }

    /**
     * The simulated transcripts don't intersect, so there's only one transcript framing the end.
     */
    @Test
    fun associatedTranscripts2End() {
        val result = Transcripts.associatedTranscriptsTwo(
            Location(chromosome.length - 1, chromosome.length, chromosome),
            limit = chromosome.length
        )
        assertEquals(
            result.size, 1, "Expected a single transcript to be associated with the chromosome end, " +
                    "but got ${result.size} transcripts"
        )
        assertEquals(chromosome.transcripts.last(), result.first())
    }

    @Test
    fun associatedTranscriptsLargeLocations() {
        for (i in 0..999) {
            val start = RANDOM.nextInt((chromosome.length * 0.99).toInt())
            val end = start + RANDOM.nextInt((chromosome.length * 0.01).toInt())
            val location = Location(start, end, chromosome)
            doTestAssociatedTranscriptsSingle(location)
            doTestAssociatedTranscriptsTwo(location)
            doTestAssociatedTranscriptsPlus(location)
        }
    }

    @Test
    fun associatedTranscriptsSmallLocations() {
        for (i in 0..999) {
            val start = RANDOM.nextInt((chromosome.length - 1))
            val end = start + 1
            val location = Location(start, end, chromosome)
            doTestAssociatedTranscriptsSingle(location)
            doTestAssociatedTranscriptsTwo(location)
            doTestAssociatedTranscriptsPlus(location)
        }
    }

    @Test
    fun associatedTranscriptsEquidistant() {
        // equidistant are only original transcripts, not additional splicing variants
        val transcripts = chromosome.transcripts.filter { it.ensemblId.endsWith("_1") }

        for (i in 0..transcripts.size - 2) {
            val firstTSS = transcripts[i].location.get5Bound()
            val secondTSS = transcripts[i + 1].location.get5Bound()
            if ((firstTSS - secondTSS) % 2 == 0) {
                val midpoint = (firstTSS + secondTSS) / 2
                val location = Location(midpoint, midpoint + 1, chromosome)
                val associatedTranscripts = Transcripts.associatedTranscriptsSingle(location, limit = chromosome.length)
                assertEquals(
                    2, associatedTranscripts.size,
                    "Expected to get both TSS equidistant to $midpoint, but got ${associatedTranscripts.size}"
                )
                assertEquals(
                    arrayListOf(transcripts[i], transcripts[i + 1]), associatedTranscripts,
                    "Expected ${transcripts[i].ensemblId} and ${transcripts[i + 1].ensemblId}, but got " +
                            "${associatedTranscripts.map { it.ensemblId }}"
                )
                doTestAssociatedTranscriptsSingle(location)
            }
        }
    }

    @Test
    fun associatedTranscripts2Equidistant() {
        val transcripts = chromosome.transcripts
        for (i in 0..transcripts.size - 2) {
            val firstTSS = transcripts[i].location.get5Bound()
            val secondTSS = transcripts[i + 1].location.get5Bound()
            val midpoint = (firstTSS + secondTSS) / 2
            val location = Location(midpoint, midpoint + 1, chromosome)
            val associatedTranscripts = Transcripts.associatedTranscriptsTwo(location, limit = chromosome.length)
            assertEquals(
                2, associatedTranscripts.size,
                "Expected to get both TSS framing $midpoint, but got ${associatedTranscripts.size}"
            )
//            val expectedList = arrayListOf<Transcript>()
//            for (k in 0 until associatedTranscripts.size) {
//                expectedList.add(transcripts[i+k])
//            }
            assertEquals(
                arrayListOf(transcripts[i], transcripts[i + 1]), associatedTranscripts,
                "Expected ${transcripts[i].ensemblId} and ${transcripts[i + 1].ensemblId}, but got " +
                        "${associatedTranscripts.map { it.ensemblId }}"
            )
            doTestAssociatedTranscriptsTwo(location)
        }
    }

    fun doTestAssociatedTranscriptsSingle(location: Location) {
        val actual = Transcripts.associatedTranscriptsSingle(location, limit = chromosome.length)
        assert(actual.isNotEmpty())

        val transcripts = chromosome.transcripts
        val dists = transcripts.map { Transcripts.greatDistance(it, location) }
        val minDist = dists.minOrNull()!!

        val candidates = dists.mapIndexed { i, d -> if (d == minDist) transcripts[i] else null }.filterNotNull().toSet()

        assert(actual.toSet() == candidates) {
            "Actual transcripts (${actual.map { it.ensemblId }}) " +
                    "not equal to expected ones (${candidates.map { it.ensemblId }})"
        }
    }

    fun doTestAssociatedTranscriptsTwo(location: Location) {
        val actual = Transcripts.associatedTranscriptsTwo(location, limit = chromosome.length)
        assert(actual.isNotEmpty())

        val transcripts = chromosome.transcripts
        val dists = transcripts.map { Transcripts.signedGreatDistance(it, location) }
        val leftDist = dists.filter { it < 0 }.maxOrNull() ?: 0
        val rightDist = dists.filter { it > 0 }.minOrNull() ?: 0

        val candidates =
            dists.mapIndexed { i, d -> if (d in leftDist..rightDist) transcripts[i] else null }.filterNotNull().toSet()

        assert(actual.toSet() == candidates) {
            "Actual transcripts (${actual.map { it.ensemblId }}) " +
                    "not equal to expected ones (${candidates.map { it.ensemblId }})"
        }
    }


    internal fun regulatoryDomainBasal(t: Transcript): Location =
        RelativePosition.FIVE_PRIME.of(t.location, -5000, 1001)

    val regulatoryDomainsExt: List<Location> by lazy {
        val transcripts = chromosome.transcripts
        transcripts.indices.map {
            regulatoryDomainExt(it)
        }
    }

    internal fun regulatoryDomainExt(i: Int): Location {
        val transcripts = chromosome.transcripts
        val basal = regulatoryDomainBasal(transcripts[i])
        val start = Math.min(
            transcripts.map { regulatoryDomainBasal(it).endOffset }
                .filter { it < basal.endOffset }.maxOrNull() ?: chromosome.range.startOffset,
            basal.startOffset
        )
        val end = Math.max(
            transcripts.map { regulatoryDomainBasal(it).startOffset }
                .filter { it > basal.startOffset }.minOrNull() ?: chromosome.range.endOffset,
            basal.endOffset
        )
        return Location(start, end, chromosome, transcripts[i].strand)
    }

    private fun doTestAssociatedTranscriptsPlus(location: Location) {
        val actual = Transcripts.associatedTranscriptsPlus(location, distal = chromosome.length)
        assert(actual.isNotEmpty())

        val transcripts = chromosome.transcripts
        val midpoint = (location.startOffset + location.endOffset) / 2
        val candidates = transcripts.filterIndexed { i, _ ->
            regulatoryDomainsExt[i].let { midpoint in it.startOffset until it.endOffset }
        }.toSet()
        assert(actual.toSet() == candidates) {
            "Actual transcripts (${actual.map { it.ensemblId }}) " +
                    "not equal to expected ones (${candidates.map { it.ensemblId }})"
        }
    }

    @Test
    fun utr5() {
        // PLUS:
        var t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 83,
            listOf(Range(12, 15), Range(18, 25), Range(70, 90))
        )
        assertEquals(listOf(Range(12, 15), Range(18, 20)), t.utr5.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 83,
            listOf(Range(10, 15), Range(18, 20), Range(20, 30), Range(70, 90))
        )
        assertEquals(listOf(Range(10, 15), Range(18, 20)), t.utr5.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 83,
            listOf(Range(20, 30), Range(70, 90))
        )
        assertEquals(listOf(), t.utr5)

        // MINUS:
        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(18, 25), Range(70, 90), Range(92, 95))
        )
        assertEquals(listOf(Range(81, 90), Range(92, 95)), t.utr5.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(18, 80), Range(83, 90), Range(92, 95))
        )
        assertEquals(listOf(Range(83, 90), Range(92, 95)), t.utr5.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(18, 80))
        )
        assertEquals(listOf(), t.utr5.map(Location::toRange))
    }

    @Test
    fun utr3() {
        // PLUS:
        var t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 83,
            listOf(Range(18, 25), Range(70, 90), Range(92, 95))
        )
        assertEquals(listOf(Range(83, 90), Range(92, 95)), t.utr3.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 85,
            listOf(Range(18, 81), Range(83, 90), Range(92, 95))
        )
        //stop codon is: 80,83,84
        assertEquals(listOf(Range(85, 90), Range(92, 95)), t.utr3.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.PLUS),
            Range(20, 80), 83,
            listOf(Range(28, 40), Range(43, 80))
        )
        assertEquals(listOf(), t.utr3.map(Location::toRange))

        // MINUS:
        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(12, 15), Range(17, 25), Range(70, 90))
        )
        assertEquals(listOf(Range(12, 15), Range(17, 18)), t.utr3.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(12, 15), Range(18, 25), Range(70, 90))
        )
        assertEquals(listOf(Range(12, 15)), t.utr3.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 15,
            listOf(Range(10, 15), Range(16, 18), Range(20, 30), Range(70, 90))
        )
        // stop codon is: 20,17,16
        assertEquals(listOf(Range(10, 15)), t.utr3.map(Location::toRange))

        t = Transcript(
            "1", "2", "3",
            Location(10, 100, chromosome, Strand.MINUS),
            Range(21, 81), 17,
            listOf(Range(18, 30), Range(70, 90))
        )
        assertEquals(listOf(), t.utr3)

    }

    @Test
    fun organismIsValid() {
        // Just load test organism and ensure that fields initialized correctly + auto caches recreation works.
        // (e.g. if Transcript was changed, but to wasn't re-created)

        val chr = Chromosome(Genome["to1"], "chr1")
        var t = chr.transcripts.stream().filter { it.ensemblId == "ENSTSIMGENE.CHR1.0_1" }.findFirst().get()
        assertTrue(!t.isCoding, "ENSTSIMGENE.CHR1.0_1 should be non-coding")
        assertNull(t.cdsRange, "CDS of a non-coding gene should be null")

        t = chr.transcripts.stream().filter { it.ensemblId == "ENSTSIMGENE.CHR1.1_1" }.findFirst().get()
        assertTrue(t.isCoding, "ENSTSIMGENE.CHR1.1 should be coding")
        assertNotNull(t.cdsRange, "CDS of a coding gene shouldn't be null")
    }

    @Test
    fun cachedTranscriptsJsonPath() {
        // only 'to1' doesn't require annotations downloading
        val genome = Genome["to1"]

        val gtf = genome.genesGtfPath(false)
        assertEquals("genes.gtf.gz", gtf.name)

        // XXX json files for normal genomes are bundled annotations.tar.gz
        // do not change this test data without updating files in annotations.tar.gz
        assertEquals(
            "genes.gtf.json.gz",
            Transcripts.cachedTranscriptsJsonPath(gtf).name
        )

        assertEquals(
            "foo/Homo_sapiens.GRCh37.87.gtf.json.gz".toPath(),
            Transcripts.cachedTranscriptsJsonPath("foo/Homo_sapiens.GRCh37.87.gtf.gz".toPath())
        )
    }
}
