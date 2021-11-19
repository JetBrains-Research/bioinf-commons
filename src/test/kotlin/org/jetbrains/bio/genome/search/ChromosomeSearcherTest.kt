package org.jetbrains.bio.genome.search

import htsjdk.samtools.util.SequenceUtil
import junit.framework.TestCase
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.search.sa.SuffixArraySearcher
import org.jetbrains.bio.genome.sequence.Nucleotide
import org.jetbrains.bio.genome.sequence.NucleotideAlternative
import org.jetbrains.bio.genome.sequence.asNucleotideSequence
import java.util.*

class ChromosomeSearcherTest : TestCase() {
    private var randomized: Boolean = true // if false, the following can be ignored
    private var fromChromosome: Boolean = false
    private var strand: Strand? = null
    private var chomosomePos: Int = 0
    private val chromosome = Chromosome(Genome["to1"], "chr1")
    private var searcher: ChromosomeSearcher? = null
    private var pattern: String? = null
    private var mutatePattern = false
    private var actualPattern: String? = null

    @Throws(Exception::class)
    public override fun setUp() {
        super.setUp()
        searcher = ChromosomeSearcher(SuffixArraySearcher(chromosome), MAX_MISMATCH, chromosome)
    }

    private fun mutatePattern() {
        val generator = Random(System.currentTimeMillis())
        val newPattern = pattern!!.toCharArray()
        for (i in 0 until MAX_MISMATCH) {
            val pos = generator.nextInt(pattern!!.length)
            var newChar: Char
            do {
                newChar = Nucleotide.ALPHABET[generator.nextInt(Nucleotide.ALPHABET.size)]
            } while (newChar == newPattern[pos])
            newPattern[pos] = newChar
        }
        actualPattern = pattern
        pattern = String(newPattern)
    }

    private fun generatePattern() {
        if (!randomized) {
            pattern = if (!fromChromosome)
                FIXED_PATTERN
            else
                chromosome.sequence.substring(
                    chomosomePos, chomosomePos + FIXED_PATTERN.length, strand!!
                )
        } else {
            val generator = Random(System.currentTimeMillis())
            if (fromChromosome) {
                var res: String
                do {
                    chomosomePos = generator.nextInt(chromosome.length - FIXED_PATTERN.length)
                    strand = if (generator.nextBoolean()) Strand.PLUS else Strand.MINUS
                    res = chromosome.sequence.substring(
                        chomosomePos, chomosomePos + FIXED_PATTERN.length, strand!!
                    )
                } while (res.indexOf('n') != -1 || res.indexOf('N') != -1)
                pattern = res
            } else {
                val res = CharArray(FIXED_PATTERN.length)
                for (i in FIXED_PATTERN.indices)
                    res[i] = Nucleotide.ALPHABET[generator.nextInt(Nucleotide.ALPHABET.size)]
                pattern = String(res)
            }
        }
        if (mutatePattern)
            mutatePattern()
    }

    @Throws(Exception::class)
    fun testCorrectness() {
        mutatePattern = false

        randomized = false
        fromChromosome = false
        generatePattern()
//        println("testCorrectness on a fixed pattern: " + pattern!!)
        performTestCorrectness()

        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testCorrectness on a pattern from chromosome: " + pattern!!)
        println("Position: $chomosomePos $strand")
        performTestCorrectness()

        mutatePattern = true
        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testCorrectness on a changed pattern from chromosome: $pattern (was: $actualPattern)")
        println("Position: $chomosomePos $strand")
        performTestCorrectness()

        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testCorrectness on a random pattern: " + pattern!!)
        performTestCorrectness()
    }

    private fun performTestCorrectness() {
        searcher!!.find(pattern!!).forEach { match ->
            val matchingSeq = chromosome.sequence.substring(match.startOffset, match.endOffset, match.strand)
            var mismatchCount = 0
            for (i in FIXED_PATTERN.indices) {
                if (!NucleotideAlternative
                        .fromChar(pattern!![i])
                        .match(Nucleotide.valueOf(matchingSeq[i].toString().uppercase()).byte)
                )
                    ++mismatchCount
            }
//            println("Found pattern: " + matchingSeq + " at " + match.startOffset + ' ' + match.strand)
            assertTrue(
                "testCorrectness(): mismatch count exceeds maximum allowed value: $mismatchCount",
                mismatchCount <= MAX_MISMATCH
            )
        }
    }

    @Throws(Exception::class)
    fun testSpecificity() {
        mutatePattern = false
        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testSpecificity on a pattern from chromosome: " + pattern!!)
//        println("Position: $chomosomePos $strand")
        performTestSpecificity()

        mutatePattern = true
        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testSpecificity on a changed pattern from chromosome: $pattern (was: $actualPattern)")
//        println("Position: $chomosomePos $strand")
        performTestSpecificity()
    }

    private fun performTestSpecificity() {
        if (!fromChromosome) {
            return
        }
        assertTrue("testSpecificity(): the pattern was not found at the expected position",
            searcher!!.find(pattern!!).anyMatch { match ->
                match.startOffset == chomosomePos && match.strand === strand
            })
    }

    fun testDistinct() {
        mutatePattern = false

        randomized = false
        fromChromosome = false
        generatePattern()
//        println("testDistinct on a fixed pattern: " + pattern!!)
        performTestDistinct()

        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testDistinct on a pattern from chromosome: " + pattern!!)
//        println("Position: $chomosomePos $strand")
        performTestDistinct()

        mutatePattern = true
        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testDistinct on a changed pattern from chromosome: $pattern (was: $actualPattern)")
//        println("Position: $chomosomePos $strand")
        performTestDistinct()

        randomized = true
        fromChromosome = true
        generatePattern()
//        println("testDistinct on a random pattern: " + pattern!!)
        performTestDistinct()
    }

    private fun performTestDistinct() {
        val results = HashSet<Location>()
        searcher!!.find(pattern!!).forEach { m ->
            assertTrue("testDistinct(): the lazy iterator produces duplicate entries.", results.add(m))
        }
    }

    fun testUniqueSearch() {
        val searcher = SimpleSearcher(chromosome.sequence.toString())
        for (i in 0..99) {
            val s = Random().nextInt(chromosome.length - 20)
            val pattern = chromosome.sequence.substring(s, s + 20)
            if (pattern.indexOf('n') != -1) {
                continue
            }
            val location = this.searcher!!.findUnique(pattern)
            if (location != null) {
//                println("Found $pattern at $location")
                assertEquals(
                    pattern,
                    chromosome.sequence.substring(location.startOffset, location.endOffset, location.strand)
                )

                assertEquals(
                    1L, searcher.find(pattern.asNucleotideSequence()).count() +
                            searcher.find(SequenceUtil.reverseComplement(pattern).asNucleotideSequence()).count()
                )
            }
        }
    }

    companion object {
        private const val MAX_MISMATCH = 2
        private const val FIXED_PATTERN = "gcagtgagccaggtaTtc"
    }
}
