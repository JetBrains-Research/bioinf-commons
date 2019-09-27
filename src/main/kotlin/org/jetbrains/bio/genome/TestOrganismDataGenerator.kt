package org.jetbrains.bio.genome

import com.google.common.math.IntMath
import gnu.trove.list.array.TFloatArrayList
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.big.FixedStepSection
import org.jetbrains.bio.genome.sequence.Nucleotide
import org.jetbrains.bio.genome.sequence.TwoBitWriter
import org.jetbrains.bio.io.FastaRecord
import org.jetbrains.bio.io.write
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempFile
import java.util.*
import java.util.concurrent.ThreadLocalRandom
import java.util.stream.Collectors

/** A generator for the fake "to1" genome build. */
@Suppress("unused")
object TestOrganismDataGenerator {
    private val LOG = Logger.getLogger(TestOrganismDataGenerator::class.java)

    private val CHROMOSOMES_SIZES = mapOf(
            "chr1" to IntMath.pow(10, 7),
            "chr2" to IntMath.pow(10, 6),
            "chr3" to IntMath.pow(10, 6),
            "chrX" to IntMath.pow(10, 6),
            "chrM" to IntMath.pow(10, 6)
    )

    @JvmStatic
    fun main(args: Array<String>) {
        LOG.info("User home: '${System.getProperty("user.home", "n/a")}'")

        // Ensure genomesPaths initialized
        Configuration.genomesPath.createDirectories()

        LOG.info("Genomes path: ${Configuration.genomesPath}")
        val dataPath = Configuration.genomesPath / Genome.TEST_ORGANISM_BUILD
        // Ensure data path initialized
        dataPath.createDirectories()
        LOG.info("Generating test organism into $dataPath")

        LOG.info("Generating chrom.sizes")
        val build = Genome.TEST_ORGANISM_BUILD
        val chromSizesPath = Configuration.genomesPath / build / "$build.chrom.sizes"
        CSVFormat.TDF.print(chromSizesPath.bufferedWriter()).use { csvPrinter ->
            CHROMOSOMES_SIZES.forEach { name, size ->
                csvPrinter.printRecord(name, size.toString())
            }
        }

        LOG.info("Generating other files")
        val genome = Genome[Genome.TEST_ORGANISM_BUILD]
        generateSequence(genome)
        generateTranscripts(genome)
        generateCentromere(genome)
        generateCytobands(genome)
        generateRepeats(genome)
        generateCGI(genome)
        generateMapability(genome)
        LOG.info("Done")
    }

    /**
     * Generates 2bit sequences.
     *
     * @param maxGaps maximum number of gaps (sequences of consecutive Ns)
     *                on a chromosome.
     * @param maxGapLength maximum gap length.
     */
    private fun generateSequence(genome: Genome, maxGaps: Int = 10, maxGapLength: Int = 100) {
        LOG.info("Generating FASTA sequence")
        withTempFile(genome.build, ".fa") { fastaPath ->
            val r = ThreadLocalRandom.current()
            CHROMOSOMES_SIZES.asSequence().map {
                val (name, length) = it
                val sequence = r.ints(length.toLong(), 0, Nucleotide.ALPHABET.size)
                        .mapToObj { Nucleotide.ALPHABET[it].toString() }
                        .collect(Collectors.joining())
                        .toCharArray()

                val numGaps = Math.min(r.nextInt(maxGaps + 1), 3)
                for (i in 0 until numGaps) {
                    // We could do better here by generating non-overlapping
                    // ranges, but it's not worth the effort.
                    val gapLength = r.nextInt(maxGapLength) + 1
                    val gapOffset = r.nextInt(length - gapLength)
                    (gapOffset until gapOffset + gapLength).forEach {
                        sequence[it] = 'n'
                    }
                }

                FastaRecord(name, String(sequence))
            }.asIterable().write(fastaPath)

            LOG.info("Converting FASTA sequence to 2bit")
            TwoBitWriter.convert(fastaPath, genome.twoBitPath(false))
        }
    }

    /**
     * Generates JSON-serialized transcript annotations.
     *
     * @see Transcripts
     */
    private fun generateTranscripts(genome: Genome) {
        LOG.info("Generating transcript annotations")

        val r = ThreadLocalRandom.current()

        genome.genesGtfPath(false).bufferedWriter().use { w ->

            val transcripts = ArrayList<Transcript>()

            for (name in CHROMOSOMES_SIZES.keys.sorted()) {
                val chromosome = Chromosome(genome, name)
                val length = CHROMOSOMES_SIZES[name]!!

                var currentEnd = 0
                var geneNumber = 0
                while (length - currentEnd > 6e4) {
                    val geneSymbol = "simgene.$name.$geneNumber".toUpperCase()
                    val currentStart = currentEnd + 10000 + r.nextInt(40000)
                    val currentLength = 100 + Math.min(length - currentStart - 10000, 10000)
                    val strand = if (r.nextBoolean()) Strand.PLUS else Strand.MINUS
                    val transcript = generateTranscript(geneNumber, geneSymbol, chromosome, strand,
                            currentStart, currentLength)
                    transcripts.add(transcript)

                    currentEnd = currentStart + currentLength
                    geneNumber++
                }
            }

            writeGtf(w, transcripts)
        }
    }

    /** Here be dragons! */
    private fun generateTranscript(geneNumber: Int, geneSymbol: String,
                                   chromosome: Chromosome, strand: Strand,
                                   offset: Int, length: Int): Transcript {
        val ensemblId = "ENST$geneSymbol"

        val minExonLength = 3
        val minIntronLength = 2
        val r = ThreadLocalRandom.current()
        val exonCount = 1 + r.nextInt(8)
        val leftFlankingIntron = r.nextBoolean()
        val rightFlankingIntron = r.nextBoolean()
        val exonStarts = ArrayList<Int>(exonCount)
        val exonEnds = ArrayList<Int>(exonCount)
        var lengthPool = (length - minExonLength * exonCount - minIntronLength * (exonCount - 1)
                - (if (leftFlankingIntron) minIntronLength else 0)
                - if (rightFlankingIntron) minIntronLength else 0)
        exonStarts.add(0, if (leftFlankingIntron) offset + minIntronLength + r.nextInt(lengthPool) else offset)
        if (leftFlankingIntron) {
            val intronLen = exonStarts[0] - offset
            lengthPool -= intronLen - minIntronLength // length pool already takes in consideration min lengths
        }
        for (i in 0 until exonCount - 1) {
            exonEnds.add(i, exonStarts[i] + minExonLength + r.nextInt(lengthPool))
            val exonLen = exonEnds[i] - exonStarts[i]
            lengthPool -= exonLen - minExonLength // length pool already takes in consideration min lengths
            exonStarts.add(i + 1, exonEnds[i] + minIntronLength + r.nextInt(lengthPool))
            val intronLen = exonStarts[i + 1] - exonEnds[i]
            lengthPool -= intronLen - minIntronLength
        }

        exonEnds.add(exonCount - 1,
                if (rightFlankingIntron)
                    exonStarts[exonCount - 1] + minExonLength + r.nextInt(lengthPool)
                else
                    offset + length)

        val exonRanges = (0 until exonCount).map { Range(exonStarts[it], exonEnds[it]) }.sorted()

        val mRNALength = exonEnds.sum() - exonStarts.sum()
        val cdsLength = if (geneNumber % 5 == 0 || mRNALength < 3) 0 else r.nextInt(mRNALength / 3) * 3

        val cdsRange: Range?
        val utr3End5: Int
        if (cdsLength != 0) {
            var cdsRelativeStartOffset = r.nextInt(mRNALength - cdsLength)
            var cdsStartExon = 0
            while (cdsRelativeStartOffset >= exonEnds[cdsStartExon] - exonStarts[cdsStartExon]) {
                cdsRelativeStartOffset -= exonEnds[cdsStartExon] - exonStarts[cdsStartExon]
                ++cdsStartExon
            }

            val cdsStart = exonStarts[cdsStartExon] + cdsRelativeStartOffset

            while (cdsRelativeStartOffset + cdsLength > exonEnds[cdsStartExon] - exonStarts[cdsStartExon]) {
                cdsRelativeStartOffset -= exonEnds[cdsStartExon] - exonStarts[cdsStartExon]
                ++cdsStartExon
            }

            val cdsEnd = exonStarts[cdsStartExon] + cdsRelativeStartOffset + cdsLength
            cdsRange = Range(cdsStart, cdsEnd)


            utr3End5 = GtfReader.determineUTR3End5(cdsRange.on(chromosome, strand), exonRanges, ensemblId)
        } else {
            cdsRange = null
            utr3End5 = -1
        }
        return Transcript(ensemblId, "ENSG$geneSymbol",
                geneSymbol,
                Location(offset, offset + length, chromosome, strand),
                cdsRange, utr3End5,
                exonRanges)
    }

    /**
     * Generates centromere annotations in UCSC format.
     *
     * @see Gaps
     */
    private fun generateCentromere(genome: Genome) {
        LOG.info("Generating centromere annotations")
        genome.gapsPath.bufferedWriter().use {
            val printer = Gaps.FORMAT.print(it)
            for (chromosome in genome.chromosomes) {
                printer.printRecord(100,
                        chromosome.name,
                        chromosome.length / 2 - 1000,
                        chromosome.length / 2 + 1000,
                        "", "", 2000, "centromere", "")
            }
        }
    }

    /**
     * Generates cytoband annotations in UCSC format.
     *
     * @see CytoBands
     */
    private fun generateCytobands(genome: Genome) {
        LOG.info("Generating cytoband annotations")
        genome.cytobandsPath!!.bufferedWriter().use {
            val printer = CytoBands.FORMAT.print(it)
            for ((i, chromosome) in genome.chromosomes.withIndex()) {
                printer.printRecord(chromosome.name,
                        chromosome.length / 2 - 1000,
                        chromosome.length / 2 + 1000,
                        "q$i.${i % 3 + 1}",
                        "unknown region tag")
            }
        }
    }

    /**
     * Generates (stub !) repeat annotations in UCSC format.
     *
     * @see Repeats
     */
    private fun generateRepeats(genome: Genome) {
        LOG.info("Generating repeat annotations (stub)")
        genome.repeatsPath!!.bufferedWriter().use { w ->
            // This is the first line from mm9 annotations.
            w.write("607\t687\t174\t0\t0\tchr1\t3000001\t3000156\t-194195276\t-" +
                    "\tL1_Mur2\tLINE\tL1\t-4310\t1567\t1413\t1")
        }
    }

    /**
     * Generates (stub !) cpg islands annotations in UCSC format.
     *
     * @see Repeats
     */
    private fun generateCGI(genome: Genome) {
        LOG.info("Generating repeat annotations (stub)")
        genome.cpgIslandsPath!!.bufferedWriter().use { w ->
            // This is the first line from hg19 annotations.
            w.write("""
                |585	chr1	28735	29810	CpG: 116	1075	116	787	21.6	73.2	0.83
                |586	chr1	135124	135563	CpG: 30	439	30	295	13.7	67.2	0.64
            """.trimMargin().trim())
        }
    }

    /**
     * Mapability is a wiggle track with 0 for non-mapable nucleotides and 1 for mapable ones.
     *
     * It's generated in the same folder as "chrom.sizes" with a filename "mapability.bigWig".
     *
     * We only generate mapability for chrX.
     * This is done to test that the genome mean substitution for no-data chromosome works correctly.
     */
    private fun generateMapability(genome: Genome) {
        LOG.info("Generating mapability bigWig")
        val path = genome.chromSizesPath.parent / "mapability.bigWig"
        val gq = genome.toQuery()
        val chrX = gq["chrX"]!!
        val random = ThreadLocalRandom.current()
        val section = FixedStepSection(
            chrX.name,
            start = 0,
            values = TFloatArrayList(
                (0 until chrX.length).map { if (random.nextInt(5) == 4) 0.0f else 1.0f }.toFloatArray()
            )
        )
        BigWigFile.write(listOf(section), gq.get().map { it.name to it.length }, path)
    }
}