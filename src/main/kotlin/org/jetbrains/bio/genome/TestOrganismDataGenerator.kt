package org.jetbrains.bio.genome

import com.google.common.math.IntMath
import gnu.trove.list.array.TFloatArrayList
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.big.FixedStepSection
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.format.FastaRecord
import org.jetbrains.bio.genome.format.TwoBitWriter
import org.jetbrains.bio.genome.format.write
import org.jetbrains.bio.genome.search.sa.SuffixArray
import org.jetbrains.bio.genome.sequence.Nucleotide
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.DataOutputStream
import java.io.IOException
import java.nio.file.Files
import java.nio.file.Path
import java.util.*
import java.util.stream.Collectors
import kotlin.math.abs
import kotlin.math.min

/** A generator for the fake "to1" genome build. */
@Suppress("unused")
object TestOrganismDataGenerator {
    private val LOG = LoggerFactory.getLogger(TestOrganismDataGenerator::class.java)

    val RANDOM_SEED: Long = 42
    private val RANDOM = Random().apply {
        setSeed(RANDOM_SEED)
    }

    private val CHROMOSOMES_SIZES = mapOf(
        "chr1" to IntMath.pow(10, 7),
        "chr2" to IntMath.pow(10, 6),
        "chr3" to IntMath.pow(10, 6),
        "chrX" to IntMath.pow(10, 6),
        "chrM" to IntMath.pow(10, 4)
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

        LOG.info("Check files...")
        val build = Genome.TEST_ORGANISM_BUILD
        val genomePath = Configuration.genomesPath / build
        val chromSizesPath = genomePath / "$build.chrom.sizes"
        val genome = Genome[Genome.TEST_ORGANISM_BUILD]
        val twoBitPath = genome.twoBitPath(downloadIfMissing = false)
        val genesGtfPath = genome.genesGtfPath(downloadIfMissing = false)
        val genesDescriptionPath = genome.genesDescriptionsPath
        val gapsPath = genome.gapsPath!!
        val cytobandsPath = genome.cytobandsPath!!
        val repeatsPath = genome.repeatsPath!!
        val cpgIslandsPath = genome.cpgIslandsPath!!
        val mappabilityPath = genomePath / "mapability.bigWig"
        val saPath = genomePath / "sa"
        if (genomePath.exists && chromSizesPath.exists && twoBitPath.exists &&
            genesGtfPath.exists && genesDescriptionPath.exists && gapsPath.exists && cytobandsPath.exists &&
            repeatsPath.exists && cpgIslandsPath.exists && mappabilityPath.exists && saPath.exists
        ) {
            LOG.info("All files present. Nothing to generate.")
            return
        }

        LOG.info("Generating genome $build files...")
        generateChromSizes(chromSizesPath)
        generateSequence(genome, twoBitPath)
        generateTranscripts(genome, genesGtfPath, genesDescriptionPath)
        generateCentromere(genome, gapsPath)
        generateCytobands(genome, cytobandsPath)
        generateRepeats(repeatsPath)
        generateCGI(cpgIslandsPath)
        generateMapability(genome, mappabilityPath)
        generateSA(genome, saPath)
        LOG.info("Done")
    }

    /**
     * Creates genome.chrom.sizes file.
     */
    private fun generateChromSizes(chromSizesPath: Path) {
        //XXX: no randomization here

        LOG.info("Generating chrom.sizes path")
        CSVFormat.TDF.print(chromSizesPath.bufferedWriter()).use { csvPrinter ->
            CHROMOSOMES_SIZES.forEach { (name, size) ->
                csvPrinter.printRecord(name, size.toString())
            }
        }
    }

    /**
     * Generates 2bit sequences.
     *
     * @param maxGaps maximum number of gaps (sequences of consecutive Ns)
     *                on a chromosome.
     * @param maxGapLength maximum gap length.
     */
    private fun generateSequence(genome: Genome, twoBitPath: Path, maxGaps: Int = 10, maxGapLength: Int = 100) {
        RANDOM.setSeed(RANDOM_SEED) // allows changing several generate* functions w/o affecting each other

        LOG.info("Generating FASTA sequence")
        withTempFile(genome.build, ".fa") { fastaPath ->
            CHROMOSOMES_SIZES.map { (name, length) ->
                val sequence = RANDOM.ints(length.toLong(), 0, Nucleotide.ALPHABET.size)
                    .mapToObj { Nucleotide.ALPHABET[it].toString() }
                    .collect(Collectors.joining())
                    .toCharArray()

                val numGaps = min(RANDOM.nextInt(maxGaps + 1), 3)
                for (i in 0 until numGaps) {
                    // We could do better here by generating non-overlapping
                    // ranges, but it's not worth the effort.
                    val gapLength = RANDOM.nextInt(maxGapLength) + 1
                    val gapOffset = RANDOM.nextInt(length - gapLength)
                    (gapOffset until gapOffset + gapLength).forEach {
                        sequence[it] = 'n'
                    }
                }

                FastaRecord(name, String(sequence))
            }.write(fastaPath)

            LOG.info("Converting FASTA sequence to 2bit")
            TwoBitWriter.convert(fastaPath, twoBitPath)
        }
    }

    /**
     * Generates JSON-serialized transcript annotations.
     *
     * @see Transcripts
     */
    private fun generateTranscripts(genome: Genome, genesGtfPath: Path, genesDescriptionPath: Path) {
        RANDOM.setSeed(RANDOM_SEED) // allows changing several generate* functions w/o affecting each other

        LOG.info("Generating transcript annotations")

        val transcripts = arrayListOf<Transcript>()

        genesGtfPath.bufferedWriter().use { w ->

            for (name in CHROMOSOMES_SIZES.keys.sorted()) {
                val chromosome = Chromosome(genome, name)
                val length = CHROMOSOMES_SIZES[name]!!

                var currentEnd = 0
                var geneNumber = 0
                while (length - currentEnd > 6e4) {
                    val geneSymbol = "simgene.$name.$geneNumber".uppercase()
                    val currentStart = currentEnd + 10000 + RANDOM.nextInt(40000)
                    val currentLength = 100 + min(length - currentStart - 10000, 10000)
                    val strand = if (RANDOM.nextBoolean()) Strand.PLUS else Strand.MINUS
                    val transcript = generateTranscript(
                        geneNumber, geneSymbol, chromosome, strand,
                        currentStart, currentLength
                    )
                    transcripts.add(transcript)

                    currentEnd = currentStart + currentLength
                    geneNumber++
                }
            }

            writeGtf(w, transcripts)
        }

        GeneDescription.serializeFromId2DescriptionMapping(
            transcripts.mapIndexed { idx, t ->
                t.ensemblGeneId to "Test Gene ${t.geneSymbol} [Source:HGNC Symbol;Acc:HGNC:${idx + 1}]"
            },
            genesDescriptionPath
        )
    }

    /** Here be dragons! */
    private fun generateTranscript(
        geneNumber: Int, geneSymbol: String,
        chromosome: Chromosome, strand: Strand,
        offset: Int, length: Int
    ): Transcript {
        val ensemblId = "ENST$geneSymbol"

        val minExonLength = 3
        val minIntronLength = 2
        val exonCount = 1 + RANDOM.nextInt(8)
        val leftFlankingIntron = RANDOM.nextBoolean()
        val rightFlankingIntron = RANDOM.nextBoolean()
        val exonStarts = ArrayList<Int>(exonCount)
        val exonEnds = ArrayList<Int>(exonCount)
        var lengthPool = (length - minExonLength * exonCount - minIntronLength * (exonCount - 1)
                - (if (leftFlankingIntron) minIntronLength else 0)
                - if (rightFlankingIntron) minIntronLength else 0)
        exonStarts.add(0, if (leftFlankingIntron) offset + minIntronLength + RANDOM.nextInt(lengthPool) else offset)
        if (leftFlankingIntron) {
            val intronLen = exonStarts[0] - offset
            lengthPool -= intronLen - minIntronLength // length pool already takes in consideration min lengths
        }
        for (i in 0 until exonCount - 1) {
            exonEnds.add(i, exonStarts[i] + minExonLength + RANDOM.nextInt(lengthPool))
            val exonLen = exonEnds[i] - exonStarts[i]
            lengthPool -= exonLen - minExonLength // length pool already takes in consideration min lengths
            exonStarts.add(i + 1, exonEnds[i] + minIntronLength + RANDOM.nextInt(lengthPool))
            val intronLen = exonStarts[i + 1] - exonEnds[i]
            lengthPool -= intronLen - minIntronLength
        }

        exonEnds.add(
            exonCount - 1,
            if (rightFlankingIntron)
                exonStarts[exonCount - 1] + minExonLength + RANDOM.nextInt(lengthPool)
            else
                offset + length
        )

        val exonRanges = (0 until exonCount).map { Range(exonStarts[it], exonEnds[it]) }.sorted()

        val mRNALength = exonEnds.sum() - exonStarts.sum()
        val cdsLength = if (geneNumber % 5 == 0 || mRNALength < 3) 0 else RANDOM.nextInt(mRNALength / 3) * 3

        val cdsRange: Range?
        val utr3End5: Int
        if (cdsLength != 0) {
            var cdsRelativeStartOffset = RANDOM.nextInt(mRNALength - cdsLength)
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
        return Transcript(
            ensemblId, "ENSG$geneSymbol",
            geneSymbol,
            Location(offset, offset + length, chromosome, strand),
            cdsRange, utr3End5,
            exonRanges
        )
    }

    /**
     * Generates centromere annotations in UCSC format.
     *
     * @see Gaps
     */
    private fun generateCentromere(genome: Genome, gapsPath: Path) {
        //XXX: no randomization here

        LOG.info("Generating centromere annotations")
        gapsPath.bufferedWriter().use {
            val printer = Gaps.FORMAT.print(it)
            for (chromosome in genome.chromosomes) {
                printer.printRecord(
                    100,
                    chromosome.name,
                    chromosome.length / 2 - 1000,
                    chromosome.length / 2 + 1000,
                    "", "", 2000, "centromere", ""
                )
            }
        }
    }

    /**
     * Generates cytoband annotations in UCSC format.
     *
     * @see CytoBands
     */
    private fun generateCytobands(genome: Genome, cytobandsPath: Path) {
        //XXX: no randomization here

        LOG.info("Generating cytoband annotations")
        cytobandsPath.bufferedWriter().use {
            val printer = CytoBands.FORMAT.print(it)
            for ((i, chromosome) in genome.chromosomes.withIndex()) {
                printer.printRecord(
                    chromosome.name,
                    chromosome.length / 2 - 1000,
                    chromosome.length / 2 + 1000,
                    "q$i.${i % 3 + 1}",
                    "unknown region tag"
                )
            }
        }
    }

    /**
     * Generates (stub !) repeat annotations in UCSC format.
     *
     * @see Repeats
     */
    private fun generateRepeats(repeatsPath: Path) {
        //XXX: no randomization here

        LOG.info("Generating repeat annotations (stub)")
        repeatsPath.bufferedWriter().use { w ->
            // This is the first line from mm9 annotations.
            w.write(
                "607\t687\t174\t0\t0\tchr1\t3000001\t3000156\t-194195276\t-" +
                        "\tL1_Mur2\tLINE\tL1\t-4310\t1567\t1413\t1"
            )
        }
    }

    /**
     * Generates (stub !) cpg islands annotations in UCSC format.
     *
     * @see Repeats
     */
    private fun generateCGI(cpgIslandsPath: Path) {
        //XXX: no randomization here

        LOG.info("Generating repeat annotations (stub)")
        cpgIslandsPath.bufferedWriter().use { w ->
            // This is the first line from hg19 annotations.
            w.write(
                """
                |585	chr1	28735	29810	CpG: 116	1075	116	787	21.6	73.2	0.83
                |586	chr1	135124	135563	CpG: 30	439	30	295	13.7	67.2	0.64
            """.trimMargin().trim()
            )
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
    private fun generateMapability(genome: Genome, mappabilityPath: Path) {
        RANDOM.setSeed(RANDOM_SEED) // allows changing several generate* functions w/o affecting each other

        LOG.info("Generating mapability bigWig")
        val gq = genome.toQuery()
        val chrX = gq["chrX"]!!
        val section = FixedStepSection(
            chrX.name,
            start = 0,
            values = TFloatArrayList(
                (0 until chrX.length).map { if (RANDOM.nextInt(5) == 4) 0.0f else 1.0f }.toFloatArray()
            )
        )
        BigWigFile.write(listOf(section), gq.get().map { it.name to it.length }, mappabilityPath)
    }

    private fun generateSA(genome: Genome, saPath: Path) {
        RANDOM.setSeed(RANDOM_SEED) // allows changing several generate* functions w/o affecting each other

        LOG.info("Processing SA indexes and FASTQ mismatched reads")
        for (chromosome in genome.chromosomes) {
            SuffixArray.create(chromosome)
            val sequence = chromosome.sequence
            // Generate 35-bp FASTQ reads with 0 to 2 mismatches
            DataOutputStream(
                Files.newOutputStream(saPath / "${chromosome.name}.fastq")
            ).use { outputStream ->
                var j = 0
                while (j < 1e3) {
                    writeFastq(outputStream, j, sequence, 35, j % 3)
                    j++
                }
                outputStream.flush()
            }
        }
    }

    @Suppress("SameParameterValue")
    @Throws(IOException::class)
    private fun writeFastq(
        file: DataOutputStream,
        number: Int,
        sequence: NucleotideSequence,
        length: Int,
        maxMismatchCount: Int
    ) {
        val strand = if (RANDOM.nextBoolean()) Strand.PLUS else Strand.MINUS
        var pos: Int
        var read: String
        do {
            pos = RANDOM.nextInt(sequence.length - length)
            read = sequence.substring(pos, pos + length, strand)
        } while (read.indexOf('n') != -1)

        val alphabet = Nucleotide.ALPHABET
        val readCharArray = read.toCharArray()
        val mismatches = RANDOM.ints().map { i -> abs(i % length) }.distinct().limit(maxMismatchCount.toLong())
        mismatches.forEach { mismatchPos -> readCharArray[mismatchPos] = alphabet[RANDOM.nextInt(alphabet.size)] }

        read = String(readCharArray).uppercase()
        val header = "TESTDATA.$number TAKENFROM:$pos:$strand:$maxMismatchCount length=$length\n"
        file.writeBytes("@$header")
        file.writeBytes(read + '\n')
        file.writeBytes("+$header")
        for (i in 0 until length) {
            file.writeByte('B'.code)
        }
        file.writeByte('\n'.code)
    }
}