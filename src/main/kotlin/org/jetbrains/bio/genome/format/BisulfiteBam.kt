package org.jetbrains.bio.genome.format

import com.google.common.collect.Iterators
import gnu.trove.list.array.TByteArrayList
import gnu.trove.list.array.TCharArrayList
import htsjdk.samtools.*
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.cram.ref.ReferenceSource
import kotlinx.support.jdk7.use
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.methylome.CytosineContext
import org.jetbrains.bio.genome.methylome.Methylome
import org.jetbrains.bio.genome.methylome.MethylomeBuilder
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import org.jetbrains.bio.genome.sequence.asNucleotideSequence
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.awaitAll
import org.jetbrains.bio.util.name
import org.jetbrains.bio.util.parallelismLevel
import org.slf4j.LoggerFactory
import java.io.File
import java.nio.file.Path
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.Executors

/**
 * A parser for raw BS-seq alignments in BAM or CRAM formats.
 *
 * The only counts the reads for the strand with a cytosine. For example
 * the following
 *
 *   +5C
 *   |
 *   ACAGATCGG
 *   |
 *   -2C
 *
 * results in 5C count for the position being considered.
 *
 * @author Sergei Lebedev
 */
object BisulfiteBamParser {
    fun parse(path: Path, genomeQuery: GenomeQuery, minBasePhred: Byte): Methylome {
        // This is of course an upper bound, but it's better than nothing.
        val total = genomeQuery.get().parallelStream().mapToLong {
            var acc = 0L
            for (i in 0 until it.length) {
                val b = it.sequence.charAt(i)
                if (b == 'c' || b == 'g') {
                    acc++
                }
            }

            acc
        }.sum()

        val builder = Methylome.builder(genomeQuery)
        val progress = Progress {
            title = path.name
        }.bounded(total)

        val executor = Executors.newWorkStealingPool(parallelismLevel())
        executor.awaitAll(genomeQuery.get().map {
            Callable {
                val sequence = it.sequence.toString().asNucleotideSequence()
                parse(path, it, sequence, builder, progress, minBasePhred)
            }
        })
        check(executor.shutdownNow().isEmpty())

        progress.done()
        return builder.build()
    }

    fun parse(
        path: Path, chromosome: Chromosome, sequence: NucleotideSequence,
        builder: MethylomeBuilder, progress: Progress,
        minBasePhred: Byte
    ) {

        getPiler(path, chromosome, minBasePhred).use { piler ->
            loop@ while (piler.hasNext()) {
                val pc = piler.next()
                // Reference position is 1-based in SAM.
                val offset = pc.position - 1
                val ch = sequence.charAt(offset)
                val strand = when (ch) {
                    // Skip records corresponding to non-cytosines.
                    'a', 't', 'n' -> continue@loop
                    'c' -> Strand.PLUS
                    'g' -> Strand.MINUS
                    else -> error("unexpected character in sequence: $ch")
                }

                val context = CytosineContext.determine(sequence, offset, strand)
                piler.handleLocus(builder, chromosome, strand, offset, context, pc)
                progress.report()
            }
        }
    }

    private fun getPiler(path: Path, chromosome: Chromosome, minBasePhred: Byte): SamPiler<out PilerColumn> {
        // 'htsjdk' doesn't allow concurrent queries on 'BAMFileReader'
        // thus we have to re-create 'SamReader' for each chromosome.
        val samReader = SamReaderFactory.makeDefault()
            .referenceSource(
                WrappedReferenceSource(
                    chromosome.name,
                    chromosome.sequence
                )
            )
            .open(path.toFile())

        check(samReader.hasIndex()) { "$path must have index available" }
        return BowtieSamPiler(samReader, chromosome, minBasePhred)
    }
}

/**
 * Piler piles up SAM and BAM records.
 *
 * See http://thefreedictionary.com/piler if you don't believe the KDoc.
 *
 * @author Sergei Lebedev
 */
private abstract class SamPiler<T : PilerColumn>(
    private val samReader: SamReader,
    chromosome: Chromosome,
    private val bismarkFormat: Boolean = false
) :
    Iterator<T>, AutoCloseable {

    /** An iterator for SAM records with the SAME reference index. */
    private val samIterator = Iterators.peekingIterator<SAMRecord>(
        samReader.query(chromosome.name, 0, 0, false)
    )

    /** A queue for loci waiting for more SAM records coming.  */
    protected val waitingQueue = ArrayList<T>()

    /** A queue for completed loci, which won't get any more SAM records.  */
    private val completedQueue = ArrayDeque<T>()

    init {
        val samHeader = samReader.fileHeader
        if (samHeader.sequenceDictionary.getSequence(chromosome.name) == null) {
            throw NoSuchElementException(chromosome.name)
        }

        val sortOrder = samHeader.sortOrder
        if (sortOrder == null || sortOrder == SortOrder.unsorted) {
            LOG.warn(
                "SAM sort order is unspecified. "
                        + "Assuming SAM is coordinate sorted, but exceptions may occur if it is not."
            )
        } else if (sortOrder != SortOrder.coordinate) {
            throw IllegalArgumentException(
                "Cannot operate on a SAM file that is not coordinate sorted."
            )
        }
    }

    protected abstract fun updateWaiting(record: SAMRecord)
    abstract fun handleLocus(
        builder: MethylomeBuilder,
        chromosome: Chromosome, strand: Strand,
        offset: Int, context: CytosineContext?,
        pc: PilerColumn
    )

    override fun hasNext(): Boolean {
        prefetch()
        return completedQueue.isNotEmpty()
    }

    override fun next(): T {
        check(completedQueue.isNotEmpty()) { "no data" }
        return completedQueue.removeFirst()
    }

    override fun close() = samReader.close()

    private fun prefetch() {
        while (completedQueue.isEmpty() && samIterator.hasNext()) {
            val record = samIterator.peek()
            if (record.readUnmappedFlag
                || record.isSecondaryOrSupplementary
                || record.duplicateReadFlag
                || record.mappingQuality == 0
            ) {  // BWA multi-alignment evidence
                samIterator.next()
                continue
            }

            completeWaiting(record)
            samIterator.next()
            updateWaiting(record)
        }

        if (completedQueue.isEmpty() && !samIterator.hasNext()) {
            while (waitingQueue.isNotEmpty()) {
                val pc = waitingQueue.removeAt(0)
                if (pc.isNotEmpty()) {
                    completedQueue.addLast(pc)
                }
            }
        }
    }

    private fun completeWaiting(record: SAMRecord) {
        // Complete piled up columns preceding alignment start.
        val alignmentStart = record.alignmentStart
        while (waitingQueue.isNotEmpty() && waitingQueue[0].position < alignmentStart) {
            val pc = waitingQueue.removeAt(0)
            if (pc.isNotEmpty()) {
                completedQueue.addLast(pc)
            }
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(SamPiler::class.java)
    }
}

private class BowtieSamPiler(
    samReader: SamReader,
    chromosome: Chromosome,
    private val minBasePhred: Byte = 0
) : SamPiler<BowtiePilerColumn>(samReader, chromosome) {

    override fun updateWaiting(record: SAMRecord) {
        val cigar = record.cigar ?: return

        val alignmentStart = record.alignmentStart
        var readBase = 1
        var refBase = alignmentStart
        for (c in 0 until cigar.numCigarElements()) {
            val e = cigar.getCigarElement(c)
            val length = e.length
            when (e.operator!!) {
                CigarOperator.H, CigarOperator.P -> {
                }  // ignore hard clips and pads
                CigarOperator.S -> readBase += length   // soft clip read bases
                CigarOperator.N -> refBase += length    // reference skip
                CigarOperator.D -> refBase += length
                CigarOperator.I -> readBase += length
                CigarOperator.M, CigarOperator.EQ, CigarOperator.X -> {
                    for (i in 0 until length) {
                        // 1-based reference position that the current base aligns to
                        val refPos = refBase + i

                        // 0-based offset from the aligned position of the first base in
                        // the read to the aligned position of the current base.
                        val refOffset = refPos - alignmentStart

                        // Ensure there are columns up to and including this position.
                        waitingQueue.ensureCapacity(refOffset - waitingQueue.size)
                        for (j in waitingQueue.size..refOffset) {
                            waitingQueue.add(BowtiePilerColumn(alignmentStart + j))
                        }

                        // 0-based offset into the read of the current base
                        val readOffset = readBase + i - 1
                        if (record.baseQualities[readOffset] >= minBasePhred) {
                            waitingQueue[refOffset].add(readOffset, record, "")
                        }
                    }

                    readBase += length
                    refBase += length
                }
            }
        }
    }

    override fun handleLocus(
        builder: MethylomeBuilder,
        chromosome: Chromosome, strand: Strand,
        offset: Int, context: CytosineContext?,
        pc: PilerColumn
    ) {
        var countA = 0
        var countT = 0
        var countC = 0
        var countG = 0

        if (strand.isPlus()) {
            for (i in 0 until pc.size()) {
                // The strand of the record should match the reference strand,
                // because otherwise we aren't looking at a cytosine.
                if (pc.getStrand(i).isMinus()) {
                    continue
                }

                when (pc.getReadBase(i).toChar()) {
                    'A' -> countA++
                    'T' -> countT++
                    '=',  // fall-through.
                    'C' -> countC++
                    'G' -> countG++
                }
            }
        } else {
            for (i in 0 until pc.size()) {
                if (pc.getStrand(i).isPlus()) {
                    continue
                }

                when (pc.getReadBase(i).toChar()) {
                    'A' -> countT++
                    'T' -> countA++
                    'C' -> countG++
                    '=',  // fall-through.
                    'G' -> countC++
                }
            }
        }

        val countATCG = countA + countT + countC + countG
        // XXX for the minus strand we add up the counts using
        // reverse-complementary nucleotides, so to get 'methylatedCount'
        // we need to use 'countC' and NOT 'countG'.
        builder.add(chromosome, strand, offset, context, countC, countATCG)
    }

}

private class BismarkSamPiler(
    samReader: SamReader,
    chromosome: Chromosome,
    private val minBasePhred: Byte = 0
) : SamPiler<BismarkPilerColumn>(samReader, chromosome) {

    override fun updateWaiting(record: SAMRecord) {
        val alignmentStart = record.alignmentStart

        val bismarkAlignment = record.getAttribute("XM") as String
        for (readOffset in bismarkAlignment.indices) {
            // 1-based reference position that the current base aligns to
            val refPos = alignmentStart + readOffset

            // 0-based offset from the aligned position of the first base in
            // the read to the aligned position of the current base.
            val refOffset = refPos - alignmentStart

            // Ensure there are columns up to and including this position.
            waitingQueue.ensureCapacity(refOffset - waitingQueue.size)
            for (j in waitingQueue.size..refOffset) {
                waitingQueue.add(BismarkPilerColumn(alignmentStart + j))
            }

            // 0-based offset into the read of the current base
            if (record.baseQualities[readOffset] >= minBasePhred) {
                waitingQueue[refOffset].add(readOffset, record, bismarkAlignment)
            }
        }
    }

    override fun handleLocus(
        builder: MethylomeBuilder,
        chromosome: Chromosome, strand: Strand,
        offset: Int, context: CytosineContext?,
        pc: PilerColumn
    ) {
        var countA = 0
        var countT = 0
        var countC = 0
        var countG = 0

        for (i in 0 until pc.size()) {
            if (pc.getStrand(i) != strand) {
                continue
            }
            val baseAnn = pc.getBaseAnnotation(i)
            when (baseAnn) {
                'Z', 'X', 'H', 'U' -> countC++ //methylated
                'z', 'x', 'h', 'u' -> countT++ //unmethylated
            }
            if (baseAnn != '.') {
                val bisContext = when (baseAnn) {
                    'Z', 'z' -> CytosineContext.CG
                    'X', 'x' -> CytosineContext.CHG
                    'H', 'h' -> CytosineContext.CHH
                    'U', 'u' -> CytosineContext.ANY
                    else -> error("unexpected character in annotation: $baseAnn")
                }
                if (bisContext != context) {
                    println("${strand.char}${chromosome.name}:$offset bismark $bisContext != $context")
                }
            }

        }

        val countATCG = countA + countT + countC + countG
        // XXX for the minus strand we add up the counts using
        // reverse-complementary nucleotides, so to get 'methylatedCount'
        // we need to use 'countC' and NOT 'countG'.
        builder.add(chromosome, strand, offset, context, countC, countATCG)
    }

}

/**
 * A single piled up column. We deliberately store _only_ the information
 * required for [BisulfiteBamParser] to reduce allocation rate.
 *
 * @author Sergei Lebedev
 */
interface PilerColumn {
    val position: Int
    fun add(offset: Int, record: SAMRecord, alignmentStr: String)
    fun getStrand(i: Int): Strand
    fun getReadBase(i: Int): Byte
    fun getBaseAnnotation(i: Int): Char
    fun size(): Int
    fun isNotEmpty(): Boolean
}

private data class BowtiePilerColumn(override val position: Int) : PilerColumn {
    private val bases = TByteArrayList()

    override fun add(offset: Int, record: SAMRecord, alignmentStr: String) {
        val base = record.readBases[offset]
        assert(base > 0) { "read base should be positive." }
        bases.add(if (record.readNegativeStrandFlag) (-base).toByte() else base)
    }

    override fun getStrand(i: Int): Strand {
        return if (bases[i] < 0) Strand.MINUS else Strand.PLUS
    }

    override fun getReadBase(i: Int): Byte {
        val base = bases[i]
        return if (base < 0) (-base).toByte() else base
    }

    override fun getBaseAnnotation(i: Int) = getReadBase(i).toChar()

    override fun size() = bases.size()

    override fun isNotEmpty() = !bases.isEmpty
}

private data class BismarkPilerColumn(override val position: Int) : PilerColumn {
    private val bases = TByteArrayList()
    private val annotations = TCharArrayList()
    private val readNegativeStrandFlags = BitSet()

    override fun add(offset: Int, record: SAMRecord, bismarkBsAnnotation: String) {
        val base = record.readBases[offset]
        bases.add(if (record.readNegativeStrandFlag) (-base).toByte() else base)

        // todo: here we can ignore '.' in bismark
        val baseAnn: Char
        if (record.readNegativeStrandFlag) {
            val currPosition = annotations.size()
            readNegativeStrandFlags.set(currPosition)
            baseAnn = bismarkBsAnnotation[bismarkBsAnnotation.length - 1 - offset]
        } else {
            baseAnn = bismarkBsAnnotation[offset]
        }
        annotations.add(baseAnn)
    }

    override fun getStrand(i: Int) = if (readNegativeStrandFlags.get(i)) Strand.MINUS else Strand.PLUS

    override fun getReadBase(i: Int): Byte {
        val base = bases[i]
        return if (base < 0) (-base).toByte() else base
    }

    override fun getBaseAnnotation(i: Int) = annotations[i]

    override fun size() = annotations.size()

    override fun isNotEmpty() = !annotations.isEmpty
}

/**
 * A fake [htsjdk.samtools.cram.ref.ReferenceSource] which stores bytes for a single [Chromosome].
 */
internal class WrappedReferenceSource(
    private val name: String,
    sequence: NucleotideSequence
) : ReferenceSource(null as File?) {

    // XXX please keep lazy to reduce memory consumption.
    private val bytes: ByteArray by lazy(LazyThreadSafetyMode.PUBLICATION) {
        sequence.toString().toByteArray()
    }

    override fun getReferenceBases(
        record: SAMSequenceRecord,
        tryNameVariants: Boolean
    ): ByteArray? {
        if (tryNameVariants) {
            for (variant in getVariants(record.sequenceName)) {
                if (variant == name) {
                    return bytes
                }
            }
        } else {
            if (record.sequenceName == name) {
                return bytes
            }
        }

        return null
    }

    private fun getVariants(name: String): List<String> {
        val variants = ArrayList<String>()
        when {
            name == "M" -> variants.add("MT")
            name == "MT" -> variants.add("M")
            name.startsWith("chr") -> variants.add(name.substring(3))
            else -> variants.add("chr$name")
        }

        if (name == "chrM") {
            // chrM case:
            variants.add("MT")
        }

        return variants
    }
}