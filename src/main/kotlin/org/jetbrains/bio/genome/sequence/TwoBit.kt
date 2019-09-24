// Remove this once KT-11820 is done.
@file:Suppress("UsePropertyAccessSyntax")

package org.jetbrains.bio.genome.sequence

import com.google.common.base.Preconditions.checkElementIndex
import com.google.common.collect.ImmutableSet
import com.google.common.math.IntMath
import com.google.common.primitives.Ints
import gnu.trove.list.array.TIntArrayList
import gnu.trove.map.TObjectIntMap
import gnu.trove.map.hash.TObjectIntHashMap
import org.jetbrains.bio.io.FastaReader
import org.jetbrains.bio.util.size
import java.io.DataOutput
import java.io.DataOutputStream
import java.io.IOException
import java.io.RandomAccessFile
import java.math.RoundingMode
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.nio.channels.Channels
import java.nio.channels.FileChannel
import java.nio.file.Path
import java.util.*
import java.util.stream.Collectors

/**
 * Compressed DNA sequence.
 *
 * The compression works as follows:
 *
 *   1. Each of the four nucleotides is packed to 2 bits
 *      (00 - T, 01 - C, 10 - A, 11 - G).
 *   2. N-blocks are masked by the nucleotide preceding the block, for
 *      example "ACNNNG" will be masked as "ACCCCG". Thus, the length
 *      of the compressed sequence _exactly_ matches the length of the
 *      uncompressed sequence. Unfortunately UCSC spec. doesn't stress
 *      this enough.
 *   3. The coordinates of the N-blocks are extracted into two arrays:
 *      block starts and block sizes.
 *
 * The format also allows to include repeat masks, but our implementation
 * doesn't use them at the moment.
 *
 * @author Sergei Lebedev
 */
class TwoBitSequence(
        /** Total length of the DNA sequence including Ns. */
        override val length: Int,
        /** N-block coordinates. */
        private val nBlockStarts: IntArray, private val nBlockSizes: IntArray,
        /** 2-bit packed DNA sequence. */
        private val packedDna: IntArray) : NucleotideSequence {

    private val nBlockLut = BinaryLut.of(nBlockStarts, 8)

    init {
        require(length > 0) { "invalid length" }
        require(nBlockStarts.size == nBlockSizes.size) { "missing N blocks" }
    }

    override fun charAt(pos: Int): Char {
        checkElementIndex(pos, length, "index")
        val block = getBlock(pos)
        if (inNBlock(pos, block)) {
            return Nucleotide.N
        }

        return Nucleotide.getChar(getByte(pos))
    }

    private fun inNBlock(index: Int, block: Int): Boolean {
        return block >= 0 && block < nBlockStarts.size &&
               index - nBlockStarts[block] < nBlockSizes[block]
    }

    private fun getByte(index: Int): Byte {
        val i = index / NUCLEOTIDE_PER_INTEGER  // offset in packed DNA array.
        // Bits in a pack are ordered rtl, while DNA is packed ltr,
        // for example, "TCAG" is packed as "00011011". Thus we have to wrap
        // around the bit index.
        val j = NUCLEOTIDE_PER_INTEGER - 1 - index % NUCLEOTIDE_PER_INTEGER
        return (packedDna[i] and NUCLEOTIDE_MASKS[j] ushr BITS_PER_NUCLEOTIDE * j).toByte()
    }

    /**
     * Returns the index of the block containing a given position or -1 otherwise.
     */
    private fun getBlock(pos: Int): Int {
        if (nBlockStarts.size == 0 || pos < nBlockStarts[0]) {
            return -1
        }

        val i = nBlockLut.binarySearch(nBlockStarts, pos)
        return if (i < 0) i.inv() - 1 else i  // return preceding block if not found.
    }

    @Throws(IOException::class)
    fun write(output: DataOutput) = with(output) {
        writeInt(length)

        writeInt(nBlockStarts.size)
        for (nBlockStart in nBlockStarts) {
            writeInt(nBlockStart)
        }
        for (nBlockSize in nBlockSizes) {
            writeInt(nBlockSize)
        }

        writeInt(0)  // maskBlockCount.
        writeInt(0)  // reserved.
        for (pack in packedDna) {
            writeInt(pack)
        }
    }

    override fun equals(other: Any?): Boolean = when {
        this === other ->  true
        other == null || javaClass != other.javaClass -> false
        else -> {
            val tbs = other as TwoBitSequence
            length == tbs.length &&
            Arrays.equals(nBlockSizes, tbs.nBlockSizes) &&
            Arrays.equals(nBlockStarts, tbs.nBlockStarts) &&
            Arrays.equals(packedDna, tbs.packedDna)
        }
    }

    override fun hashCode() = Arrays.deepHashCode(
            arrayOf(length, nBlockStarts, nBlockSizes, packedDna))

    override fun toString() = substring(0, length)

    companion object {
        /** Number of bits per nucleotide.  */
        val BITS_PER_NUCLEOTIDE = 2

        /** Nucleotides to be encoded per `int` value.  */
        val NUCLEOTIDE_PER_INTEGER = Integer.SIZE / BITS_PER_NUCLEOTIDE

        /**
         * An array of nucleotide masks.
         *
         * @see getByte
         */
        private val NUCLEOTIDE_MASKS = IntArray(NUCLEOTIDE_PER_INTEGER)

        init {
            NUCLEOTIDE_MASKS[0] = (1 shl BITS_PER_NUCLEOTIDE) - 1
            for (i in 1..NUCLEOTIDE_PER_INTEGER - 1) {
                NUCLEOTIDE_MASKS[i] = NUCLEOTIDE_MASKS[i - 1] shl BITS_PER_NUCLEOTIDE
            }
        }

        private fun String.pack(start: Int, end: Int): Int {
            val l = end - start
            require(l <= NUCLEOTIDE_PER_INTEGER) {
                "Cannot pack $l nucleotides into ${Integer.SIZE} bits."
            }

            val mask = (1 shl BITS_PER_NUCLEOTIDE) - 1
            var acc = 0
            for (i in start..end - 1) {
                acc = acc shl BITS_PER_NUCLEOTIDE
                acc += Nucleotide.getByte(this[i]).toInt() and mask
            }

            // Ensure acc is correctly aligned.
            return acc shl (BITS_PER_NUCLEOTIDE * (NUCLEOTIDE_PER_INTEGER - l))
        }

        /**
         * Encodes a human-readable genomic sequence in 2bit format.
         */
        fun encode(s: String): TwoBitSequence {
            val dnaSize = s.length
            require(dnaSize > 0) { "no data" }

            val nBlockStarts = TIntArrayList()
            val nBlockSizes = TIntArrayList()
            var start = -1
            for (i in 0..dnaSize - 1) {
                if (Nucleotide.getByte(s[i]) == Nucleotide.ANY_NUCLEOTIDE_BYTE) {
                    if (start < 0) {
                        start = i
                    }
                } else if (start >= 0) {
                    nBlockStarts.add(start)
                    nBlockSizes.add(i - start)
                    start = -1
                }
            }

            if (start >= 0) {
                nBlockStarts.add(start)
                nBlockSizes.add(dnaSize - start)
            }

            val packedDna = TIntArrayList(dnaSize / NUCLEOTIDE_PER_INTEGER)
            var i = 0
            while (i < dnaSize) {
                val j = Math.min(i + NUCLEOTIDE_PER_INTEGER, dnaSize)
                packedDna.add(s.pack(i, j))
                i += NUCLEOTIDE_PER_INTEGER
            }

            return TwoBitSequence(dnaSize, nBlockStarts.toArray(), nBlockSizes.toArray(),
                                  packedDna.toArray())
        }
    }
}

/**
 * A reader for 2bit compressed sequence format.
 *
 * See http://genome.ucsc.edu/FAQ/FAQformat.html.format7
 *
 * @author Sergei Lebedev
 */
object TwoBitReader {
    /** Magic number used for determining [ByteOrder].  */
    internal const val MAGIC = 0x1A412743

    /**
     * Extracts available sequence names from the 2bit file specified
     * by `path`.
     *
     * @throws IllegalStateException if the 2bit file at `path`
     *         doesn't conform to 2bit format specification.
     */
    @Throws(IOException::class)
    @JvmStatic fun names(path: Path): ImmutableSet<String> {
        val buf = getBuffer(path)
        val index = getIndex(buf)
        return ImmutableSet.copyOf(index.keySet())
    }

    /**
     * Extracts length of the sequence corresponding to the given
     * `name` from the 2bit file specified by `path`.
     *
     * @throws IllegalStateException if the 2bit file at `path`
     *         doesn't conform to 2bit format specification.
     * @throws IllegalArgumentException if the 2bit index doesn't
     *         contain `name`.
     */
    @Throws(IOException::class)
    @JvmStatic fun length(path: Path, name: String): Int {
        val buf = getBuffer(path)
        val index = getIndex(buf)
        require(index.containsKey(name)) { "unknown sequence: $name" }
        buf.position(index[name])
        return buf.getInt()
    }

    /**
     * Reads a sequence corresponding to the given `name` from
     * the 2bit file specified by `path`.
     *
     * @throws IllegalStateException if the 2bit file at `path`
     *         doesn't conform to 2bit format specification.
     * @throws IllegalArgumentException if the 2bit index doesn't
     *         contain `name`.
     */
    @Throws(IOException::class)
    @JvmStatic fun read(path: Path, name: String): TwoBitSequence {
        val buf = getBuffer(path)
        val index = getIndex(buf)
        require(index.containsKey(name)) { "unknown sequence: $name" }
        buf.position(index[name])
        return read(buf)
    }

    /**
     * Reads a mapping from sequence names to offsets in the 2bit file.
     */
    internal fun getIndex(buf: ByteBuffer): TObjectIntMap<String> {
        val version = buf.getInt()
        val sequenceCount = buf.getInt()
        val reserved = buf.getInt()

        check(reserved == 0) { "invalid reserved value: $reserved" }
        check(version == 0) { "unexpected version: $version" }
        check(sequenceCount > 0) { "invalid number of sequences: $sequenceCount" }

        val index = TObjectIntHashMap<String>(sequenceCount)
        for (i in 0..sequenceCount - 1) {
            val chunk = ByteArray(buf.get().toInt())
            buf.get(chunk)
            index.put(String(chunk), buf.getInt())
        }

        return index
    }

    private fun ByteOrder.flip() = when (this) {
        ByteOrder.BIG_ENDIAN -> ByteOrder.LITTLE_ENDIAN
        ByteOrder.LITTLE_ENDIAN -> ByteOrder.BIG_ENDIAN
        else -> throw IllegalStateException()
    }

    /**
     * Determines the byte order in the 2bit file at `path`.
     */
    private fun getBuffer(path: Path): ByteBuffer {
        return FileChannel.open(path).use { fc ->
            val buf = fc.map(FileChannel.MapMode.READ_ONLY, 0, path.size.bytes)
            buf.order(ByteOrder.nativeOrder())
            val nativeMagic = buf.getInt()
            if (nativeMagic != MAGIC) {
                val reversedMagic = Integer.reverseBytes(nativeMagic)
                check(reversedMagic == MAGIC) { "bad signature" }
                buf.order(buf.order().flip())
            }

            buf
        }
    }

    fun read(buf: ByteBuffer): TwoBitSequence {
        val view = buf.asIntBuffer()
        val dnaSize = view.get()
        val nBlockCount = view.get()
        val nBlockStarts = IntArray(nBlockCount)
        view.get(nBlockStarts)

        val nBlockSizes = IntArray(nBlockCount)
        view.get(nBlockSizes)

        // skip repeat masks.
        val maskBlockCount = view.get()
        view.position(view.position() + maskBlockCount + maskBlockCount)

        val reserved = view.get()
        check(reserved == 0) { "invalid reserved value: $reserved" }

        val packsCount = dnaSize / TwoBitSequence.NUCLEOTIDE_PER_INTEGER

        // 'dnaSize' might _not_ be a multiple of 'Integer.SIZE', so we
        // have to make sure we read any leftover bases.
        val leftover = dnaSize % TwoBitSequence.NUCLEOTIDE_PER_INTEGER
        val packedDna = IntArray(packsCount + (if (leftover > 0) 1 else 0))
        // we're forced to use 'buf' here because DNA is laid out independent
        // of the byte order used.
        buf.position(buf.position() + view.position() * Integer.BYTES)
        buf.order(ByteOrder.BIG_ENDIAN).asIntBuffer().get(packedDna, 0, packsCount)

        if (leftover > 0) {
            buf.position(buf.position() + packsCount * Integer.BYTES)
            val pack = ByteArray(Integer.BYTES)
            buf.get(pack, 0, IntMath.divide(leftover * TwoBitSequence.BITS_PER_NUCLEOTIDE,
                                            java.lang.Byte.SIZE, RoundingMode.CEILING))
            packedDna[packsCount] = Ints.fromByteArray(pack)
        }

        buf.clear()
        return TwoBitSequence(dnaSize, nBlockStarts, nBlockSizes, packedDna)
    }
}

/**
 * An unfortunately named class for converting FASTA to 2bit.
 *
 * @author Sergei Lebedev
 */
object TwoBitWriter {
    @Throws(IOException::class)
    @JvmStatic fun convert(fastaPath: Path, twoBitPath: Path) {
        val chrNames = FastaReader.read(fastaPath, false)
            .map { it.description.split(' ').first() }
            .collect(Collectors.toList())
        convert(chrNames, fastaPath, twoBitPath)
    }

    @Throws(IOException::class)
    fun convert(chrNames: List<String>, fastaPath: Path, twoBitPath: Path) {
        RandomAccessFile(twoBitPath.toFile(), "rw").use { raf ->
            // XXX this always uses BE byte order.
            with(raf) {
                writeInt(TwoBitReader.MAGIC)
                writeInt(0)  // version.
                writeInt(chrNames.size)
                writeInt(0)  // reserved.

                val index = IntArray(chrNames.size)
                chrNames.forEachIndexed { i, chrName ->
                    write(chrName.length)
                    writeBytes(chrName)
                    index[i] = Ints.checkedCast(filePointer)
                    writeInt(0)
                }

                var recIdx = 0
                FastaReader.read(fastaPath, true).forEach { record ->
                    val offset = Ints.checkedCast(filePointer)
                    seek(index[recIdx].toLong())
                    writeInt(offset)
                    seek(offset.toLong())

                    val outputStream = Channels.newOutputStream(raf.channel).buffered()
                    TwoBitSequence.encode(record.sequence).write(DataOutputStream(outputStream))
                    outputStream.flush()

                    recIdx++
                }
            }
        }
    }
}