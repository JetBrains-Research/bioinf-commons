package org.jetbrains.bio.io

import com.google.common.base.Splitter
import com.google.common.collect.UnmodifiableIterator
import com.google.common.io.Closeables
import kotlinx.support.jdk7.use
import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.big.ExtendedBedEntry
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.toStrand
import org.jetbrains.bio.io.BedFormat.Companion.auto
import org.jetbrains.bio.io.BedParser.Companion.Stringency
import org.jetbrains.bio.io.BedParser.Companion.Stringency.LENIENT
import org.jetbrains.bio.io.BedParser.Companion.Stringency.STRICT
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.awt.Color
import java.io.*
import java.net.URI
import java.nio.file.Path


private val NON_DATA_LINE_PATTERN = "^(?:#|track|browser).*$".toRegex()

/**
 * Even though BED has a well-defined spec it is widely known
 * as a very non-deterministic format. [BedFormat] is here to
 * help you fight the chaos, allowing for various subsets
 * of fields (aka columns) to be included or excluded.
 *
 * We consider all BED files as having a bedN+ format, where 3 <= N <= 15. Thus, we don't pose any limit
 * on the number (and existence) of the extra fields.
 *
 * If you have a good BED file, just use [auto].
 *
 * See https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for
 * the complete spec.
 *
 * @author Oleg Shpynov
 * @author Sergei Lebedev
 */
data class BedFormat(
        val fieldsNumber: Byte = 6,
        val delimiter: Char = '\t'
) {

    val fmtStr: String

    /**
     * Create Format containing standard BedEntry fields up to given one
     */
    constructor(field: BedField)
            : this((field.field.index + 1).toByte())

    init {
        require(fieldsNumber in 3..15) {
            "Fields number in BED file is between 3 and 15, but was $fieldsNumber"
        }

        fmtStr = "bed$fieldsNumber+"
    }

    fun <T> parse(reader: Reader, src: String = "unknown source", f: (BedParser) -> T) = parse(
            reader.buffered(), src, f
    )

    fun <T> parse(path: Path, f: (BedParser) -> T) = parse(
            path.bufferedReader(), path.toAbsolutePath().toString(), f
    )

    fun parseLocations(path: Path, genome: Genome) = parse(
            path.bufferedReader(),
            path.toAbsolutePath().toString()) {
        it.map { entry ->
            entry.unpack(
                minOf(fieldsNumber, (BedField.STRAND.field.index + 1).toByte()), false, delimiter
            ).toLocation(genome)
        }
    }

    private fun <T> parse(reader: BufferedReader, src: String, f: (BedParser) -> T) =
            BedParser(reader, src, delimiter).use { f(it) }

    fun print(path: Path) = print(path.bufferedWriter())

    fun print(writer: Writer) = BedPrinter(writer.buffered(), this)

    /**
     * Returns a copy of the format with a given delimiter.
     */
    fun delimiter(delimiter: Char) = copy(delimiter = delimiter)

    operator fun contains(field: BedField): Boolean = fieldsNumber >= (field.field.index + 1)

    override fun toString() = "($fmtStr, '$delimiter')"

    companion object {
        val DEFAULT = from("bed12")
        /**
         * chromosome, start, end, name, score, strand, coverage/foldchange, -log(pvalue), -log(qvalue)
         */
        @Suppress("unused")
        val MACS2 = from("bed6+3")
        val RGB = from("bed9")

        fun from(fmtStr: String, delimiter: Char = '\t') =
                BedFormat(parseFormatString(fmtStr), delimiter = delimiter)

        fun fromString(s: String): BedFormat {
            if (s.first() == '(' && s.last() == ')') {
                val sep = s.indexOf(',')
                if (sep != -1) {
                    val fmtStr = s.substring(1, sep)
                    val delimiterStr = s.substring(sep + 1, s.length - 1).trim()
                    if (delimiterStr.length == 3 && delimiterStr.first() == '\'' && delimiterStr.last() == '\'' ) {
                        return from(fmtStr, delimiterStr[1])
                    }
                }
            }
            error("Unsupported string: $s")
        }

        fun parseFormatString(fmtStr: String): Byte {
            var s = fmtStr.toLowerCase()
            check(s.startsWith("bed")) {
                "Format name is required to start with 'bed' prefix, e.g. bed3, bed6+ or bed3+6, but was: $fmtStr"
            }
            s = s.substring(3)
            val chunks = s.split('+')
            check(chunks.size <= 2) {
                "No more than one '+' is expected, but was: $fmtStr"
            }
            return chunks[0].toByte()
        }

        fun detectDelimiter(path: Path) = detectDelimiter(path.toUri())

        fun detectDelimiter(source: URI) = source.reader().use { reader ->
            val csvFile = source.hasExt("csv")
            val maxAttempts = 1000
            var currentAttempt = 0

            var delimiter: Char? = null
            while (delimiter == null && currentAttempt < maxAttempts) {
                currentAttempt++
                val line = reader.readLine() ?: break

                if (!NON_DATA_LINE_PATTERN.matches(line)) {
                    val delimiterCandidate = when {
                        // normally *.csv files with \t are actually tab separated and
                        // comma could be from name field
                        csvFile && !line.contains('\t') -> ','
                        else -> '\t'
                    }
                    delimiter = detectDelimiterFromLine(line, delimiterCandidate)
                }
            }
            delimiter ?: (if (csvFile) ',' else '\t')
        }

        internal fun detectDelimiter(text: String, defaultDelimiter: Char): Char {
            for (line in text.lineSequence()) {
                if (NON_DATA_LINE_PATTERN.matches(line)) {
                    continue
                }

                return detectDelimiterFromLine(line, defaultDelimiter)
            }
            return defaultDelimiter
        }

        private fun detectDelimiterFromLine(line: String, delimiter: Char) = when {
            delimiter in line -> delimiter
            '\t' in line -> '\t'
            ' ' in line -> ' '
            else -> delimiter
        }

        fun auto(entry: BedEntry): BedFormat {
            val delimiter = if (entry.rest.isEmpty()) {
                '\t'
            } else {
                detectDelimiterFromLine(entry.rest, '\t')
            }
            val line = BedPrinter.toLine(entry, delimiter)
            return detectFormatFromLine(delimiter, line, null)
        }

        fun auto(path: Path) = auto(path.toUri())
        fun auto(text: String, source: String?, defaultDelimiter: Char = '\t') = auto(
                text.byteInputStream(), detectDelimiter(text, defaultDelimiter), source
        )
        fun auto(source: URI) = source.stream().use { stream ->
            auto(stream, detectDelimiter(source), source.presentablePath())
        }

        private fun auto(stream: InputStream, delimiter: Char, source: String?) =
                stream.bufferedReader().use { reader ->
                    var format = from("bed3", delimiter)
                    var headerCandidateMet = false
                    while (true) {
                        val line = reader.readLine() ?: break

                        if (NON_DATA_LINE_PATTERN.matches(line)) {
                            continue
                        }

                        // Try to parse 1st line, if failed it could be header
                        // with 'start' 'end' fields instead of integer offsets,
                        // so try again with 2nd line, if fails again => throw an
                        // error
                        format = try {
                            detectFormatFromLine(delimiter, line, source)
                        } catch (e: IllegalArgumentException) {
                            if (!headerCandidateMet) {
                                headerCandidateMet = true
                                // try again with next line
                                continue
                            } else {
                                // give up
                                throw e
                            }
                        }
                        break
                    }
                    format
                }

        private fun detectFormatFromLine(delimiter: Char, line: String, source: String?): BedFormat {
            val chunks = Splitter.on(delimiter).trimResults().omitEmptyStrings().split(line).toList()

            // https://genome.ucsc.edu/FAQ/FAQformat.html#format1:
            //
            // BED format provides a flexible way to define the data lines that are displayed
            // in an annotation track. BED lines have three required fields and nine additional
            // optional fields. The number of fields per line must be consistent throughout
            // any single set of data in an annotation track. The order of the optional fields
            // is binding: lower-numbered fields must always be populated if higher-numbered fields are used.
            //
            // Order:
            //      chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb,
            //      blockCount, blockSizes, blockStarts
            //
            // Extended BED format: bedN[+[P]]
            //   * N is between 3 and 15,
            //   * optional (+) if extra "bedPlus" fields,
            //   * optional P specifies the number of extra fields. Not required, but preferred.
            //   * Examples: -type=bed6 or -type=bed6+ or -type=bed6+3
            //   * ( see http://genome.ucsc.edu/FAQ/FAQformat.html#format1

            var ci = 0

            var blockCount = 0
            for (field in BedField.values()) {
                try {
                    val value = field.field.read(chunks[ci])

                    @Suppress("NON_EXHAUSTIVE_WHEN")
                    when (field) {
                        BedField.BLOCK_COUNT -> blockCount = value as Int
                        BedField.BLOCK_STARTS, BedField.BLOCK_SIZES -> {
                            check(when (value) {
                                null -> blockCount == 0
                                else -> blockCount <= (value as IntArray).size
                            })
                        }
                    }

                    if (ci == chunks.size) {
                        break
                    }
                    ci++
                } catch (ignored: Throwable) {
                    // Do not skip fields, just stop
                    break
                }
            }

            val fieldsNumber = ci

            require(fieldsNumber in 3..15) {
                val fileInfo = if (source != null) "Source: $source\nUnknown BED format:\n" else ""
                "$fileInfo$line\nFields number in BED file is between 3 and 15, but was $fieldsNumber"
            }

            return BedFormat(ci.toByte(), delimiter)
        }


    }
}

/**
 * Reads [BedEntry]-s from the supplied [reader].
 *
 * [BedParser] is usually used via [BedFormat.parse] method. It implements [Iterable] over [BedEntry],
 * so it's normally a receiver of [map] or [forEach] calls:
 *      bedFormat.parse(file) { parser ->
 *          parser.forEach { bedEntry -> println(bedEntry) }
 *      }
 *
 * @property reader The reader that provides lines to convert into [BedEntry]-s.
 * @property source The BED source description, used in log and exception messages.
 * @property separator The separator character for line parsing.
 *
 * @property stringency controls the parser stringency (see [Stringency] for detailed explanation).
 * It can be changed at any time:
 *      parser.stringency = Stringency.STRICT
 * It's set to [LENIENT] by default.
 *
 * @property linesFailedToParse shows how many lines were skipped due to parsing errors. Note that it's useless
 * with [stringency] == [STRICT], since the parser will throw on any line that fails to parse.
 */
class BedParser(
        internal val reader: BufferedReader,
        private val source: String,
        val separator: Char
): Iterable<BedEntry>, AutoCloseable {

    private val splitter = Splitter.on(separator).limit(4).trimResults().omitEmptyStrings()

    var stringency = LENIENT

    var linesFailedToParse: Int = 0
        private set

    override fun iterator(): UnmodifiableIterator<BedEntry> {
        return object : UnmodifiableIterator<BedEntry>() {
            private var nextEntry: BedEntry? = parseNext()

            override fun next(): BedEntry {
                val entry = nextEntry ?: throw NoSuchElementException()
                nextEntry = parseNext()
                return entry
            }

            override fun hasNext(): Boolean = nextEntry != null

        }
    }

    /**
     * Parses the next line. Depending on [stringency], either throws [BedFormatException]
     * or keeps going on parsing errors.
     * @return the parsed entry or null if there are no more lines to parse.
     */
    private fun parseNext(): BedEntry? {
        do {
            val line = readLine() ?: return null
            try {
                return parse(line)
            } catch (e: IllegalArgumentException) {
                val message = "$source: failed to parse BED line:\n$line"
                when (stringency) {
                    STRICT -> throw BedFormatException(message, e)
                    LENIENT -> {
                        linesFailedToParse++
                        LOG.debug(message, e)
                    }
                }
            }
        // this statement is normally unreachable
        } while (true)
    }

    /**
     * Returns the next line that doesn't match [NON_DATA_LINE_PATTERN], or null if end-of-file was reached.
     */
    private fun readLine(): String? {
        var line: String?
        do {
            line = reader.readLine()
        } while (line != null && NON_DATA_LINE_PATTERN.matches(line))
        return line
    }

    /**
     * Returns the parsed [BedEntry], or throws [IllegalArgumentException] if parsing failed.
     */
    private fun parse(line: String): BedEntry {
        @Suppress("UnstableApiUsage")
        val chunks = splitter.splitToList(line)
        return try {
            BedEntry(
                    chunks[0], chunks[1].toInt(), chunks[2].toInt(),
                    if (chunks.size == 3) "" else chunks[3]
            )
        } catch (e: Exception) {
            throw IllegalArgumentException("invalid BED: '$line'", e)
        }
    }

    override fun close() {
        @Suppress("UnstableApiUsage")
        Closeables.closeQuietly(reader)
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(BedParser::class.java)

        /**
         * [STRICT] parser throws a [BedFormatException] on lines it can't parse.
         * [LENIENT] parser just logs the error and continues.
         */
        enum class Stringency {
            STRICT, LENIENT
        }
    }
}

class BedPrinter(private val writer: BufferedWriter, private val format: BedFormat) : AutoCloseable {

    fun print(line: String) {
        writer.write(line)
        writer.newLine()
    }

    fun print(entry: BedEntry) {
        print(toLine(entry, format.delimiter))
    }

    fun print(entry: ExtendedBedEntry) {
        print(toLine(entry, format))
    }

    @Suppress("UnstableApiUsage")
    override fun close() = Closeables.close(writer, true)

    companion object {
        fun toLine(entry: BedEntry, delimiter: Char): String {
            val restStr = if (entry.rest.isEmpty()) "" else "$delimiter${entry.rest}"
            return "${entry.chrom}$delimiter${entry.start}$delimiter${entry.end}$restStr"
        }

        fun toLine(entry: ExtendedBedEntry, format: BedFormat): String {
            val delimiter = format.delimiter

            val prefix = "${entry.chrom}$delimiter${entry.start}$delimiter${entry.end}"

            return when (format.fieldsNumber) {
                3.toByte() -> prefix
                else -> {
                    val rest = entry.pack(format.fieldsNumber, null, delimiter).rest
                    "$prefix$delimiter$rest"
                }
            }
        }

    }
}

abstract class AbstractBedField<out T>(val name: String, val index: Byte) {
    abstract fun read(value: String): T

    override fun toString(): String = name
}

enum class BedField(val field: AbstractBedField<Any?>) {
    CHROMOSOME(object : AbstractBedField<String>("chrom", 0) {
        override fun read(value: String): String {
            // chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
            // or scaffold (e.g. scaffold10671).
            return value
        }
    }),

    START_POS(object : AbstractBedField<Int>("chromStart", 1) {
        override fun read(value: String): Int {
            // chromStart - The starting position of the feature in the chromosome
            // or scaffold. The first base in a chromosome is numbered 0.
            return value.toInt()
        }
    }),

    END_POS(object : AbstractBedField<Int>("chromEnd", 2) {
        override fun read(value: String): Int {
            // chromEnd - The ending position of the feature in the chromosome
            // or scaffold. The chromEnd base is not included in the display of
            // the feature. For example, the first 100 bases of a chromosome are
            // defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99
            return value.toInt()
        }
    }),

    NAME(object : AbstractBedField<String>("name", 3) {
        override fun read(value: String): String {
            // name - Defines the name of the BED line. This label is displayed to the left
            // of the BED line in the Genome Browser window when the track is open to full
            // display mode or directly to the left of the item in pack mode.
            if (value.length == 1) {
                val ch = value.first()
                check(ch != '+' && ch != '-') { "Name expected, but was: $value" }
            }
            return value
        }
    }),

    SCORE(object : AbstractBedField<Int>("score", 4) {
        override fun read(value: String): Int {
            // score - A score between 0 and 1000. If the track line useScore attribute is set
            // to 1 for this annotation data set, the score value will determine the level of gray
            // in which this feature is displayed (higher numbers = darker gray). This table shows
            // the Genome Browser's translation of BED score values into shades of gray ..
            return value.toInt()
        }
    }),

    STRAND(object : AbstractBedField<Char>("strand", 5) {
        override fun read(value: String): Char {
            // strand - Defines the strand - either '+' or '-'.
            check(value.length == 1)
            val ch = value.first()
            check(ch == '+' || ch == '-' || ch == '.') {
                "expected one of \"+-.\", but got: $value"
            }
            return ch
        }
    }),

    THICK_START(object : AbstractBedField<Int>("thickStart", 6) {
        override fun read(value: String): Int {
            // thickStart - The starting position at which the feature is drawn thickly
            // (for example, the start codon in gene displays). When there is no thick part,
            // thickStart and thickEnd are usually set to the chromStart position.
            return value.toInt()
        }
    }),

    THICK_END(object : AbstractBedField<Int>("thickEnd", 7) {
        override fun read(value: String): Int {
            // thickEnd - The ending position at which the feature is drawn thickly
            // (for example, the stop codon in gene displays).
            return value.toInt()
        }
    }),

    ITEM_RGB(object : AbstractBedField<Int>("itemRgb", 8) {
        override fun read(value: String): Int {
            // itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb
            // attribute is set to "On", this RBG value will determine the display color of the data
            // contained in this BED line. NOTE: It is recommended that a simple color scheme
            // (eight colors or less) be used with this attribute to avoid overwhelming the color
            // resources of the Genome Browser and your Internet browser.
            return if (value == "0" || value == ".") {
                0
            } else {
                val chunks = value.splitToInts(3)
                Color(chunks[0], chunks[1], chunks[2]).rgb
            }
        }
    }),

    BLOCK_COUNT(object : AbstractBedField<Int>("blockCount", 9) {
        override fun read(value: String): Int {
            // blockCount - The number of blocks (exons) in the BED line.
            return value.toInt()
        }
    }),

    BLOCK_SIZES(object : AbstractBedField<IntArray?>("blockSizes", 10) {
        override fun read(value: String): IntArray? {
            if (value == ".") {
                return null
            }
            // blockSizes - A comma-separated list of the block sizes. The number of items in this
            // list should correspond to blockCount.
            return value.splitToInts(-1)
        }
    }),

    BLOCK_STARTS(object : AbstractBedField<IntArray?>("blockStarts", 11) {
        override fun read(value: String): IntArray? {
            if (value == ".") {
                return null
            }
            // blockStarts - A comma-separated list of block starts. All of the blockStart positions
            // should be calculated relative to chromStart. The number of items in this list should
            // correspond to blockCount
            return value.splitToInts(-1)
        }
    }),
}

fun Location.toBedEntry() = ExtendedBedEntry(
        chromosome.name, startOffset, endOffset, strand = strand.char
)

fun ExtendedBedEntry.toLocation(genome: Genome) =
        Location(start, end, Chromosome(genome, chrom), strand.toStrand())

fun LocationsMergingList.saveWithUniqueNames(bedPath: Path) {
    val printer = BedFormat().print(bedPath)
    var i = 0
    locationIterator().forEach { (start, end, chr) ->
        printer.print(ExtendedBedEntry(chr.name, start, end, "$i", strand = '+'))
        i++
    }
    printer.close()
}

fun String.splitToInts(size: Int): IntArray {
    val s = Splitter.on(',').split(this).toList()

    // actual fields my be less that size, but not vice versa
    check(s.isNotEmpty() == (size != 0))

    val actualSize = if (size < 0) s.size else size

    val chunks = IntArray(actualSize)
    for (i in 0 until minOf(s.size, actualSize)) {
        val chunk = s[i]
        if (chunk.isNotEmpty()) {
            chunks[i] = chunk.toInt()
        }
    }
    return chunks
}

fun BedEntry.unpack(format: BedFormat, omitEmptyStrings: Boolean = false) =
        unpack(format.fieldsNumber, true, format.delimiter, omitEmptyStrings)

/**
 * Only unpack regular BED fields, omit any extra ones.
 */
fun BedEntry.unpackRegularFields(format: BedFormat, omitEmptyStrings: Boolean = false) =
        unpack(format.fieldsNumber, false, format.delimiter, omitEmptyStrings)

interface BedProvider {
    val bedSource: URI
}

class BedFormatException(message: String, cause: Throwable) : Exception(message, cause)