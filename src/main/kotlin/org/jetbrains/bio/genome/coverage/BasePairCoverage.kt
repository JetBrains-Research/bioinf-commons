package org.jetbrains.bio.genome.coverage

import com.google.common.base.MoreObjects
import gnu.trove.list.TIntList
import gnu.trove.list.array.TIntArrayList
import gnu.trove.set.hash.TIntHashSet
import kotlinx.support.jdk7.use
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode
import org.jetbrains.bio.dataframe.DataFrameMappers
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.npy.NpzFile
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.IOException
import java.nio.file.Path

/**
 * The container maintains a sorted strand-independent list of single nucleotide offsets for each chromosome.
 * In future could be updated for stranded data, but for current usecase it is not required.
 *
 * Class could be used to represent all CpG covered by 'merged' CpG methylome. Main difference with SingleEndCoverage
 * that it doesn't use `tag` shift approach, doesn't calculate fragment size and is filled by single offsets,
 * not by SE read intervals.
 *
 * Immutable. Saves data in [NpzFile] format.
 *
 * @author Roman Cherniatchik
 */
class BasePairCoverage private constructor(
    override val genomeQuery: GenomeQuery,
    val data: GenomeMap<TIntList>
) : Coverage {
    /**
     * Calculates how many offsets observed in current location. Location strand is ignored.
     *
     * @param location the Location, strand will be ignored
     * @return an integer representing the coverage value for the given location
     */
    override fun getCoverage(location: Location): Int {
        val data = data[location.chromosome]

        /* we don't really care if the offsets are outside of chromosome range,
           since this won't lead to incorrect results */

        val startOffset = location.startOffset
        val endOffset = location.endOffset
        val index = data.binarySearchLeft(startOffset)
        var size = 0
        while (index + size < data.size() && data[index + size] < endOffset) {
            size++
        }
        return size
    }

    override val depth = genomeQuery.get().map { chr -> data[chr].size().toLong() }.sum()

    override fun toString() = MoreObjects.toStringHelper(this)
        .addValue(genomeQuery).toString()

    fun filter(
        regions: LocationsMergingList,
        includeRegions: Boolean,
        progress: Boolean = false,
        ignoreRegionsOnMinusStrand: Boolean = true
    ): BasePairCoverage {
        val builder = Builder(genomeQuery, false)

        val rangeLists = regions.rangeLists

        val progressInd = if (progress) {
            Progress { title = "Filtering coverage" }
                .bounded(depth)
        } else {
            null
        }

        val strands = if (ignoreRegionsOnMinusStrand) arrayOf(Strand.PLUS) else Strand.values()
        genomeQuery.get().forEach { chr ->
            for (offset in data[chr]) {
                var included = false
                for (strand in strands) {
                    if (!rangeLists.contains(chr, strand)) {
                        continue
                    }

                    val rl = rangeLists[chr, strand]
                    if (rl.includesRange(offset)) {
                        included = true
                        break
                    }
                }

                if (includeRegions == included) {
                    // if both are false or both are true
                    builder.process(chr, offset)
                }

                progressInd?.report()
            }
        }
        progressInd?.done()

        return builder.build(false) // preserve same uniqueness state as in this
    }


    @Throws(IOException::class)
    fun saveToNpz(outputPath: Path) {
        NpzFile.write(outputPath).use { writer ->
            writer.write(Coverage.VERSION_FIELD, intArrayOf(Coverage.VERSION))
            writer.write(Coverage.COV_TYPE_FIELD, intArrayOf(Coverage.CoverageType.BASEPAIR.ordinal))
            writer.write(BASEPAIR_VERSION_FIELD, intArrayOf(BASEPAIR_VERSION))

            for (chromosome in genomeQuery.get()) {
                val key = chromosome.name
                writer.write(key, data[chromosome].toArray())
            }
        }
    }

    @Throws(IOException::class)
    fun saveToTSV(outputPath: Path, offsetIsOneBased: Boolean = true) {
        DataFrameMappers.TSV.format.print(outputPath.bufferedWriter()).use { csvPrinter ->
            for (chromosome in genomeQuery.get()) {
                val chrName = chromosome.name
                for (offset in data[chromosome].toArray()) {
                    csvPrinter.printRecord(chrName, offset + (if (offsetIsOneBased) 1 else 0))
                }
            }
        }
    }

    class Builder(val genomeQuery: GenomeQuery, val offsetIsOneBased: Boolean) {
        val data: GenomeMap<TIntList> = genomeMap(genomeQuery) { TIntArrayList() }

        private var basePairsCount = 0L


        /**
         * Add an offset to the coverage being built.
         */
        fun process(chr: Chromosome, offset: Int): Builder {
            if (offsetIsOneBased) {
                require(offset >= 1) { "One-based offset should be >= 1, but was: $offset" }
                require(offset <= chr.length) { "One-based offset should be <= ${chr.length} (${chr.name} size), but was: $offset" }
            } else {
                require(offset >= 0) { "Zero-based offset should be >= 1, but was: $offset" }
                require(offset < chr.length) { "One-based offset should be < ${chr.length} (${chr.name} size), but was: $offset" }
            }

            val shift = if (offsetIsOneBased) 1 else 0
            data[chr].add(offset - shift)
            basePairsCount++
            return this
        }


        /**
         * Generate a [BasePairCoverage] object.
         * [unique] controls whether duplicate tags should be preserved ([unique] == false)
         * or squished into one tag ([unique] == true).
         * Only offsets at the exact same offset on the exact same strand
         * are considered duplicate.
         */
        fun build(unique: Boolean): BasePairCoverage {
            for (chr in data.genomeQuery.get()) {
                if (unique) {
                    // XXX we can do linear time de-duplication on
                    //     the sorted sequence.
                    val tags = data[chr]
                    data[chr] = TIntArrayList(TIntHashSet(tags))
                }
                data[chr].sort()
            }

            val cov = BasePairCoverage(genomeQuery, data = data)
            val coverageDepth = cov.depth
            if (unique) {
                assert(coverageDepth <= basePairsCount) { "Expected $basePairsCount >= $coverageDepth" }
            } else {
                assert(coverageDepth == basePairsCount) { "Expected: $basePairsCount but was $coverageDepth" }
            }
            return cov
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(BasePairCoverage::class.java)
        const val BASEPAIR_VERSION = 3
        const val BASEPAIR_VERSION_FIELD = "basepair_version"

        fun builder(genomeQuery: GenomeQuery, offsetIsOneBased: Boolean) = Builder(genomeQuery, offsetIsOneBased)

        /**
         * Loads from TAB separated file with at least 2 columns:
         *  - 1st column: chromosome name
         *  - 2nd column: genome offset
         *
         *  @param path: File path
         *  @param offsetIsOneBased: If true offsets are considered 1-based, otherwise 0-based
         *  @param header Is first line of file header or not
         *  @param progress Whether to show progress or not
         */
        fun loadFromTSV(
            gq: GenomeQuery,
            path: Path, offsetIsOneBased: Boolean,
            header: Boolean = false,
            progress: Boolean = false,
            failOnMissingChromosomes: Boolean = true
        ): BasePairCoverage {

            val format = CSVFormat.TDF.withQuoteMode(QuoteMode.MINIMAL).withCommentMarker('#')
                .withRecordSeparator("\n")!!

            LOG.info("Loading coverage $path ${path.size}")
            val linesNumber = if (progress) {
                path.bufferedReader().use {
                    it.lines().mapToInt { line -> if (line[0] != format.commentMarker) 1 else 0 }.sum()
                }
            } else {
                0
            }

            val rowsNumber = if (linesNumber == 0) 0 else linesNumber - (if (header) 1 else 0)
            if (progress) {
                LOG.info("#Rows in [${path.fileName}]: ${rowsNumber.formatLongNumber()}")
            }

            val progressInd = if (progress) {
                Progress { title = "Reading coverage from $path" }
                    .bounded(rowsNumber.toLong())
            } else {
                null
            }

            val chromosomeNamesMap = gq.genome.chromosomeNamesMap

            val builder = Builder(gq, offsetIsOneBased)
            path.bufferedReader().use {
                for (row in format.withSkipHeaderRecord(header).parse(it)) {
                    val chrName = row[0]
                    val offset = row[1].toInt()

                    val chr = chromosomeNamesMap[chrName]
                    if (chr == null) {
                        if (failOnMissingChromosomes) {
                            throw IllegalArgumentException("Unknown chromosome '$chrName' for genome: '${gq.genome.presentableName()}'")
                        }
                        continue
                    }
                    builder.process(chr, offset)
                    progressInd?.report()
                }

                progressInd?.done()
            }
            return builder.build(true)
        }

        internal fun load(
            npzReader: NpzFile.Reader,
            path: Path,
            genomeQuery: GenomeQuery,
            failOnMissingChromosomes: Boolean
        ): BasePairCoverage {
            val covTypeByte = npzReader[Coverage.COV_TYPE_FIELD].asIntArray().single()
            check(covTypeByte == Coverage.CoverageType.BASEPAIR.ordinal) {
                "$path attempting to read other coverage (type=$covTypeByte) cache file"
            }
            try {
                val version = npzReader[BASEPAIR_VERSION_FIELD].asIntArray().single()
                check(version == BASEPAIR_VERSION) {
                    "$path basepair coverage version is $version instead of ${BASEPAIR_VERSION}"
                }
            } catch (e: IllegalStateException) {
                throw IllegalStateException("$path basepair coverage version is missing", e)
            }

            val data: GenomeMap<TIntList> = genomeMap(genomeQuery) { TIntArrayList() }
            for (chromosome in genomeQuery.get()) {
                try {
                    val npyArray = npzReader[chromosome.name]
                    data[chromosome] = TIntArrayList.wrap(npyArray.asIntArray())
                } catch (e: NullPointerException) { // JDK11 doesn't throw ISE in case of missing zip entry
                    val msg = "File $path doesn't contain data for ${chromosome.name}."
                    LOG.trace(msg)
                    if (failOnMissingChromosomes) {
                        throw java.lang.IllegalStateException(msg, e)
                    }
                } catch (e: IllegalStateException) {
                    val msg = "File $path doesn't contain data for ${chromosome.name}."
                    LOG.trace(msg)
                    if (failOnMissingChromosomes) {
                        throw java.lang.IllegalStateException(msg, e)
                    }
                }
            }
            return BasePairCoverage(genomeQuery, data = data)
        }

    }
}