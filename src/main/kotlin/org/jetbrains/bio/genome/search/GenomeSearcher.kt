package org.jetbrains.bio.genome.search

import com.google.common.collect.ImmutableList
import htsjdk.samtools.fastq.FastqReader
import joptsimple.BuiltinHelpFormatter
import joptsimple.OptionParser
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.search.sa.SuffixArray
import org.jetbrains.bio.genome.search.sa.SuffixArraySearcher
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.util.*
import java.util.concurrent.TimeUnit
import java.util.stream.Stream
import java.util.stream.StreamSupport

/**
 * Genome searcher, mismatches and extended alphabet patterns supported.
 * @author: Alexey Dievsky
 * @author Oleg Shpynov
 */
class GenomeSearcher(genomeQuery: GenomeQuery, mismatches: Int) {

    private val searchers: List<ChromosomeSearcher>

    init {
        genomeQuery.get().forEach {
            val path = SuffixArray.getPath(it)
            if (path.notExists || !path.isReadable) {
                SuffixArray.create(it)
            }
        }
        searchers = genomeQuery.get().map { ChromosomeSearcher(SuffixArraySearcher(it), mismatches, it) }
    }

    fun find(text: String): Stream<Location> {
        return searchers.parallelStream().flatMap { it.find(text) }
    }

    fun findUnique(text: String): Location? {
        val iterator = searchers.mapNotNull { it.findUnique(text) }.iterator()
        if (!iterator.hasNext()) {
            return null
        }
        val location = iterator.next()
        return if (iterator.hasNext()) null else location
    }

    /**
     * Returns one of [UNIQUE], [NOT_UNIQUE], [NOT_FOUND]
     */
    fun test(text: String): Int {
        var foundUnique = false
        for (searcher in searchers) {
            when (searcher.test(text)) {
                ChromosomeSearcher.NOT_UNIQUE -> return ChromosomeSearcher.NOT_UNIQUE
                ChromosomeSearcher.UNIQUE -> if (!foundUnique)
                    foundUnique = true
                else
                    return ChromosomeSearcher.NOT_UNIQUE
            }
        }
        return if (foundUnique)
            ChromosomeSearcher.UNIQUE
        else
            ChromosomeSearcher.NOT_FOUND
    }

    companion object {
        @JvmStatic
        fun getPath(
            input: Path, mismatches: Int,
            reference: String, unique: Boolean
        ): Path = "${input}_${reference}_$mismatches${if (unique) "_unique" else ""}.bed".toPath()


        private val optionsParser = object : OptionParser() {
            init {
                acceptsAll(ImmutableList.of("i", "input"), "fastq or fastq.gz file")
                    .withRequiredArg()
                    .required()
                    .withValuesConvertedBy(PathConverter.exists())

                acceptsAll(
                    ImmutableList.of("r", "reference"),
                    "reference genome, possible with chromosomes, i.e. mm9[1,2,3]"
                )
                    .withRequiredArg()
                    .required()

                acceptsAll(ImmutableList.of("m", "mismatches"), "number of mismatches allowed")
                    .withRequiredArg()
                    .ofType(Int::class.java)
                    .defaultsTo(0)

                acceptsAll(ImmutableList.of("u", "unique"), "unique")

                formatHelpWith(BuiltinHelpFormatter(300, 2))
            }
        }


        @JvmStatic
        fun main(args: Array<String>) {
            optionsParser.parse(args) { options ->
                val input = options.valueOf("input") as Path
                val reference = options.valueOf("reference") as String
                val mismatches = options.valueOf("mismatches") as Int
                val unique = "unique" in options
                val searcher = GenomeSearcher(Genome[reference].toQuery(), mismatches)
                // Fastq record is 2 lines text, header and sequence
                val records: Long = input.bufferedReader().use { it.lines().count() / 2 }
                val builder = Progress {
                    title = "Align $input on $reference${if (unique) " unique" else ""} with $mismatches mismatches"
                    period = 5 to TimeUnit.MINUTES
                }

                val progress = if (records > 0) builder.bounded(records) else builder.unbounded()
                val reader = FastqReader(input.toFile())
                CSVFormat.TDF.print(getPath(input, mismatches, reference, unique).bufferedWriter()).use { printer ->
                    val sequences = StreamSupport.stream(
                        Spliterators.spliteratorUnknownSize(reader, Spliterator.SIZED), false
                    )
                        .map { it.readString }
                        .map { it.replace("n".toRegex(), "") }
                        .filter { it.isNotEmpty() }
                        .parallel()
                    val alignment = if (unique) {
                        sequences.map { seq ->
                            progress.report()
                            val location = searcher.findUnique(seq)
                            if (location != null) seq to location else null
                        }.filter { it != null }
                    } else {
                        sequences.flatMap { seq ->
                            progress.report()
                            searcher.find(seq).map { seq to it }
                        }
                    }
                    alignment.forEach {
                        synchronized(printer) {
                            val (_, l) = it!!
                            printer.printRecord(
                                l.chromosome.name,
                                l.toRange().startOffset.toString(),
                                l.toRange().endOffset.toString(),
                                "",
                                l.strand.toString()
                            )
                        }
                    }
                    progress.done()
                }
            }
        }
    }
}
