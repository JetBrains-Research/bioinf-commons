package org.jetbrains.bio.genome

import com.google.common.annotations.VisibleForTesting
import org.apache.log4j.Logger
import java.nio.file.Path
import java.util.*

/**
 * Genome query for all chromosomes (somatic & autosomal).
 * Important: unlocalized and unmapped fragments and the mitochondrial chromosome are ignored.
 *
 * @author Oleg Shpynov
 */
class GenomeQuery(val genome: Genome,
                  /** A subset of chromosomes to be considered or empty list for all chromosomes. */
                  val restriction: Set<String> = emptySet()) {

    constructor(build: String, vararg names: String) : this(Genome[build], names.toSet())

    constructor(chromSizesPath: Path) :
            this(Genome[buildByChromSizesPath(chromSizesPath), chromSizesPath])

    val build: String
        get() = genome.build

    private val filteredChromosomes: List<Chromosome> by lazy {
        genome.chromosomes.filter {
            if (restriction.isEmpty()) {
                chrDefaultChoice(it.name)
            } else {
                it.name in restriction
            }
        }
    }

    fun get(): List<Chromosome> = filteredChromosomes

    fun only(chromosomes: List<String>): GenomeQuery {
        chromosomes.forEach {
            check(it in this) {
                "Unknown chromosome name: $it for $description"
            }
        }
        if (chromosomes.sorted() == get().map { it.name }.sorted()) {
            return this
        }
        return GenomeQuery(genome, chromosomes.toSet())
    }

    val id: String
        get() = build + if (restriction.isNotEmpty()) "[${restriction.sorted().joinToString(",")}]" else ""

    val description: String
        get() {
            val contents = if (restriction.isNotEmpty()) {
                restriction.joinToString(", ")
            } else {
                "all chromosomes"
            }
            return "$build [$contents]"
        }

    operator fun get(name: String): Chromosome? {
        val result = genome.chromosomeNamesMap[name]
        if (result == null) {
            LOG.debug("Chromosome $name not found in genome ${genome.presentableName()}")
            return null
        }
        /** [Chromosome] always contain canonical name, i.e. given in chrom.sizes.path */
        val canonicalName = result.name
        if (!chrDefaultChoice(canonicalName)) {
            LOG.debug("Chromosome $name in ignored")
            return null
        }
        if (restriction.isNotEmpty() && canonicalName !in restriction) {
            LOG.debug("Chromosome $name in not in restricted genome $description")
            return null
        }
        return result
    }

    operator fun contains(name: String): Boolean = get(name) != null

    override fun toString() = id

    override fun equals(other: Any?) = when {
        other === this -> true
        other == null || other !is GenomeQuery -> false
        else -> build == other.build &&
                restriction == other.restriction
    }

    override fun hashCode() = Objects.hash(build, restriction)

    companion object {
        @VisibleForTesting
        internal val LOG = Logger.getLogger(GenomeQuery::class.java)

        private fun buildByChromSizesPath(chromSizesPath: Path): String {
            val fileName = chromSizesPath.fileName.toString()
            if (!fileName.endsWith(".chrom.sizes")) {
                val build = fileName.substringBefore(".")
                LOG.warn("Unexpected chrom sizes file name: $fileName, expected <build>.chrom.sizes. " +
                        "Detected build: $build")
                return build
            }
            val build = fileName.substringBeforeLast(".chrom.sizes")
            LOG.debug("Chrom sizes name: $fileName. Detected build: $build")
            return build
        }

        private val MAPPED_CHRS_PATTERN = "chr[0-9a-tv-zA-TV-Z]+[0-9a-zA-Z]*".toRegex()

        /**
         * By default let's ignore chrMT, unmapped contigs and alternative contigs
         */
        fun chrDefaultChoice(chrName: String) = chrName != "chrM" && chrName.matches(MAPPED_CHRS_PATTERN)
    }
}

fun Genome.toQuery() = GenomeQuery(this)

/**
 * Restores genome query from [String]
 */
fun String.toGenomeQuery(): GenomeQuery {
    if ("[" !in this) {
        return GenomeQuery(this)
    }
    val build = substringBefore('[')
    val names = substringAfter('[').replace("]", "").split(',').toTypedArray()
    return GenomeQuery(build, *names)
}