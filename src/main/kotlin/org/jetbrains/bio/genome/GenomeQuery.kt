package org.jetbrains.bio.genome

import com.google.common.annotations.VisibleForTesting
import org.apache.log4j.Logger
import java.util.*

/**
 * Genome query for all chromosomes.
 */
class GenomeQuery (val genome: Genome, vararg names: String) {
    /** A subset of chromosomes to be considered or null for all chromosomes. */
    val restriction: Set<String>? = if (names.isNotEmpty()) names.toSet() else null

    init {
        restriction?.forEach {
            check(it in genome.chromosomeNamesMap) {
                "Unknown chromosome name for $genome: $it"
            }
        }
    }

    val build: String
        get() = genome.build

    private val chromosomes: List<Chromosome> by lazy {
        genome.chromosomes.filter { accepts(it) }
    }

    fun get(): List<Chromosome> = chromosomes
    fun accepts(chromosome: Chromosome) = restriction == null || chromosome.name in restriction

    val id: String
        get() = build + if (restriction != null) "[${restriction.sorted().joinToString(",")}]" else ""

    val description: String
        get() {
            return "$build [${restriction?.sorted()?.joinToString(", ") ?: "all chromosomes"}]"
        }

    operator fun get(name: String): Chromosome? {
        val result = genome.chromosomeNamesMap[name]
        if (result == null) {
            LOG.debug("Chromosome $name not found in genome ${genome.presentableName()}")
            return null
        }
        /** [Chromosome] always contain canonical name, i.e. given in chrom.sizes.path */
        val canonicalName = result.name
        if (restriction != null && canonicalName !in restriction) {
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

        /**
         * Parses [String] as genome with possible custom chromosomes set.
         * Is designed to support serialized [GenomeQuery.id]
         *
         * @param str Genome defining string,  e.g. "hg19" or "hg19[chr1,chr2]"
         */
        fun parseGenomeQueryId(str: String): Pair<String, Array<String>> {
            if ("[" !in str) {
                return str to emptyArray()
            }
            val build = str.substringBefore('[')
            val names = str.substringAfter('[').replace("]", "").split(',').filter { it.isNotBlank() }.toTypedArray()
            check(names.isNotEmpty()) { "Empty restriction is not allowed within []" }
            return build to names
        }
    }
}

fun Genome.toQuery() = GenomeQuery(this)