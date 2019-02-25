package org.jetbrains.bio.genome

import com.google.common.collect.ImmutableListMultimap
import com.google.common.collect.ListMultimap
import com.google.common.collect.Maps
import org.apache.log4j.Logger

/**
 * You've got names? We've got genes!
 *
 * @see GeneAliasType
 * @author Sergei Lebedev
 */
object GeneResolver {
    /** A mapping of UPPERCASE "gene names" to genes for a specific organism.  */
    private val GENES_MAPS_CACHE
            = Maps.newConcurrentMap<Pair<String, GeneAliasType>, ListMultimap<String, Transcript>>()

    private val LOG = Logger.getLogger(GeneResolver::class.java)

    fun getAnyGene(build: String, anyAlias: String): Gene? {
        val geneIds = matching(build, anyAlias).asSequence()
                .map { it.ensemblGeneId }
                .distinct().toList()

        if (geneIds.size > 1) {
            LOG.warn("$anyAlias resolved to ${geneIds.size} genes: ${geneIds.joinToString(",")}")
        }

        return Genome[build].genes.find { it.ensemblGeneId in geneIds}
    }


    /**
     * Returns a transcript with a given name.
     *
     * If the name is ambiguous then any of the matching genes might
     * be returned.
     */
    fun getAny(build: String, anyAlias: String): Transcript? {
        val it = matching(build, anyAlias).iterator()
        val match = if (it.hasNext()) it.next() else null
        if (it.hasNext()) {
            val matches = listOf(match) + it
            LOG.warn("$anyAlias resolved to ${matches.size} genes: $matches")
        }

        return match
    }

    fun getAny(build: String, alias: String, aliasType: GeneAliasType): Transcript? {
        val it = matching(build, alias, aliasType).iterator()
        val match = if (it.hasNext()) it.next() else null
        if (it.hasNext()) {
            val matches = listOf(match) + it
            LOG.warn("$alias ($aliasType) resolved to ${matches.size} genes: $matches")
        }

        return match
    }

    /**
     * Returns all genes with a given name.
     */
    fun get(build: String, anyAlias: String): Sequence<Transcript> {
        return matching(build, anyAlias)
    }

    fun get(build: String, alias: String, aliasType: GeneAliasType): Sequence<Transcript> {
        return matching(build, alias, aliasType)
    }

    private fun matching(build: String, anyAlias: String): Sequence<Transcript> {
        return GeneAliasType.values().asSequence()
                .flatMap { aliasType -> matching(build, anyAlias, aliasType) }
    }

    private fun matching(build: String, alias: String, aliasType: GeneAliasType): Sequence<Transcript> {
        val genesMap = genesMapFor(build, aliasType)
        return genesMap[alias.toUpperCase()].asSequence()
    }

    private fun genesMapFor(build: String, aliasType: GeneAliasType): ListMultimap<String, Transcript> {
        return GENES_MAPS_CACHE.computeIfAbsent(build to aliasType) {
            val genesMap = ImmutableListMultimap.builder<String, Transcript>()
            for (gene in Genome[build].transcripts) {
                val name = gene.names[aliasType] ?: ""
                if (name.isNotEmpty()) {
                    genesMap.put(name.toUpperCase(), gene)
                }
            }

            genesMap.build()
        }
    }
}
