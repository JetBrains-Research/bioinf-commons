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
    private val TRANSCRIPTS_MAPS_CACHE
            = Maps.newConcurrentMap<Pair<Genome, GeneAliasType>, ListMultimap<String, Transcript>>()

    private val LOG = Logger.getLogger(GeneResolver::class.java)

    fun getAnyGene(genome: Genome, anyAlias: String): Gene? {
        val geneIds = matching(genome, anyAlias).asSequence()
                .map { it.ensemblGeneId }
                .distinct().toList()

        if (geneIds.size > 1) {
            LOG.warn("$anyAlias resolved to ${geneIds.size} genes: ${geneIds.joinToString(",")}")
        }

        return genome.genes.find { it.ensemblGeneId in geneIds}
    }


    /**
     * Returns a transcript with a given name.
     *
     * If the name is ambiguous then any of the matching genes might
     * be returned.
     */
    fun getAny(genome: Genome, anyAlias: String): Transcript? {
        val it = matching(genome, anyAlias).iterator()
        val match = if (it.hasNext()) it.next() else null
        if (it.hasNext()) {
            val matches = listOf(match) + it
            LOG.warn("$anyAlias resolved to ${matches.size} genes: $matches")
        }

        return match
    }

    fun getAny(genome: Genome, alias: String, aliasType: GeneAliasType): Transcript? {
        val it = matching(genome, alias, aliasType).iterator()
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
    fun get(genome: Genome, anyAlias: String): Sequence<Transcript> {
        return matching(genome, anyAlias)
    }

    fun get(genome: Genome, alias: String, aliasType: GeneAliasType): Sequence<Transcript> {
        return matching(genome, alias, aliasType)
    }

    private fun matching(genome: Genome, anyAlias: String): Sequence<Transcript> {
        return GeneAliasType.values().asSequence()
                .flatMap { aliasType -> matching(genome, anyAlias, aliasType) }
    }

    private fun matching(genome: Genome, alias: String, aliasType: GeneAliasType): Sequence<Transcript> {
        val genesMap = transcriptsMapFor(genome, aliasType)
        return genesMap[alias.toUpperCase()].asSequence()
    }

    private fun transcriptsMapFor(genome: Genome, aliasType: GeneAliasType): ListMultimap<String, Transcript> {
        return TRANSCRIPTS_MAPS_CACHE.computeIfAbsent(genome to aliasType) {
            val transcriptsMap = ImmutableListMultimap.builder<String, Transcript>()
            for (transcript in genome.transcripts) {
                val name = transcript.names[aliasType] ?: ""
                if (name.isNotEmpty()) {
                    transcriptsMap.put(name.toUpperCase(), transcript)
                }
            }

            transcriptsMap.build()
        }
    }
}
