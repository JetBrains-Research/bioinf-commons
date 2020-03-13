package org.jetbrains.bio.genome.query

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GeneClass
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Transcript

class ChromosomeTranscriptsQuery(private val genesClass: GeneClass)
    : Query<Chromosome, Collection<Transcript>> {

    override fun apply(input: Chromosome): Collection<Transcript> {
        return if (genesClass === GeneClass.ALL) {
            input.transcripts
        } else {
            input.transcripts.filter { it in genesClass }
        }
    }

    override val id: String get() = genesClass.id

    override val description: String get() {
        return "Chromosome genes ${genesClass.description}"
    }
}

class TranscriptsQuery @JvmOverloads constructor(private val genesClass: GeneClass = GeneClass.ALL)
    : Query<GenomeQuery, List<Transcript>> {

    override fun apply(genomeQuery: GenomeQuery): List<Transcript> {
        return genomeQuery.get().asSequence().flatMap { it.transcripts.asSequence() }
                .filter { it in genesClass }
                .toList()
    }

    override val id: String get() = genesClass.id

    override val description: String get() = genesClass.description
}
