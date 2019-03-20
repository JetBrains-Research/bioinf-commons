package org.jetbrains.bio.genome


import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.collect.Iterators
import com.google.common.collect.PeekingIterator
import org.apache.commons.csv.CSVRecord
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.util.checkOrRecalculate
import java.nio.file.Path
import java.util.*


object GeneDescription {
    private val CACHE: Cache<Genome, Map<String, String>> = CacheBuilder.newBuilder()
            .softValues()
            .initialCapacity(1)
            .build<Genome, Map<String, String>>()


    fun getMap(genome: Genome): Map<String, String> {
        return CACHE.get(genome) {
            val descriptionPath = genome.genesDescriptionsPath
            descriptionPath.checkOrRecalculate("Genes") { output ->
                val pairs = downloadAnnotation(genome).toList()
                DataFrame()
                        .with("name", pairs.map { it.first }.toTypedArray())
                        .with("description", pairs.map { it.second }.toTypedArray())
                        .save(output.path)
            }

            load(descriptionPath)
        }
    }

    private fun load(descriptionPath: Path): Map<String, String> {
        val dataFrame = DataFrame.load(descriptionPath)
        val names = dataFrame.sliceAsObj<String>("name")
        val descriptions = dataFrame.sliceAsObj<String>("description")
        return names.zip(descriptions).toMap()
    }


    fun downloadAnnotation(genome: Genome): Map<String, String> {

        val genesDescriptionMap = TreeMap<String, String>()
        val attributes = listOf("ensembl_gene_id", "description")
        genome.annotationsConfig?.mart?.query(attributes) { pipe ->
            val block = ArrayList<CSVRecord>(1)
            val it = Iterators.peekingIterator(pipe.iterator())
            while (it.hasNext()) {
                // Even though we query Biomart with `unique = true`, it
                // might sometimes return multiple entries for the same
                // transcript ID. And these entries often contain
                // complementary data. To work this around we combine
                // duplicates into blocks.
                it.nextBlock(block) { it["ensembl_gene_id"] }

                val ensemblGeneId = block["ensembl_gene_id"]

                val description = block["description"]
                if (description == "") {
                    continue
                }

                genesDescriptionMap[ensemblGeneId] = description
            }
        }
        return genesDescriptionMap
    }

    fun getDescription(transcript: Transcript): String? {
        val genome: Genome = transcript.chromosome.genome
        return getMap(genome)[transcript.ensemblGeneId]
    }

}


private inline fun <T> PeekingIterator<T>.nextBlock(
        into: ArrayList<T>, transform: (T) -> Any) {
    into.clear()
    do {
        into.add(next())
    } while (hasNext() && transform(peek()) == transform(into.first()))
}

private operator fun List<CSVRecord>.get(field: String): String {
    for (row in this) {
        val value = row[field]
        if (value.isNotEmpty()) {
            return value
        }
    }

    return ""
}
