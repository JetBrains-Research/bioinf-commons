package org.jetbrains.bio.genome

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.text.SimpleDateFormat
import java.util.*

/**
 * @author Roman.Chernyatchik
 */
object AnnotationsConfig {
    private val LOG = Logger.getLogger(AnnotationsConfig::class.java)
    const val VERSION: Int = 3

    private val mapping: Map<String, GenomeAnnotationsConfig> = load()

    val builds = mapping.keys
    operator fun get(build: String): GenomeAnnotationsConfig {
        if (build == Genome.TEST_ORGANISM_BUILD) {
            return GenomeAnnotationsConfig("Test Organism", null, "n/a", emptyMap(),
                    false, "n/a", "n/a", "n/a", "n/a",
                    "n/a", null, null, null)
        }
        check(build in builds) {
            "Unexpected build name $build, annotations are available only for ${builds.joinToString()}"
        }
        return mapping.getValue(build)
    }


    internal fun load(path: Path): Pair<Int, Map<String, GenomeAnnotationsConfig>> {
        val yaml = YamlMapper.load(path)
        val buildsConfig = yaml.genomes.keys.map { build ->
            build to yaml[build]
        }.toMap()
        return yaml.version to buildsConfig
    }

    private fun load(): Map<String, GenomeAnnotationsConfig> {
        val path = Configuration.genomesPath / "annotations.yaml"
        if (path.exists) {
            val (version, mapping) = load(path)
            if (version == VERSION) {
                return mapping
            }

            val dateFormat = SimpleDateFormat("yyyy-MM-dd")
            val hash = mapping.hashCode().toString().sha
            val suffix = "${dateFormat.format(Date())}.$hash"

            val bkPath = path.withExtension("$suffix.yaml")
            path.move(bkPath, StandardCopyOption.ATOMIC_MOVE)
            LOG.info("Outdated annotations file $path (version: $version) was backed up to $bkPath " +
                    "and replaced with recent version: $VERSION).")
        }
        path.checkOrRecalculate { output ->
            val text = AnnotationsConfig::class.java
                    .getResourceAsStream("/annotations.yaml")
                    .bufferedReader(Charsets.UTF_8)
                    .readText()
            output.path.write(text)
        }
        val (newVersion, newMapping) = AnnotationsConfig.load(path)
        require(newVersion == VERSION) {
            "Bundled annotations are expected to have latest format version $VERSION but was ${newVersion}."
        }
        return newMapping
    }

    internal class YamlMapper {
        @JvmField
        var version: Int = 0

        /**
         * A temporary object for loading weakly-typed YAML data.
         */
        @JvmField
        var genomes = LinkedHashMap<String, LinkedHashMap<String, Any>>()

        operator fun get(build: String): GenomeAnnotationsConfig {
            val genomeAttrs = genomes[build]
            if (genomeAttrs == null) {
                checkNotNull(genomeAttrs) { "Annotations not defined for $build" }
            }

            val gtf2UCSCChrNamesMapping = if ("gtf_chrs" in genomeAttrs) {
                (genomeAttrs["gtf_chrs"] as List<Map<String, String>>).map {
                    it.entries.first().let { (k, v) -> k to v }
                }.toMap()
            } else {
                emptyMap()
            }
            val biomart = genomeAttrs["biomart"]
            val mart = if (biomart != null) {
                val biomartMap = biomart as Map<String, String>
                val url = biomartMap["url"] as String
                val dataset = biomartMap["dataset"] as String
                Mart(dataset, url)
            } else {
                null
            }
            return GenomeAnnotationsConfig(
                    genomeAttrs["species"] as String,
                    genomeAttrs["alias"] as String?,
                    genomeAttrs["gtf"] as String,
                    gtf2UCSCChrNamesMapping,
                    genomeAttrs["ucsc_annotations_legacy"]?.toString()?.toBoolean() ?: false,
                    genomeAttrs["sequence"] as String,
                    genomeAttrs["chromsizes"] as String,
                    genomeAttrs["repeats"] as String,
                    genomeAttrs["cytobands"] as String?,
                    genomeAttrs["gaps"] as String,
                    genomeAttrs["centromeres"] as String?,
                    genomeAttrs["cgis"] as String?,
                    mart
            )
        }

        companion object {
            private fun createMapper(): ObjectMapper {
                val yamlFactory = YAMLFactory()
                yamlFactory.configure(YAMLGenerator.Feature.MINIMIZE_QUOTES, true)
                return ObjectMapper(yamlFactory)
            }

            /** Loads configuration from a YAML file with id as file name. */
            fun load(path: Path) = path.bufferedReader().use { reader ->
                createMapper().readValue(reader, YamlMapper::class.java)
            }!!

        }
    }
}

/**
 * @param species Species name
 * @param alias Genome build alternative name
 * @param gtfUrl GTF genes annotations url
 * @param gtfChrsMapping GTF to UCSC chr names custom mapping
 * @param ucscAnnLegacyFormat Legacy format provides separate file for each chromosome
 * @param sequenceUrl 2bit file url
 * @param chromsizesUrl Chromosomes size file url
 * @param repeatsUrl Repeats file url
 * @param cytobandsUrl CytoBands file url
 * @param gapsUrl Gaps file url
 * @param centromeresUrl Centromeres file url. The latest human annotations defines
 *  centromeres info in this file instead of gaps
 * @param cpgIslandsUrl CpG islands file url.
 */
data class GenomeAnnotationsConfig(
        val species: String,
        val alias: String?,
        val gtfUrl: String,
        val gtfChrsMapping: Map<String, String>,
        val ucscAnnLegacyFormat: Boolean,
        val sequenceUrl: String,
        val chromsizesUrl: String,
        val repeatsUrl: String,
        val cytobandsUrl: String?,
        val gapsUrl: String,
        val centromeresUrl: String?,
        val cpgIslandsUrl: String?,
        val mart: Mart?
)
