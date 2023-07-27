package org.jetbrains.bio.genome

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator
import org.jetbrains.bio.genome.AnnotationsConfigLoader.saveYaml
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.text.SimpleDateFormat
import java.util.*

/**
 * Generic genome information. Used to define [Genome].
 *
 * @param species Species name
 * @param ucscAlias UCSC name for the build
 * @param names All known build names (release name and other aliases, including the UCSC alias) as ordered set
 * @param description About Genome build
 * @param gtfUrl GTF genes annotations url
 * @param chrAltName2CanonicalMapping Mapping of alternative chr names to canonical names from *.chrom.sizes.
 *  E.g. MT -> chrM for GTF parsing with UCSC genomes
 * @param ucscAnnLegacyFormat Legacy format provides separate file for each chromosome
 * @param sequenceUrl 2bit file url
 * @param chromsizesUrl Chromosomes size file url
 * @param repeatsUrl Repeats file url
 * @param cytobandsUrl CytoBands file url
 * @param gapsUrl Gaps file url
 * @param centromeresUrl Centromeres file url. The latest human annotations defines
 *  centromeres info in this file instead of gaps
 * @param cpgIslandsUrl CpG islands file url.
 *
 * @author Roman.Chernyatchik
 */
data class GenomeAnnotationsConfig(
    val species: String,
    val ucscAlias: String?,
    val names: List<String>,
    val description: String,
    val gtfUrl: String,
    val chrAltName2CanonicalMapping: Map<String, String>,
    val ucscAnnLegacyFormat: Boolean,
    val sequenceUrl: String,
    val chromsizesUrl: String,
    val repeatsUrl: String?,
    val cytobandsUrl: String?,
    val gapsUrl: String?,
    val centromeresUrl: String?,
    val cpgIslandsUrl: String?,
    val mart: Biomart?
)

object AnnotationsConfigLoader {
    const val VERSION: Int = 4
    private val LOG = LoggerFactory.getLogger(AnnotationsConfigLoader::class.java)

    private var pathAndConfig: Pair<Path, Map<String, GenomeAnnotationsConfig>>? = null

    /**
     * Parse YAML file and load configuration. If YAML file doesn't exist if would be created with
     * default bundled content using [bioinf-commons/src/main/resources/annotations.yaml]
     */
    fun init(yamlPath: Path, checkAlreadyInitialized: Boolean = true) {
        if (!initialized || !checkAlreadyInitialized) {
            pathAndConfig = yamlPath to loadOrCreateAnnotationsConfig(yamlPath)
        } else {
            val path = pathAndConfig!!.first
            error("Already initialized with $path, couldn't be overridden with $yamlPath.")
        }
    }

    val initialized: Boolean
        get() = pathAndConfig != null

    /**
     * Ensure annotations config loaded, use [AnnotationsConfigLoader.init]
     */
    val yamlPath: Path
        get() {
            require(initialized) { "Annotations config not initialized, call [AnnotationsConfig.init] first" }
            return pathAndConfig!!.first
        }

    /**
     * Ensure annotations config loaded, use [AnnotationsConfigLoader.init]
     */
    val builds: Set<String>
        get() {
            require(initialized) { "Annotations config not initialized, call [AnnotationsConfig.init] first" }
            return pathAndConfig!!.second.keys
        }


    operator fun get(build: String): GenomeAnnotationsConfig {
        require(initialized) { "Annotations config not initialized, call [AnnotationsConfig.init] first" }
        check(build in builds) {
            "Unexpected build name $build, annotations are available only for ${builds.joinToString()}"
        }
        return pathAndConfig!!.second[build]!!
    }

    fun parseYaml(yamlPath: Path, currentVersion: Int):
            Triple<Int, Map<String, GenomeAnnotationsConfig>?, AnnotationConfigYAMLSerializable> {
        val yaml = AnnotationConfigYAMLSerializable.load(yamlPath)
        val yamlVersion = yaml.version
        if (yamlVersion != currentVersion) {
            return Triple(yamlVersion, null, yaml)
        }

        return Triple(yamlVersion, yaml.genomes.mapValues {
            val deserializedMap = yaml.genomes[it.key]
            checkNotNull(deserializedMap) { "Annotations not defined for $it" }
            fromMap(it.key, deserializedMap)
        }, yaml)
    }

    fun saveYaml(yamlPath: Path, map: Map<String, GenomeAnnotationsConfig>) {
        val convertedMap = LinkedHashMap<String, LinkedHashMap<String, Any>>()
        map.forEach { (build, genomeAnnotationsConfig) ->
            convertedMap[build] = toMap(build, genomeAnnotationsConfig)
        }
        val toSerialize = AnnotationConfigYAMLSerializable().apply {
            version = VERSION
            genomes = convertedMap
        }
        AnnotationConfigYAMLSerializable.save(yamlPath, toSerialize)
    }


    private fun loadOrCreateAnnotationsConfig(yamlConfig: Path): Map<String, GenomeAnnotationsConfig> {
        if (yamlConfig.exists) {
            // check if outdated
            val (version, mapping, yaml) = parseYaml(yamlConfig, VERSION)
            if (mapping != null) {
                return mapping
            }

            val dateFormat = SimpleDateFormat("yyyy-MM-dd")
            val hash = yaml.genomes.hashCode().toString().sha
            val suffix = "${dateFormat.format(Date())}.$hash"

            val bkPath = yamlConfig.withExtension("$suffix.yaml")
            yamlConfig.move(bkPath, StandardCopyOption.ATOMIC_MOVE)
            LOG.error(
                "Outdated annotations file $yamlConfig (version: $version) was backed up to $bkPath " +
                        "and replaced with recent version: $VERSION)."
            )
        }
        // create if not exist or was outdated
        yamlConfig.checkOrRecalculate { output ->
            val text = AnnotationsConfigLoader::class.java
                .getResourceAsStream("/annotations.yaml")!!
                .bufferedReader(Charsets.UTF_8)
                .readText()
            output.path.write(text)
        }
        // ensure version
        val (newVersion, newMapping, _) = parseYaml(yamlConfig, VERSION)
        require(newVersion == VERSION) {
            "Bundled annotations are expected to have latest format version $VERSION but was $newVersion."
        }
        return newMapping!!
    }

    /**
     * This is a technical class that is used for serialization / deserialization in YAML
     */
    class AnnotationConfigYAMLSerializable {
        @JvmField
        var version: Int = 0

        /**
         * A temporary object for loading weakly-typed YAML data.
         */
        @JvmField
        var genomes = LinkedHashMap<String, LinkedHashMap<String, Any>>()

        companion object {
            /**
             * Used during serialization / deserialization procedure
             */
            private fun createObjectMapper(): ObjectMapper {
                val yamlFactory = YAMLFactory()
                yamlFactory.configure(YAMLGenerator.Feature.MINIMIZE_QUOTES, true)
                return ObjectMapper(yamlFactory)
            }

            /**
             * Deserialization. See [saveYaml] for serialization
             */
            fun load(yamlPath: Path): AnnotationConfigYAMLSerializable = yamlPath.bufferedReader().use { reader ->
                createObjectMapper().readValue(reader, AnnotationConfigYAMLSerializable::class.java)
            }!!

            /**
             * Serialization. See [saveYaml] for deserialization
             */
            fun save(yamlPath: Path, serializable: AnnotationConfigYAMLSerializable) {
                yamlPath.bufferedWriter().use { writer ->
                    createObjectMapper().writeValue(writer, serializable)
                }
            }

        }
    }

    private const val UCSC_ALIAS_FIELD = "ucsc_alias"
    private const val SPECIES_FIELD = "species"
    private const val DESCRIPTION_FIELD = "description"
    private const val CHR_ALT_NAME_TO_CANONICAL_FIELD = "chr_alt_name_to_canonical"
    private const val ALIASES_FIELD = "aliases"
    private const val GTF_FIELD = "gtf"
    private const val UCSC_ANNOTATIONS_LEGACY_FIELD = "ucsc_annotations_legacy"
    private const val SEQUENCE_FIELD = "sequence"
    private const val CHROMSIZES_FIELD = "chromsizes"
    private const val REPEATS_FIELD = "repeats"
    private const val CYTOBANDS_FIELD = "cytobands"
    private const val GAPS_FIELD = "gaps"
    private const val CENTROMERES_FIELD = "centromeres"
    private const val CGIS_FIELD = "cgis"
    private const val BIOMART_FIELD = "biomart"
    private const val BIOMART_URL_FIELD = "url"
    private const val BIOMART_DATASET_FIELD = "dataset"

    fun fromMap(build: String, deserializedMap: Map<String, Any>): GenomeAnnotationsConfig {
        /** Formats fix: example below deserializes into list of maps
         * chr_alt_name_to_canonical:
         *    MT: chrM
         */
        val chrAltName2CanonicalMapping = if (CHR_ALT_NAME_TO_CANONICAL_FIELD in deserializedMap) {
            (deserializedMap[CHR_ALT_NAME_TO_CANONICAL_FIELD] as List<*>).associate {
                (it as Map<*, *>).entries.first().let { (k, v) -> k as String to v as String}
            }
        } else {
            emptyMap()
        }

        val ucscAlias = deserializedMap[UCSC_ALIAS_FIELD] as? String
        /* names = build + ucsc_alias + aliases, with duplicates removed later */
        val names = arrayListOf(build)
        ucscAlias?.let { names.add(it) }
        names.addAll(when (val aliases = deserializedMap[ALIASES_FIELD]) {
            is String -> listOf(aliases)
            is List<*> -> aliases.map { it.toString() }
            else -> emptyList()
        })

        val biomart = deserializedMap[BIOMART_FIELD]
        val mart = if (biomart != null) {
            val biomartMap = biomart as Map<*, *>
            val url = biomartMap[BIOMART_URL_FIELD] as String
            val dataset = biomartMap[BIOMART_DATASET_FIELD] as String
            Biomart(dataset, url)
        } else {
            null
        }
        try {
            return GenomeAnnotationsConfig(
                deserializedMap[SPECIES_FIELD] as String,
                ucscAlias,
                names.distinct(),
                deserializedMap[DESCRIPTION_FIELD] as String,
                deserializedMap[GTF_FIELD] as String,
                chrAltName2CanonicalMapping,
                deserializedMap[UCSC_ANNOTATIONS_LEGACY_FIELD]?.toString()?.toBoolean() ?: false,
                deserializedMap[SEQUENCE_FIELD] as String,
                deserializedMap[CHROMSIZES_FIELD] as String,
                deserializedMap[REPEATS_FIELD] as String?,
                deserializedMap[CYTOBANDS_FIELD] as String?,
                deserializedMap[GAPS_FIELD] as String?,
                deserializedMap[CENTROMERES_FIELD] as String?,
                deserializedMap[CGIS_FIELD] as String?,
                mart
            )
        } catch (e: Exception) {
            throw RuntimeException("Cannot parse: [${deserializedMap}]", e)
        }
    }

    fun toMap(build: String, genomeAnnotationsConfig: GenomeAnnotationsConfig): LinkedHashMap<String, Any> {
        val aliases = genomeAnnotationsConfig.names.filter {
            it != build && it != genomeAnnotationsConfig.ucscAlias
        }
        val result = linkedMapOf(
            SPECIES_FIELD to genomeAnnotationsConfig.species,
            DESCRIPTION_FIELD to genomeAnnotationsConfig.description,
            // Should be serialized as list or single value (in most cases)
            ALIASES_FIELD to if (aliases.size == 1) aliases.first() else aliases,
            GTF_FIELD to genomeAnnotationsConfig.gtfUrl,
            SEQUENCE_FIELD to genomeAnnotationsConfig.sequenceUrl,
            CHROMSIZES_FIELD to genomeAnnotationsConfig.chromsizesUrl
        )
        if (genomeAnnotationsConfig.gapsUrl != null) {
            result[GAPS_FIELD] = genomeAnnotationsConfig.gapsUrl
        }
        if (genomeAnnotationsConfig.ucscAlias != null) {
            result[UCSC_ALIAS_FIELD] = genomeAnnotationsConfig.ucscAlias
        }
        if (genomeAnnotationsConfig.chrAltName2CanonicalMapping.isNotEmpty()) {
            // Transform Map<String, String> to List<Map<String, String>>, see example while deserializing
            result[CHR_ALT_NAME_TO_CANONICAL_FIELD] =
                genomeAnnotationsConfig.chrAltName2CanonicalMapping.entries.map { (k, v) -> mapOf(k to v) }
        }
        if (genomeAnnotationsConfig.ucscAnnLegacyFormat) {
            result[UCSC_ANNOTATIONS_LEGACY_FIELD] = true
        }
        if (genomeAnnotationsConfig.repeatsUrl != null) {
            result[REPEATS_FIELD] = genomeAnnotationsConfig.repeatsUrl
        }
        if (genomeAnnotationsConfig.cytobandsUrl != null) {
            result[CYTOBANDS_FIELD] = genomeAnnotationsConfig.cytobandsUrl
        }
        if (genomeAnnotationsConfig.centromeresUrl != null) {
            result[CENTROMERES_FIELD] = genomeAnnotationsConfig.centromeresUrl
        }
        if (genomeAnnotationsConfig.cpgIslandsUrl != null) {
            result[CGIS_FIELD] = genomeAnnotationsConfig.cpgIslandsUrl
        }
        if (genomeAnnotationsConfig.mart != null) {
            result[BIOMART_FIELD] = linkedMapOf(
                BIOMART_DATASET_FIELD to genomeAnnotationsConfig.mart.dataset,
                BIOMART_URL_FIELD to genomeAnnotationsConfig.mart.url
            )
        }
        return result
    }


}
