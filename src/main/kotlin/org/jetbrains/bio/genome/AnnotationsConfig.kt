package org.jetbrains.bio.genome

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.text.SimpleDateFormat
import java.util.*
import java.util.zip.GZIPInputStream

/**
 * Generic genome information.
 *
 * @param species Species name
 * @param ucscAlias UCSC name for the build
 * @param names All known build names (release name and other aliases, including the UCSC alias) as ordered set
 * @param description About Genome build
 * @param gtfUrl GTF genes annotations url
 * @param chrAltName2CanonicalMapping Mapping of alternative chr names to canonical names from *.chrom.sizes. E.g. MT -> chrM for GTF parsing with UCSC genomes
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
        val description: String?,
        val gtfUrl: String,
        val chrAltName2CanonicalMapping: Map<String, String>,
        val ucscAnnLegacyFormat: Boolean,
        val sequenceUrl: String,
        val chromsizesUrl: String,
        val repeatsUrl: String?,
        val cytobandsUrl: String?,
        val gapsUrl: String,
        val centromeresUrl: String?,
        val cpgIslandsUrl: String?,
        val mart: Biomart?
)

object AnnotationsConfigLoader {
    const val VERSION: Int = 4
    private val LOG = LoggerFactory.getLogger(AnnotationsConfigLoader::class.java)

    @Volatile
    private var pathAndYamlConfig: Pair<Path, Map<String, GenomeAnnotationsConfig>>? = null

    /**
     * Parse YAML file and load configuration. If YAML file doesn't exist if would be created with
     * default bundled content using [bioinf-commons/src/main/resources/annotations.yaml]
     */
    fun init(yamlConfig: Path) {
        // 'Double checked synchronization' to configure pathAndYamlConfig singleton
        if (pathAndYamlConfig == null) {
            synchronized(this) {
                // if 2-nd thread enters section after 1-st thread already left it: do not reassign the value
                if (pathAndYamlConfig == null) {
                    pathAndYamlConfig = yamlConfig to loadOrCreateBuild2ConfigMapping(yamlConfig)
                }
            }
        }
        val path = pathAndYamlConfig!!.first
        require(yamlConfig == path) {
            "Already initialized with $path, couldn't be overridden with $yamlConfig."
        }
    }

    val initialized: Boolean
        get() = pathAndYamlConfig != null

    /**
     * Ensure annotations config loaded, use [AnnotationsConfigLoader.init]
     */
    val yamlConfig: Path
        get() {
            if (pathAndYamlConfig == null) {
                error("Annotations config not initialized, call [AnnotationsConfig.init] first")
            }
            return pathAndYamlConfig!!.first
        }
    /**
     * Ensure annotations config loaded, use [AnnotationsConfigLoader.init]
     */
    val builds: Set<String>
        get() {
            if (pathAndYamlConfig == null) {
                error("Annotations config not initialized, call [AnnotationsConfig.init] first")
            }
            return pathAndYamlConfig!!.second.keys
        }

    /**
     * Ensure annotations config loaded, use [AnnotationsConfigLoader.init]
     */
    operator fun get(build: String): GenomeAnnotationsConfig {
        if (pathAndYamlConfig == null) {
            error("Annotations config not initialized, call [AnnotationsConfig.ensureInitialized] first")
        }

        check(build in builds) {
            "Unexpected build name $build, annotations are available only for ${builds.joinToString()}"
        }
        return pathAndYamlConfig!!.second.getValue(build)
    }

    internal fun parseYaml(yamlConfig: Path, currentVersion: Int): Triple<Int, Map<String, GenomeAnnotationsConfig>?, YamlMapper> {
        val yaml = YamlMapper.load(yamlConfig)
        val yamlVers = yaml.version
        if (yamlVers != currentVersion) {
            return Triple(yamlVers, null, yaml)
        }

        val buildsConfig = yaml.genomes.keys.map { build ->
            build to yaml[build]
        }.toMap()

        return Triple(yamlVers, Collections.unmodifiableMap(buildsConfig), yaml)
    }

    private fun loadOrCreateBuild2ConfigMapping(yamlConfig: Path): Map<String, GenomeAnnotationsConfig> {
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
            LOG.info("Outdated annotations file $yamlConfig (version: $version) was backed up to $bkPath " +
                    "and replaced with recent version: $VERSION).")
        }
        // create if not exist or was outdated
        yamlConfig.checkOrRecalculate { output ->
            val text = AnnotationsConfigLoader::class.java
                    .getResourceAsStream("/annotations.yaml")
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

            checkNotNull(genomeAttrs) { "Annotations not defined for $build" }

            val chrAltName2CanonicalMapping = if ("chr_alt_name_to_canonical" in genomeAttrs) {
                (genomeAttrs["chr_alt_name_to_canonical"] as List<Map<String, String>>).map {
                    it.entries.first().let { (k, v) -> k to v }
                }.toMap()
            } else {
                emptyMap()
            }

            /* names = build + ucsc_alias + aliases, with duplicates removed later */
            val names = arrayListOf(build)
            val ucscAlias = genomeAttrs["ucsc_alias"] as? String
            ucscAlias?.let { names.add(it) }
            names.addAll(when (val aliases = genomeAttrs["aliases"]) {
                is String -> listOf(aliases)
                is List<*> -> aliases.map { it.toString() }
                else -> emptyList()
            })

            val biomart = genomeAttrs["biomart"]
            val mart = if (biomart != null) {
                val biomartMap = biomart as Map<String, String>
                val url = biomartMap["url"] as String
                val dataset = biomartMap["dataset"] as String
                Biomart(dataset, url)
            } else {
                null
            }
            try {
                return GenomeAnnotationsConfig(
                    genomeAttrs["species"] as String,
                    ucscAlias,
                    names.distinct(),
                    genomeAttrs["description"] as String,
                    genomeAttrs["gtf"] as String,
                    chrAltName2CanonicalMapping,
                    genomeAttrs["ucsc_annotations_legacy"]?.toString()?.toBoolean() ?: false,
                    genomeAttrs["sequence"] as String,
                    genomeAttrs["chromsizes"] as String,
                    genomeAttrs["repeats"] as String?,
                    genomeAttrs["cytobands"] as String?,
                    genomeAttrs["gaps"] as String,
                    genomeAttrs["centromeres"] as String?,
                    genomeAttrs["cgis"] as String?,
                    mart
                )
            } catch (e: Exception) {
                throw RuntimeException("Cannot parse: [${genomeAttrs}]", e)
            }
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

    fun updateAnnotations(markupPath: Path) {
        val path = markupPath / "annotations.tar.gz.txt"
        if (path.exists) {
            var version = -1
            path.bufferedReader().use {
                try {
                    version = it.readLine().toInt()
                } catch (t: Throwable) {
                    LOG.error("Failed to read version number from $path", t)
                }
            }
            if (version != VERSION) {
                LOG.info("Outdated cache file $path (version: $version) recent version: $VERSION).")
                path.delete()
            }
        }
        val bundled = AnnotationsConfigLoader::class.java.getResource("/annotations.tar.gz")
        checkNotNull(bundled) { "Bundled annotations not available" }

        if (path.notExists) {
            Tar.decompress(GZIPInputStream(bundled.openStream()), markupPath.toFile(), false)
            path.bufferedWriter().use {
                it.write(VERSION.toString())
            }
        }
    }

}
