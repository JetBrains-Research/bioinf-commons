package org.jetbrains.bio.genome.data

import com.fasterxml.jackson.annotation.JsonIgnoreProperties
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.GenomeQuery.Companion.parseGenomeQueryId
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.io.Reader
import java.io.StringWriter
import java.io.Writer
import java.nio.file.Path
import java.util.*

@JsonIgnoreProperties("id", "genomeQuery", "tracksMap")
class DataConfig {

    companion object {
        val LOG: Logger = LoggerFactory.getLogger(DataConfig::class.java)

        const val TECHNICAL_REPLICATE_PATTERN = "r\\d+"

        const val FORMAT = """YAML configuration for biological data.
---
Genome examples:
- mm10
- hg19[chr1,chr2,chr3]
---
Tracks:
Each condition is allowed to have multiple replicates. Replicates
can be either implicitly labeled by their position within the
condition or have explicit human-readable labels as in the example below.

With explicit labels:
    <condition>:
        <replicate>: path/to/replicate/data
        <replicate>: path/to/replicate/data

With explicit labels and meta information:
    <condition>:
        <replicate>:
            path: path/to/replicate/data
            failed: false
        <replicate>:
            path: path/to/replicate/data
            failed: true
            comment: failed

Without labels:
    <condition>:
    - path/to/replicate/data
    - path/to/replicate/data

---
Aux:
A data-type-specific information map can be provided for each data type.
The format is:
aux:
  <data type>:
    <key1>: <value1>
    <key2>: <value2>
"""

        private fun createMapper(): ObjectMapper {
            val yamlFactory = YAMLFactory()
            yamlFactory.configure(YAMLGenerator.Feature.MINIMIZE_QUOTES, true)
            return ObjectMapper(yamlFactory)
        }

        /** Loads configuration from [Reader]. */
        fun load(reader: Reader, id: String): DataConfig =
            createMapper().readValue(reader, DataConfig::class.java).apply {
                this.id = id
                // Check genome correctness
                genomeQuery
                // Force lazy loading to check structure correctness
                tracksMap
            }


        /** Loads configuration from a YAML file with id as file name. */
        fun load(path: Path) = path.bufferedReader().use { load(it, path.stem) }
    }

    /**
     * Human-readable identifier of the dataset.
     */
    var id: String = "unknown"

    @JvmField
    var genome = "unknown"

    /**
     * A temporary object for loading weakly-typed YAML data.
     */
    @JvmField
    var tracks = LinkedHashMap<String, LinkedHashMap<String, Any>>()

    @JvmField
    var aux = LinkedHashMap<String, LinkedHashMap<String, Any>>()

    /**
     * A genome query, which specifies genome build and chromosome restriction.
     * NOTE: we don't use lateinit var here to guarantee read-only access
     */
    val genomeQuery: GenomeQuery by lazy {
        parseGenomeQueryId(genome).let { (build, chromosomes)
            ->
            GenomeQuery(Genome[build], *chromosomes)
        }
    }

    /**
     * Configured tracks.
     * NOTE: we don't use lateinit var here to guarantee read-only access
     */
    val tracksMap: LinkedHashMap<Key, Section> by lazy {
        val map = LinkedHashMap<Key, Section>()
        for ((dataType, inner) in tracks) {
            for ((condition, replicates) in inner.mapKeys { Cell(it.key) }) {
                val section = when (replicates) {
                    /** Implicitly replicated section. */
                    is List<*> ->
                        replicates.mapIndexed { index, r -> "r$index" to ReplicateData.of(r) }
                    /** A section with per-replicate labels. */
                    is Map<*, *> ->
                        replicates.map { (r, v) -> r.toString() to ReplicateData.of(v) }
                    else -> throw IllegalStateException("Unknown replicates structure: $replicates")
                }

                map[Key(dataType, condition)] = section
            }
        }
        return@lazy map
    }

    /** Saves configuration in a YAML file. */
    fun save(writer: Writer) {
        createMapper().writeValue(writer, this)
    }


    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is DataConfig) return false
        if (id != other.id) return false
        if (genome != other.genome) return false
        if (tracks != other.tracks) return false
        return true
    }

    override fun hashCode(): Int {
        return Objects.hash(id, genome, tracks)
    }

    override fun toString(): String {
        val writer = StringWriter()
        save(writer)
        return writer.toString()
    }

    fun cells() = tracksMap.keys.map(Key::cellId).distinct()

    fun dataTypes() = tracksMap.keys.map(Key::dataType).distinct()

    fun auxDataTypes() = aux.keys.distinct()

    fun getCellIds(dataType: DataType): Collection<Cell> =
        tracksMap.keys.filter { it.dataType.toDataType() == dataType }.map(Key::cellId).distinct()

    data class Key(val dataType: String, val cellId: Cell)

    operator fun <T> get(dataType: String, key: ReplicateDataKey<T>): T =
        key.parse(aux[dataType]?.get(key.name))

}

data class ReplicateDataKey<T>(val name: String, private val parser: (Any?) -> T) {
    fun parse(value: Any?): T = parser(value)
}

data class ReplicateData constructor(val path: Path, internal val meta: Map<String, Any> = emptyMap()) {

    init {
        check(path.exists && path.isReadable) {
            "Failed to access file: $path"
        }
    }

    val failedTrack: Boolean
        get() = this[FAILED_TRACK]

    operator fun <T> get(key: ReplicateDataKey<T>): T = key.parse(meta[key.name])

    companion object {

        private const val PATH_KEY = "path"
        val FAILED_TRACK = ReplicateDataKey("failed") {
            it == true
        }

        fun of(any: Any?): ReplicateData = when (any) {
            is String -> ReplicateData(any.toPath())
            is Map<*, *> -> {
                check(PATH_KEY in any) {
                    "Key $PATH_KEY is required"
                }
                ReplicateData(any[PATH_KEY]!!.toString().toPath(), any as Map<String, Any>)
            }
            else ->
                throw IllegalArgumentException("Unknown replicates structure: $any")
        }
    }
}

typealias Section = List<Pair<String, ReplicateData>>
