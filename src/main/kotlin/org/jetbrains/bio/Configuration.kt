package org.jetbrains.bio

import org.jetbrains.bio.util.*
import java.io.IOException
import java.nio.file.Path
import java.util.*
import kotlin.reflect.KProperty

/**
 * A wrapper around system-wide properties.
 *
 * Upon first access [Configuration] loads a property file located
 * at `$HOME/.epigenome/config.properties`.
 *
 * @author Sergei Lebedev
 */

object Configuration {
    const val GENOME_PATHS_PROPERTY = "genomes.path"
    private const val EXPERIMENT_PROPERTY = "experiments.path"
    private const val RAW_DATA_PATH_PROPERTY = "raw.data.path"

    /** Path to genome-specific annotation, e.g. chrom.sizes file. */
    var genomesPath: Path by OverridableLazyInitPathDelegate(this)

    /** Path to raw data, e.g. raw reads in FASTQ format. */
    var rawDataPath: Path by OverridableLazyInitPathDelegate(this)

    /** Path to the data computed by an [Experiment]. */
    var experimentsPath: Path by OverridableLazyInitPathDelegate(this)

    /** Path to GEO data. */
    val geoSamplesPath: Path
        get() = rawDataPath / "geo-samples"

    /** Path to cache data. */
    val cachePath: Path
        get() = experimentsPath / "cache"

    val defaultConfigPath: Path =
            System.getProperty("user.home", "").toPath() / ".epigenome" / "config.properties"

    private val LOCK = Object()

    @Volatile
    private var initialized: Boolean = false
        set(value) {
            require(value) { "Cannot de-initialize config" }
            field = true
        }

    private fun initialize() {
        if (initialized) {
            return
        }
        try {

            val properties = System.getProperties()

            // XXX: Logger normally not initialized at this moment, use System.out instead
            if (defaultConfigPath.exists) {
                println("Loading default config: '$defaultConfigPath'...")

                if (!defaultConfigPath.isReadable) {
                    System.err.println("Cannot read: '$defaultConfigPath'. Using default settings.")
                } else {
                    print("""
                        |----------------  config.properties -----------------
                        |${defaultConfigPath.bufferedReader().readText()}
                        |-----------------------------------------------------
                        |""".trimMargin()
                    )

                    val defaultProperties = Properties()
                    try {
                        defaultProperties.load(defaultConfigPath.inputStream())
                    } catch (e: IOException) {
                        System.err.println("Error while loading $defaultConfigPath: ${e.message}")
                        e.printStackTrace()
                    }

                    // Merge with current runtime properties, but runtime overrides defaults
                    for ((key, value) in defaultProperties) {
                        if (!properties.containsKey(key)) {
                            properties[key] = value
                        }
                    }
                }

            }

            // Init here delegated properties which hasn't been already overridden. To override just
            // invoke setter before initialization, e.g. before read access

            GENOME_PATHS_PROPERTY.let {
                if (properties.containsKey(it)) {
                    genomesPath = Configuration[it]
                }
            }
            RAW_DATA_PATH_PROPERTY.let {
                if (properties.containsKey(it)) {
                    rawDataPath = Configuration[it]
                }
            }
            EXPERIMENT_PROPERTY.let {
                if (properties.containsKey(it)) {
                    experimentsPath = Configuration[it]
                }
            }
        } finally {
            // Fix dirty state after exception, do not allow to override it with clean step
            // so as not to hide error
            initialized = true
        }
    }

    private operator fun get(property: String): Path {
        val value = checkNotNull(System.getProperty(property)) {
            "missing property $property"
        }

        val path = value.trim().toPath()
        path.createDirectories()
        check(path.isDirectory) {
            "$property is not a directory or does not exist: '$path'"
        }
        return path
    }

    fun setExperimentWorkingDir(workDir: Path) {
        Configuration.experimentsPath = workDir

        // In case of missing configuration:
        if (!Configuration.defaultConfigPath.exists) {
            val properties = System.getProperties()
            if (!properties.containsKey(GENOME_PATHS_PROPERTY)) {
                Configuration.genomesPath = workDir
            }
            if (properties.containsKey(RAW_DATA_PATH_PROPERTY)) {
                Configuration.rawDataPath = workDir
            }
        }
    }

    class OverridableLazyInitPathDelegate(val config: Configuration) {
        var field: Path? = null

        operator fun getValue(thisRef: Any?, prop: KProperty<*>): Path {
            synchronized(LOCK) {
                Configuration.initialize()
                requireNotNull(field) { "Path '${prop.name}' not initialized" }
                return field!!
            }
        }

        operator fun setValue(thisReef: Any?, prop: KProperty<*>, value: Path) {
            synchronized(LOCK) {
                require(!initialized) {
                    "Path '${prop.name}' already initialized. Cannot change '$field' to '$value'"
                }
                // Perform set only if isn't overridden, e.g. if is null:
                if (field == null) {
                    field = value
                }
            }
        }
    }
}
