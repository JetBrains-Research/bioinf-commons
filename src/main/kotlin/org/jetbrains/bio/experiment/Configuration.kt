@file:Suppress("unused")

package org.jetbrains.bio.experiment

import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.io.IOException
import java.nio.file.Path
import java.util.*
import kotlin.reflect.KProperty

/**
 * A wrapper around system-wide properties.
 *
 * Properties can be configured by Java command line options
 *   -D<property_name=...
 * or in batch mode by providing path to properties file
 *   -Dconfig.path=...
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 */

object Configuration {


    /**
     * .properties file can be configured using this id
     */
    const val CONFIG_PATH_PROPERTY = "config.path"

    /**
     * These properties will override values provided in .properties file
     */
    private const val GENOME_PATHS_PROPERTY = "genomes.path"
    private const val RAW_DATA_PATH_PROPERTY = "raw.data.path"
    private const val EXPERIMENTS_PATH_PROPERTY = "experiments.path"

    /** Path to genome-specific annotation, e.g. chrom.sizes file. */
    var genomesPath: Path by PropertyPathDelegate(GENOME_PATHS_PROPERTY)

    val genomesPathConfigured: Boolean
        get() = System.getProperties().containsKey(GENOME_PATHS_PROPERTY)


    /** Path to raw data, e.g. raw reads in FASTQ format. */
    var rawDataPath: Path by PropertyPathDelegate(RAW_DATA_PATH_PROPERTY)

    val rawDataPathConfigured: Boolean
        get() = System.getProperties().containsKey(RAW_DATA_PATH_PROPERTY)

    /** Path to the work/experiments folder. */
    var experimentsPath: Path by PropertyPathDelegate(EXPERIMENTS_PATH_PROPERTY)

    val experimentsPathConfigured: Boolean
        get() = System.getProperties().containsKey(EXPERIMENTS_PATH_PROPERTY)

    /** Path to cache data. */
    val cachePath: Path
        get() = experimentsPath / "cache"

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
            val configPath = System.getProperty(CONFIG_PATH_PROPERTY, null)
                    ?.replace("~", properties.getProperty("user.home"))?.toPath()
            if (configPath != null && configPath.exists) {
                LOG.info("Loading config: '$configPath'...")
                if (!configPath.isReadable) {
                    LOG.error("Cannot read: '$configPath'.")
                } else {
                    LOG.info("""
                        |----------------  config ----------------
                        |${configPath.bufferedReader().readText()}
                        |-----------------------------------------
                        |""".trimMargin()
                    )
                    val loadedProperties = Properties()
                    try {
                        loadedProperties.load(configPath.inputStream())
                    } catch (e: IOException) {
                        LOG.error("Error while loading $configPath.", e)
                    }

                    // Merge with current runtime properties, but runtime overrides defaults
                    for ((key, value) in loadedProperties) {
                        if (!properties.containsKey(key)) {
                            properties[key] = value
                        }
                    }
                }
            }

            // Configure fields
            if (properties.containsKey(GENOME_PATHS_PROPERTY)) {
                genomesPath = properties.getPath(GENOME_PATHS_PROPERTY)
            }
            if (properties.containsKey(RAW_DATA_PATH_PROPERTY)) {
                rawDataPath = properties.getPath(RAW_DATA_PATH_PROPERTY)
            }
            if (properties.containsKey(EXPERIMENTS_PATH_PROPERTY)) {
                experimentsPath = properties.getPath(EXPERIMENTS_PATH_PROPERTY)
            }
        } finally {
            // Mark as initialized even if something went wrong
            initialized = true
        }
    }

    /**
     * Property with the following invariant:
     * 1. all write access allowed before the first read access
     * 2. single write access allowed
     * 3. check initialized before read access
     */
    class PropertyPathDelegate(private val propertyName: String) {
        var field: Path? = null

        operator fun getValue(thisRef: Any?, prop: KProperty<*>): Path {
            synchronized(Configuration) {
                initialize()
                requireNotNull(field) {
                    "Path '${prop.name}' not initialized. " +
                            "Use -D$propertyName= or -D$CONFIG_PATH_PROPERTY="
                }
                return field!!
            }
        }

        operator fun setValue(thisReef: Any?, prop: KProperty<*>, value: Path) {
            synchronized(Configuration) {
                require(!initialized || field == value) {
                    "Path '${prop.name}' already initialized. Cannot change '$field' to '$value'"
                }
                // Perform set only if isn't overridden, e.g. if is null:
                if (field == null) {
                    field = value
                }
            }
        }
    }

    private val LOG = LoggerFactory.getLogger(Configuration::class.java)
}

private fun Properties.getPath(name: String): Path {
    val path = getProperty(name).trim().replace("~", getProperty("user.home")).toPath()
    path.createDirectories()
    check(path.isDirectory) {
        "$name is not a directory or does not exist: '$path'"
    }
    return path
}
