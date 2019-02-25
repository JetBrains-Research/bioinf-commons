package org.jetbrains.bio.experiment

import org.apache.log4j.Logger
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.util.*
import java.nio.file.Path

abstract class DataConfigExperiment(folder: String, val configuration: DataConfig) :
        Experiment("configs/${configuration.id}/$folder") {

    init {
        LOG.debug("Check configuration file")
        val configPath = experimentPath / "${configuration.id}.yaml"
        if (!configPath.isDirectory && configPath.isReadable) {
            val oldConfig: DataConfig?
            try {
                oldConfig = DataConfig.load(configPath)
            } catch (t: Throwable) {
                LOG.error("Failed to load config at $configPath", t)
                throw IllegalStateException("Failed to load config at $configPath", t)
            }
            check(configuration == oldConfig) {
                "Config file for $name $configPath changed.\nOLD config:\n$oldConfig\nNEW config:\n$configuration"
            }
            LOG.info("Config file already exists $name: $configPath")
        } else {
            saveConfig(configPath)
        }
    }

    private fun saveConfig(configPath: Path) {
        configPath.bufferedWriter().use { configuration.save(it) }
        LOG.info("Saved config file for $name: $configPath")
    }

    companion object {
        private val LOG = Logger.getLogger(DataConfigExperiment::class.java)

        fun loadDataConfig(input: String, quiet: Boolean = false): DataConfig {
            if (quiet) {
                Logs.quiet()
            }
            // Let's try load from file
            val path = input.toPath()
            if (path.exists && path.isReadable) {
                if (!quiet) {
                    LOG.info("Loading data config from $path")
                }
                return DataConfig.load(path)
            }
            throw IllegalArgumentException("Failed to load dataset for $input")
        }
    }
}
