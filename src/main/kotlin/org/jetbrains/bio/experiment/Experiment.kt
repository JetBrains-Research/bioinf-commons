package org.jetbrains.bio.experiment

import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.resolve
import org.jetbrains.bio.util.time
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Experiment is a named computation procedure with description and fixed location.
 * Main entry point is [run] method.
 * Any experiment class should implements method [doCalculations], which is being executed.
 *
 * @author Oleg Shpynov
 */
abstract class Experiment @JvmOverloads constructor(
    /** Folder, the data produced by this experiment should be stored. */
    open val experimentFolder: String,
    /** Prefix of the log-file. */
    experimentName: String? = null
) {

    // XXX we don't do this in the constructor because `.javaClass` is
    // not available there.
    val name: String = experimentName ?: javaClass.simpleName

    /** A human-readable description of this experiment. */
    var description: String? = null

    val experimentPath: Path get() = Configuration.experimentsPath / experimentFolder

    /**
     * Main entry point of each experiment, called from the [run] method.
     */
    @Throws(Exception::class)
    protected abstract fun doCalculations()

    // Please use the property above in Kotlin.
    fun getExperimentPath(vararg chunks: String): Path {
        return Configuration.experimentsPath.resolve(experimentFolder, *chunks)
    }

    fun run() {
        if (description != null) {
            LOG.info("Description:\n$description\n")
        }
        LOG.time(message = "Run experiment $experimentFolder for $name...") {
            experimentPath.createDirectories()
            doCalculations()
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(Experiment::class.java)
    }
}
