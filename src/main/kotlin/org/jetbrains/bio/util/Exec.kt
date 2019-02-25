package org.jetbrains.bio.util

import org.apache.log4j.Logger
import java.nio.file.Files
import java.nio.file.Path

/**
 * A sweeter extended version of the [Exec] API.
 *
 * @param log if `true` *stdout* and *stderr* will be logged,
 *            otherwise these will be redirect to [System.out]
 *            and [System.err] of the current process.
 * @param init a custom initializer for [ProcessBuilder].
 */
fun String.exec(vararg args: Any, log: Boolean = true,
                init: ProcessBuilder.() -> Unit = {}) =
        Exec.exec(this, *args,
                output = if (log) OutputType.LOG else OutputType.IGNORE,
                init = init)


fun Path.exec(vararg args: Any, log: Boolean = true,
              init: ProcessBuilder.() -> Unit = {}) = this.toString().exec(*args, log = log, init = init)

enum class OutputType {
    LOG,
    TEXT,
    IGNORE
}

private object Exec {
    private val LOG = Logger.getLogger(Exec::class.java)

    internal fun exec(executable: String,
                      vararg args: Any,
                      output: OutputType,
                      init: ProcessBuilder.() -> Unit = {}): String? {
        val pb = ProcessBuilder(executable, *args.map(Any::toString).toTypedArray())
        when (output) {
            OutputType.IGNORE -> pb.inheritIO()
            else -> pb.redirectErrorStream(true)
        }

        pb.init()

        LOG.info(pb.command().joinToString(" "))

        val process = pb.start()

        if (output == OutputType.LOG) {
            process.inputStream.bufferedReader().use {
                for (line in it.lineSequence()) {
                    LOG.info(line)
                }
            }
        }

        val resultText = when (output) {
            OutputType.TEXT -> process.inputStream.bufferedReader().readText()
            else -> null
        }

        val exitCode = process.waitFor()
        if (exitCode != 0) {
            throw RuntimeException("Process stopped with exit code $exitCode.\n" +
                    "Output\n${process.inputStream.bufferedReader().readText()}\n" +
                    "Error\n${process.errorStream.bufferedReader().readText()}")
        }

        return resultText
    }
}

/**
 * A helper for fetching package resources in a safe manner.
 *
 * Example:
 *
 *     private class Proxy
 *
 *     withResource(Proxy::class.java, ...) { ... }
 *
 * @param proxy class used for fetching resources
 * @param name path to resource without the starting '/', e.g. "run_sleuth.R".
 * @param block DWYW.
 */
inline fun <T> withResource(proxy: Class<*>, name: String,
                            checkNotNull: Boolean = true,
                            block: (Path) -> T): T {
    val resource = proxy.getResource("/$name")
    if (checkNotNull) {
        checkNotNull(resource) { "Resource '$name' not found for class: $proxy" }
    } else {
        Logger.getRootLogger().warn("Resource '$name' not found for class: $proxy")
    }

    return withTempDirectory(proxy.simpleName) { dir ->
        val path = dir / name
        if (resource != null) {
            path.parent.createDirectories()
            Files.newOutputStream(path).use {
                resource.openStream().copyTo(it)
            }

        }
        block(path)
    }
}
