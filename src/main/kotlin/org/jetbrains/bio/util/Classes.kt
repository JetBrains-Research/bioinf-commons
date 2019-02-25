package org.jetbrains.bio.util

import org.apache.log4j.Logger
import java.net.URL
import java.net.URLClassLoader
import java.nio.file.Path
import java.util.function.Consumer
import java.util.jar.JarFile

class ClassPathEntry(private val path: Path) {

    fun process(consumer: Consumer<String>) {
        if (path.notExists) {
            return
        }
        if (path.isDirectory) {
            processFolder("", path, consumer)
        } else {
            processJar(path, consumer)
        }
    }

    companion object {
        private val CLASS_FILE_SUFFIX = ".class"

        private fun processFolder(curPath: String, parent: Path, consumer: Consumer<String>) {
            val prefix = if (curPath.isEmpty()) "" else "$curPath."
                for (f in parent.list()) {
                    if (f.name.endsWith(CLASS_FILE_SUFFIX)) {
                        consumer.accept(prefix + f.name.substringBefore(CLASS_FILE_SUFFIX))
                    } else if (f.isDirectory) {
                        processFolder(prefix + f.name, f, consumer)
                    }
            }
        }

        private fun processJar(path: Path, consumer: Consumer<String>) {
            JarFile(path.toFile()).use { jarFile ->
                val entries = jarFile.entries()
                while (entries.hasMoreElements()) {
                    val jarEntry = entries.nextElement()
                    val name = jarEntry.name
                    if (name.endsWith(CLASS_FILE_SUFFIX)) {
                        consumer.accept(name.substringBefore(CLASS_FILE_SUFFIX).replace("/|\\\\".toRegex(), "."))
                    }
                }
            }
        }
    }
}

object ClassProcessor {
    private val LOG = Logger.getLogger(ClassProcessor::class.java)
    private val classLoaders = hashSetOf<ClassLoader>(ClassLoader.getSystemClassLoader())

    /**
     * This method tries to initialize given class with constructor(), getInstance() or INSTANCE or ourInstance fields and methods
     * @return Object instance if succeeded, null otherwise
     */
    @JvmStatic
    fun tryToInstantiate(clazzz: Class<*>): Any? {
        try {
            return clazzz.newInstance()
        } catch (e: Exception) {
            // Constructor is private or exception in it
            LOG.debug(e)
        }
        try {
            val field = clazzz.getField("INSTANCE")
            if (field != null) {
                return field.get(clazzz)
            }
        } catch (e: Exception) {
            // Cannot call getInstance
            LOG.debug(e)
        }
        LOG.warn("Failed to create instance for class: ${clazzz.name}")
        return null
    }

    @JvmStatic
    fun processClasses(consumer: Consumer<String>) {
        classPaths().forEach { ClassPathEntry(it.toURI().toPath()).process(consumer) }
    }

    fun classPaths(): List<URL> {
        return classLoaders.flatMap { classLoader ->
            (classLoader as URLClassLoader).urLs.filter { it.protocol == "file" }.map { it }
        }
    }
}

