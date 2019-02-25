package org.jetbrains.bio.util

const val JOPTSIMPLE_SUPPRESS_EXIT = "joptsimple.suppressExit"
fun checkOrFail(condition: Boolean, block: () -> String) {
    if (!condition) {
        System.err.println("ERROR: ${block()}")
        val suppressExit = System.getProperty(JOPTSIMPLE_SUPPRESS_EXIT)
        if (suppressExit?.toBoolean() != true) {
            System.exit(1)
        }
    }
}
