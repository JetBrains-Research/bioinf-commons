package org.jetbrains.bio.util

import joptsimple.*
import org.jetbrains.bio.io.BedFormat
import org.jline.terminal.TerminalBuilder
import java.io.StringWriter
import java.nio.file.Path
import java.util.*
import kotlin.math.max

/**
 * NOTE:
 * Do not implement get[], because is has ambiguous meaning: [OptionSet#valueOf] or [OptionSet#valuesOf]
 */
operator fun OptionSet.contains(option: String): Boolean = has(option)

/**
 * @param args Cmdline arguments
 * @param description Tool description
 * @param acceptNonOptionArguments Pass true in order to allow positional arguments 
 * @param block Additional cmdline args description
 */
fun OptionParser.parse(
        args: Array<String>,
        description: String? = null,
        acceptNonOptionArguments : Boolean = false,
        block: (OptionSet) -> Unit
) {
    try {
        formatHelpWith(BuiltinHelpFormatter(
                max(80, try {
                    TerminalBuilder.terminal().width
                } catch (e: Exception) {
                    System.err.println("Warning: Cannot defined terminal width: ${e.message}")
                    0
                }),
                2
        ))

        acceptsAll(listOf("h", "?", "help"), "Show help").forHelp()

        val options = parse(*args)
        if ("help" in options) {
            if (description != null) {
                System.err.println(description)
                System.err.println("")
                System.err.println("")
            }
            System.err.print("Arguments: ")
            System.err.println(Arrays.toString(args))
            printHelpOn(System.err)
            val suppressExit = System.getProperty(JOPTSIMPLE_SUPPRESS_EXIT)
            if (suppressExit == null || !suppressExit.toBoolean()) {
                System.exit(0)
            }
        }

        if (!acceptNonOptionArguments && options.nonOptionArguments().isNotEmpty()) {
            fail("Unrecognized options: ${options.nonOptionArguments()}", args)
        }

        block(options)
    } catch (e: Throwable) {
        // Motivation:
        // 1. Exceptions:
        //    * OptionException: is expected exception when checking given arguments, validates it's value range,
        //      type, etc, not stack trace should be shown to user. So show error + help about expected argument
        //      values.
        //
        //    * In case of uncaught exception - log it with stacktrace. Optionally show help - it could helpful,
        //      I'm not insisting on printing Help.
        //
        //      Also we could just rethrow exception here, but it will have discussable side effect. At the moment we
        //      cannot test, that main() fails with expected uncaught exception and actually we shouldn't do it.
        //      The only choice is to check that error was handled and logged properly. On the other hand we cannot
        //      test that in some scenario exception doesn't occur in Foo.main(args) because it was caught and
        //      logged properly.
        //
        // 2. Block could throw OptionException, not only parse() method
        //
        // 3. System.exit(1) in checkOrFail() is used in order not to return "0" exist code in case of error.
        //
        // 4. See `OptionParserExtensionsTest` and `joptsimple.suppressExit` (JOPTSIMPLE_SUPRRESS_EXIT)
        //    sys property in checkOrFail() method which is useful in method checking CLI behaviour.
        if (e is OptionException) {
            if (e.cause != null && e.cause is ValueConversionException) {
                fail(e.cause!!.message!!, args)
            } else {
                fail(e.message!!, args)
            }
        } else {
            e.printStackTrace()
            fail(e.cause?.message ?: e.message!!, args)
        }
    }
}

fun OptionParser.fail(message: String, args: Array<String>? = null) {
    val help = StringWriter()

    if (args != null) {
        System.err.print("Arguments: ")
        System.err.println(Arrays.toString(args))
    }
    printHelpOn(help)
    checkOrFail(false) { "$message\n$help" }
}

/**
 * A converter for files or directories.
 *
 * @author Sergei Lebedev
 */
abstract class PathConverter : ValueConverter<Path> {
    @Throws(ValueConversionException::class)
    abstract fun check(path: Path)

    @Throws(ValueConversionException::class)
    override fun convert(value: String): Path {
        val path = value.toPath().toAbsolutePath()
        check(path)
        return path
    }

    override fun valueType() = Path::class.java

    override fun valuePattern(): String? = null

    companion object {
        fun exists(ext: String? = null): PathConverter = object : PathConverter() {
            @Throws(ValueConversionException::class)
            override fun check(path: Path) {
                if (path.notExists) {
                    throw ValueConversionException("Path $path does not exists")
                }
                if (ext != null && path.extension.toLowerCase() != ext.toLowerCase()) {
                    throw ValueConversionException("Expected *.$ext file, but was ${path.fileName}")
                }
            }
        }

        fun noCheck(ext: String? = null): PathConverter = object : PathConverter() {
            @Throws(ValueConversionException::class)
            override fun check(path: Path) {
                if (ext != null && path.extension.toLowerCase() != ext.toLowerCase()) {
                    throw ValueConversionException("Expected *.$ext file, but was ${path.fileName}")
                }
            }
        }

        fun notExists(ext: String? = null): PathConverter = object : PathConverter() {
            @Throws(ValueConversionException::class)
            override fun check(path: Path) {
                if (path.exists) {
                    throw ValueConversionException("Path $path already exist")
                }
                if (ext != null && path.extension.toLowerCase() != ext.toLowerCase()) {
                    throw ValueConversionException("Expected *.$ext file, but was ${path.fileName}")
                }
            }
        }

        fun bedtoolsValidFile(ext: String? = null): PathConverter = object : PathConverter() {
            @Throws(ValueConversionException::class)
            override fun check(path: Path) {
                if (path.notExists) {
                    throw ValueConversionException("Path $path does not exist")
                }
                if (ext != null && path.extension.toLowerCase() != ext.toLowerCase()) {
                    throw ValueConversionException("Expected *.$ext file, but was ${path.fileName}")
                }

                val bedFormat = BedFormat.auto(path)
                if (bedFormat.delimiter != '\t') {
                    throw ValueConversionException("Expected TAB separated file, but separator is [${bedFormat.delimiter}]")
                }
                if (bedFormat.fieldsNumber < 3) {
                    throw ValueConversionException("Expected at least first 3 BED fields, but format is [${bedFormat.fmtStr}]")
                }
                if (bedFormat.fieldsNumber < 6 && (bedFormat.extraFieldsNumber != null && bedFormat.extraFieldsNumber > 0)) {
                    throw ValueConversionException("Bedtools will fail recognize strand column, your format is [${bedFormat.fmtStr}]")
                }
            }
        }

    }
}
