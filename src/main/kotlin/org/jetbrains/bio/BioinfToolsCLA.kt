package org.jetbrains.bio

import org.jetbrains.bio.gse.EnrichmentInRegions
import org.jetbrains.bio.gse.OverlapRegionsWithEachLoci
import org.jetbrains.bio.util.Logs
import org.slf4j.event.Level

object BioinfToolsCLA {
    init {
        // Load build properties
        val resource = BioinfToolsCLA::class.java.getResource("/bioinf.properties")
        resource?.openStream()?.use { System.getProperties().load(it) }
    }

    @JvmStatic
    fun main(args: Array<String>) {
        val helpMsg = helpMessage()

        if (args.isEmpty()) {
            System.err.println("ERROR: No command given.")
            System.err.println(helpMsg)
        } else {
            val arg = args[0]
            when (arg) {
                "-?", "-h", "--help" -> println(helpMsg)
                "-v", "--version" -> println(version())

                else -> {
                    var tool: Tools? = null

                    for (t in Tools.values()) {
                        if (t.command == arg) {
                            tool = t
                            break
                        }
                    }

                    if (tool != null) {
                        val restArgs = args.copyOfRange(1, args.size)
                        tool(restArgs)
                    } else {
                        System.err.println("ERROR: Unknown command: '$arg'")
                        System.err.println(helpMsg)
                    }
                }
            }
        }
    }

    fun version() =
        "${System.getProperty("bioinf.commons.tool.build.version", "@VERSION@.@build@")} " +
                "built on ${System.getProperty("bioinf.commons.tool.build.date", "@DATE@")}"

    private fun helpMessage(): String {
        val sections = ArrayList<Pair<String, String>>()
        sections.add("Option" to "Description")
        sections.add("-------------------------" to "-----------")
        sections.add("-?, -h, --help" to "Show help")
        sections.add("-v, --version" to "Show version")
        sections.add("" to "")
        sections.add("Commands:" to "")
        Tools.values().forEach { sections.add(it.command to it.description) }

        return sections.joinToString(separator = "\n") { (cmd, help) ->
            val fill = StringBuilder()
            (0..(25 - cmd.length)).forEach { _ -> fill.append(' ') }
            "$cmd$fill  $help"
        }
    }

    enum class Tools(
        val command: String, val description: String,
        private val f: (Array<String>) -> Unit
    ) {
        LOCI_ENRICHMENT_IN_REGIONS(
            "enrichmentInRegions",
            "Loci of interest enrichment in region sets compared to similar simulated loci.", { args ->
                EnrichmentInRegions.main(args)
            }),

        OVERLAP_REGIONS_WITH_EACH_LOCI(
            "overlapPerLocus",
            "For each locus from loci and of regions calculate: overlap(regions, i-th locus location)", { args ->
                OverlapRegionsWithEachLoci.main(args)
            });

        operator fun invoke(args: Array<String>) = f(args)
    }

    /**
     * Configure logging to console / file. Returns path to log file.
     */
        internal fun configureLogging(quiet: Boolean, debug: Boolean) {
        if (quiet) {
            Logs.quiet()
        }
        // Anyway we configure logging to file.
        Logs.addConsoleAppender(if (debug) Level.DEBUG else Level.INFO)
    }
}