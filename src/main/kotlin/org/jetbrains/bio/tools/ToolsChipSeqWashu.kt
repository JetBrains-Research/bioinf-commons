package org.jetbrains.bio.tools

import org.apache.log4j.Logger
import org.jetbrains.bio.dataset.ChipSeqTarget
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.util.*
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.util.*


/**
 * Main entry point to the toolsWashu repository scripts
 * See https://github.com/JetBrains-Research/toolsWashu
 */
@Suppress("Unused")
class ToolsChipSeqWashu(private val path: Path = DEFAULT_PATH) {

    init {
        check(path.isAccessible() && path.isDirectory) {
            "FAILED to find repository at: $path\n" +
                    "Repository: https://github.com/JetBrains-Research/toolsWashu"
        }
    }

    companion object {
        val DEFAULT_PATH = "/mnt/stripe/toolsWashu".toPath()
        internal val LOG = Logger.getLogger(ToolsChipSeqWashu::class.java)
    }

    /**
     * @return absolute path for given [relativePath]
     */
    private fun script(relativePath: String): String {
        val scriptPath = "$path/$relativePath".toPath()
        check(scriptPath.exists && scriptPath.isReadable) {
            "Missing file $scriptPath in ${this.path}"
        }
        return scriptPath.toAbsolutePath().toString()
    }

    fun runMACS2(
            genome: Genome,
            output: Path,
            fdr: Double,
            broadPeaks: Boolean = false,
            extParams: List<String> = listOf()
    ) {
        LOG.info("Working directory: $output")

        val suffix = (if (broadPeaks) {
            "broad"
        } else {
            "fdr"
        }) + "$fdr"

        val additionalMacs2Args = arrayListOf("-B")
        if (broadPeaks) {
            additionalMacs2Args.add("--broad")
            additionalMacs2Args.add("--broad-cutoff")
            additionalMacs2Args.add(fdr.toString())
        } else {
            additionalMacs2Args.add("-q")
            additionalMacs2Args.add(fdr.toString())
        }
        additionalMacs2Args.addAll(extParams)

        val args = arrayListOf(
                genome.build,
                "-", // Don't pass CHROM_SIZES argument to prevent big file creation
                suffix,
                "'${additionalMacs2Args.joinToString(separator = " ")}'",
                output.toString())

        val scriptPath = output / "run_MACS2_$suffix.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("parallel/macs2.sh")} ${args.joinToString(separator = " ")}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }


    fun runSicer(
            genome: Genome,
            output: Path,
            fdr: Double,
            extParams: List<String> = listOf()
    ) {
        LOG.info("Working directory: $output")

        val args = mutableListOf(output.toString(),
                genome.build,
                genome.chromSizesPath.toString(),
                fdr.toString())

        args += extParams

        val scriptPath = output / "run_Sicer.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("parallel/sicer.sh")} ${args.joinToString(separator = " ")}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runDiffReps(
            genome: Genome,
            readsPath: Path,
            output: Path,
            name: String,
            broadPeaks: Boolean
    ) {
        val chromSize = genome.chromSizesPath.toString()
        val scriptPath = output / "run_diffreps.sh"

        val args = if (broadPeaks) "--mode block --nsd broad" else "-me nb"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_diffreps.sh")} $name $chromSize $readsPath \"$args\"")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runDiffBind(
            output: Path,
            name: String,
            csvConfig: Path
    ) {
        val scriptPath = output / "run_diffbind_diff.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_diffbind.sh")} $name $csvConfig")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }


    fun runIntersection(output: Path, vararg beds: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$path && " +
                "bash ${script("bed/intersect.sh")} ${beds.map { it.toAbsolutePath() }.joinToString(" ")} > $output") {
            directory(output.parent.toFile())
        }
    }

    fun runUnion(output: Path, vararg beds: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$path && " +
                "bash ${script("bed/union.sh")} ${beds.map { it.toAbsolutePath() }.joinToString(" ")} > $output") {
            directory(output.parent.toFile())
        }
    }

    fun runReads2BW(readsPath: Path, chromSizesPath: Path, output: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$path && " +
                "bash ${script("scripts/reads2bw.sh")} $readsPath $chromSizesPath $output") {
            directory(output.parent.toFile())
        }
    }

    fun runReads2TagsBW(readsPath: Path, fragmentSize: Int, chromSizesPath: Path, output: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$path && " +
                "bash ${script("downstream/signals/bam2tagsbw.sh")} $readsPath $fragmentSize $chromSizesPath $output") {
            directory(output.parent.toFile())
        }
    }

    fun runPoolReads(paths: List<Path>, chromSizesPath: Path, outputBam: Path) {
        val scriptPath = outputBam.parent / "pool_reads.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("BAMS=()")
            writer.newLine()
            paths.forEach {
                writer.write("BAM=$(bash ${script("scripts/reads2bam.sh")} $it $chromSizesPath); BAMS+=(\"\$BAM\")")
                writer.newLine()
            }
            writer.write("samtools merge $outputBam \${BAMS[@]}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(outputBam.toFile())
        }
    }

    fun runDiffMacsPooled(name: String, output: Path) {
        val scriptPath = output / "run_macs_diff_pooled.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_macs_diff_pooled.sh")} $name $output")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }

    }

    fun runChipDiff(name: String, genome: Genome, output: Path, macsPooled: Path) {
        val scriptPath = output / "run_chipdiff_pooled.sh"

        val chromSize = genome.chromSizesPath.toString()

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_chipdiff.sh")} $name $output $chromSize $macsPooled")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runMacsBdgdiff(name: String, output: Path, macsPooled: Path) {
        val scriptPath = output / "run_macs_bdgdiff.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_macs_bdgdiff.sh")} $name $output $macsPooled")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runChipDiffPostProcess() {
        "python".exec("$path/downstream/diff/chipseq_diff_postprocess.py")
    }

    fun runRip(bamPath: Path, peakPath: Path) {
        val scriptPath = peakPath.parent / "run_rip.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$path")
            writer.newLine()
            writer.write("bash ${script("scripts/rip.sh")} $bamPath $peakPath")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(peakPath.parent.toFile())
        }
    }

    /**
     * Nucleotide consensus with regions supported by at least a half of all files.
     */
    fun medianNucleotideConsensus(resultFile: Path, tracks: List<Path>) {
        "bash".exec("-c",
                "bash $path/bed/consensus.sh -p 50 ${tracks.joinToString(separator = " ")} > $resultFile")
    }

    /**
     * Nucleotide consensus with regions supported by at least two files.
     */
    fun weakNucleotideConsensus(resultFile: Path, tracks: List<Path>) {
        "bash".exec("-c",
                "bash $path/bed/consensus.sh -c 2 ${tracks.joinToString(separator = " ")} > $resultFile")
    }
}
