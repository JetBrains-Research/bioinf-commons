package org.jetbrains.bio.genome.query

import com.google.common.base.Joiner
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.format.BisulfiteBamParser
import org.jetbrains.bio.genome.methylome.Methylome
import org.jetbrains.bio.util.*
import java.io.IOException
import java.nio.file.Path

/**
 * @author Roman.Chernyatchik
 */
abstract class MethylomeQuery protected constructor(
    val genomeQuery: GenomeQuery,
    protected val dataSetId: String,
    val cellId: String,
    private vararg val properties: String
) :
    CachingInputQuery<Methylome>() {

    /** Reads a methylome using an **unrestricted** i.e. full [GenomeQuery]. */
    @Throws(IOException::class)
    protected abstract fun read(genomeQuery: GenomeQuery): Methylome

    companion object {
        /**
         * Creates a query for a given file path.
         *
         * At the moment only HDF5 and BAM files are supported.
         * @minBasePhred Min Phred quality threshold: '0' not filtration, '20' eRRBS guidelines recommendation
         */
        fun forFile(
            genomeQuery: GenomeQuery,
            cellId: String,
            path: Path,
            dataSetId: String = "File",
            verboseDescription: Boolean = false,
            minBasePhred: Byte = 20
        ): MethylomeQuery = FileMethylomeQuery(genomeQuery, cellId, path, dataSetId, verboseDescription, minBasePhred)
    }

    override fun getUncached(): Methylome {
        val binaryPath = Configuration.cachePath / "methylome" / "$id.npz"
        binaryPath.checkOrRecalculate("Methylome") { output ->
            read(genomeQuery).save(output.path)
        }

        return Methylome.lazy(genomeQuery, binaryPath)
    }

    override val id: String
        get() = "${dataSetId}_${genomeQuery.build}_$cellId" + Joiner.on('_').join(properties).let {
            if (it.isNotBlank()) "_$it" else ""
        }

    override val description: String
        get() {
            val propDesc = if (properties.isEmpty()) "" else "(${Joiner.on('_').join(properties)})"
            return "Methylome, dataset $dataSetId$propDesc for $cellId cells line genome ${genomeQuery.description}"
        }
}

private class FileMethylomeQuery(
    genomeQuery: GenomeQuery, cellId: String,
    private val path: Path,
    dataSetId: String,
    private val verboseDescription: Boolean,
    private val minBasePhred: Byte
) :
    MethylomeQuery(genomeQuery, dataSetId, cellId, path.stem) {

    override fun getUncached() = when (path.extension) {
        "npz" -> Methylome.lazy(genomeQuery, path)
        else -> super.getUncached()
    }


    override fun read(genomeQuery: GenomeQuery) = when (path.extension) {
        "npz" -> Methylome.lazy(genomeQuery, path)
        "bam" -> BisulfiteBamParser.parse(path, genomeQuery, minBasePhred)
        else -> error("unsupported extension: ${path.extension}")
    }

    // We use only phred filter if specified and the file name, because in practice filename already
    // contains cell ID and genome build.
    override val id: String get() = path.stem + if (minBasePhred == 0.toByte()) "" else ".$minBasePhred" + path.sha

    override val description: String
        get() {
            return if (verboseDescription) {
                path.name
            } else {
                "Methylome for ${genomeQuery.description}: $cellId, dataset $dataSetId(${path.name})"
            }
        }
}