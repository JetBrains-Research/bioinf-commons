package org.jetbrains.bio.genome

import com.google.common.cache.CacheBuilder
import com.google.common.collect.ListMultimap
import com.google.common.collect.Maps
import com.google.common.collect.Multimaps
import com.google.gson.GsonBuilder
import com.google.gson.reflect.TypeToken
import org.jetbrains.bio.genome.containers.minus
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.nio.file.Path
import java.util.*
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

/**
 * Useful gene groups (aka classes).
 *
 * @author Sergei Lebedev
 */
enum class GeneClass(
    val id: String, val description: String,
    private val predicate: (Transcript) -> Boolean
) {
    CODING("coding", "coding genes", { it.isCoding }),
    NON_CODING("non-coding", "non-coding genes", { !it.isCoding }),
    ALL("all", "all genes", { true });

    operator fun contains(transcript: Transcript) = predicate(transcript)
}

/**
 * Gene identifiers available in UCSC Genome Browser.
 *
 * @author Sergei Lebedev
 */
enum class GeneAliasType(val description: String) {
    GENE_SYMBOL("Gene Symbol"),
    ENSEMBL_ID("Ensembl Transcript ID"),   // unique.
    ENSEMBL_GENE_ID("Ensembl Gene ID");

    // TODO: load gene Synonyms and prev names?
    // hgnc <- read.delim(url("
    // http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&hgnc_dbtag=on&col=gd_app_sym&col=gd_aliases&status=Approved&status=Entry+Withdrawn&status_opt=3&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit&.cgifields=&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag"))
    //  https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit  //  //
    // use tsv data from
    // https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit//  //
    // see https://www.genenames.org/download/custom/
    // biomart:
    // https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__approved_symbol_1010%2Chgnc_gene__approved_name_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104%2Chgnc_gene__hgnc_alias_symbol__alias_symbol_108%2Chgnc_gene__hgnc_previous_symbol__previous_symbol_1012%2Chgnc_gene__hgnc_alias_name__alias_name_107%2Chgnc_gene__hgnc_previous_name__previous_name_1011ranscriptranscriptranscript
    // https://biomart.genenames.org/martform
}

/**
 * A zero-overhead [EnumMap] specialized to [GeneAliasType].
 *
 * Note that no container
 *
 * @author Sergei Lebedev
 */
class GeneAliasMap internal constructor(private val transcript: Transcript) : Map<GeneAliasType, String> {
    override val keys: Set<GeneAliasType> get() = setOf(*KEY_UNIVERSE)

    override val values: Collection<String> get() = KEY_UNIVERSE.map { get(it) }

    override val entries: Set<Map.Entry<GeneAliasType, String>>
        get() {
            return KEY_UNIVERSE.mapTo(LinkedHashSet()) { Maps.immutableEntry(it, get(it)) }
        }

    override fun get(key: GeneAliasType) = when (key) {
        GeneAliasType.GENE_SYMBOL -> transcript.geneSymbol
        GeneAliasType.ENSEMBL_ID -> transcript.ensemblId
        GeneAliasType.ENSEMBL_GENE_ID -> transcript.ensemblGeneId
    }.uppercase()

    override fun containsKey(key: GeneAliasType) = true

    override fun containsValue(value: String) = value in values

    override fun isEmpty() = false

    override val size: Int get() = KEY_UNIVERSE.size

    companion object {
        private val KEY_UNIVERSE = GeneAliasType.values()
    }
}

/** This is actually a transcript, not a gene. */
class Transcript(
    /** Ensembl transcript ID. */
    val ensemblId: String,
    /** Ensembl gene ID. */
    val ensemblGeneId: String,
    /** Gene symbol associated with this transcript. */
    val geneSymbol: String,
    /** Transcript location. */
    override val location: Location,
    /** Coding sequence range or `null` for non-coding transcripts. */
    val cdsRange: Range?,
    private val utr3End5: Int, // Could be out of transcript bound if not defined
    /** A list of sorted exon ranges. Empty for non-coding transcripts. */
    private val exonRanges: List<Range>
) : LocationAware {

    init {
        require(geneSymbol.isNotEmpty()) { "missing gene symbol" }
        require(utr3End5 != 0) { "UTR3 5' end expected to be positive or -1" }
        require(cdsRange == null || cdsRange.isNotEmpty()) {
            "CDS range expected to be null in non-coding and not empty in coding genes"
        }
    }

    val names: Map<GeneAliasType, String> get() = GeneAliasMap(this)

    val chromosome: Chromosome get() = location.chromosome
    val strand: Strand get() = location.strand

    val isCoding: Boolean get() = cdsRange != null

    /**
     * In ensembl genes gtf - same behaviour, CDS is just exons parts
     */
    val cds: List<Location>
        get() = when {
            !isCoding -> emptyList()
            else -> collectExonsIntersecting(cdsRange!!).map { it.on(chromosome, strand) }
        }

    private fun collectExonsIntersecting(bounds: Range) = exonRanges.mapNotNull { r ->
        if (r.startOffset >= bounds.startOffset && r.endOffset <= bounds.endOffset) {
            // most exons included in CDS bounds, let's do not generate extra
            // GC:
            r
        } else {
            val part = r.intersection(bounds)
            if (part == Range.EMPTY) null else part
        }
    }

    val exons: List<Location> get() = exonRanges.map { it.on(chromosome, strand) }

    val introns: List<Location> get() = (location.toRange() - exonRanges).map { it.on(chromosome, strand) }

    /**
     * Everything after the TSS but before the CDS.
     */
    val utr5: List<Location>
        get() = when {
            !isCoding -> emptyList()
            else -> {
                val utr5Bounds = if (strand.isPlus()) {
                    Range(location.startOffset, cdsRange!!.startOffset)
                } else {
                    Range(cdsRange!!.endOffset, location.endOffset)
                }
                collectExonsIntersecting(utr5Bounds).map { it.on(chromosome, strand) }
            }
        }

    /**
     * Everything after the CDS but before the transcription end site.
     */
    val utr3: List<Location>
        get() = when {
            !isCoding || utr3End5 == -1 -> emptyList()
            else -> {
                val utr3Bounds = when {
                    strand.isPlus() -> Range(utr3End5, location.endOffset)
                    else -> Range(location.startOffset, utr3End5 + 1)
                }
                collectExonsIntersecting(utr3Bounds).map { it.on(chromosome, strand) }
            }
        }

    fun length() = location.length()

    override fun toString(): String {
        return "$geneSymbol $location [exons: ${exons.size}, introns: ${introns.size}]"
    }

    override fun equals(other: Any?) = when {
        this === other -> true
        other !is Transcript -> false
        else -> ensemblId == other.ensemblId
    }

    override fun hashCode() = ensemblId.hashCode()
}

/**
 * Simplified notion of gene as a collection of [Transcript]
 */
data class Gene(
    val ensemblGeneId: String,
    val geneSymbol: String,
    val transcripts: Collection<Transcript>
) {
    // Gene location can be considered as longest isoform (according to the talk at ECCB-2018) or union of isoforms.
    val location: Location

    init {

        val locations = transcripts.map(Transcript::location)
        val ranges = locations.map { it.toRange() }

        location = ranges.fold(ranges.first(), Range::union)
            .on(locations.first().chromosome, locations.first().strand)
    }
}

fun groupTranscripts(transcripts: Iterable<Transcript>): List<Gene> {
    val groups = Multimaps.index(transcripts) {
        it!!.ensemblGeneId
    }.asMap()

    return groups.map { entry ->
        val (geneId, geneTranscripts) = entry
        Gene(geneId, geneTranscripts.first().geneSymbol, geneTranscripts)
    }
}


object Transcripts {
    private val LOG = LoggerFactory.getLogger(Transcripts::class.java)

    private const val FORMAT_VERSION = 1

    private val TRANSCRIPTS_CACHE = cache<Transcript>()
    private val BOUND5_INDEX_CACHE = CacheBuilder.newBuilder()
        .softValues()
        .initialCapacity(1)
        .build<Pair<Genome, Boolean>, Map<Chromosome, Triple<List<Transcript>, IntArray, BinaryLut>>>()

    /** Visible only for [TestOrganismDataGenerator]. */
    val GSON = GsonBuilder()
        .registerTypeAdapter(Range::class.java, Range.ADAPTER)
        .registerTypeAdapter(Location::class.java, Location.ADAPTER)
        .create()

    fun bound5Index(
        genome: Genome,
        onlyCodingGenes: Boolean
    ) = BOUND5_INDEX_CACHE.get(genome to onlyCodingGenes) {
        val allTranscripts = all(genome)

        allTranscripts.keySet().associate { chr ->
            val transcripts = allTranscripts[chr].filter { !onlyCodingGenes || it.isCoding }

            val bounds5 = transcripts.map { it.location.get5Bound() }.toIntArray()
            chr to (Triple(transcripts, bounds5, BinaryLut.of(bounds5, 24)))
        }
    }

    internal fun all(genome: Genome): ListMultimap<Chromosome, Transcript> {
        return TRANSCRIPTS_CACHE.get(genome) {
            try {
                loadTranscripts(genome, false)
            } catch (e: Exception) {
                LOG.warn("Failed to load transcripts for ${genome.build}. Trying to clear caches. Error:", e)
                loadTranscripts(genome, true)
            }
        }
    }

    /**
     * Just container for 2 vars. By some reason Pair<A,B> is being deserialized into LinkedTreeMap instead of Pair.
     */
    data class JsonTranscriptome(val vers: Int, val transcripts: List<Transcript>)

    fun cachedTranscriptsJsonPath(genesGtfPath: Path): Path {
        // JSON cached genes are loaded in ~2 secs, instead of 14-23 secs parsing from gtf
        return genesGtfPath.withExtension("json.gz")
    }

    private fun loadTranscripts(genome: Genome, cleanCache: Boolean): ListMultimap<Chromosome, Transcript> {
        val gtfFile = genome.genesGtfPath(false)
        val cachedTranscripts = cachedTranscriptsJsonPath(gtfFile)

        if (cleanCache) {
            cachedTranscripts.deleteIfExists()
        }

        cachedTranscripts.checkOrRecalculate("Processing genes") { (path) ->
            // ensure gtf file exists
            genome.genesGtfPath(true)

            val transcripts = LOG.time(level = Level.INFO, message = "Parsing genes $gtfFile") {
                gtfFile.bufferedReader().use { reader ->
                    GtfReader(reader, genome).readTranscripts()
                        // sort by 5' offset then by ensemblId (to be deterministic),
                        // it is required for bound5Index() correct work,
                        // please don't change it to sort by start offset
                        .sortedWith(compareBy({ it.location.get5Bound() }, { it.ensemblId }))
                }
            }
            path.bufferedWriter().use { GSON.toJson(JsonTranscriptome(FORMAT_VERSION, transcripts), it) }
        }

        val transcripts = LOG.time(
            level = Level.INFO,
            message = "Loading genes $gtfFile"
        ) {
            loadFromJson(cachedTranscripts)
        }
        return Multimaps.index(transcripts) { it!!.chromosome }
    }

    fun loadFromJson(cachedTranscripts: Path): List<Transcript> {
        val (vers, transcripts) = cachedTranscripts.bufferedReader().use {
            GSON.fromJson<JsonTranscriptome>(it, object : TypeToken<JsonTranscriptome>() {}.type)
        }
        require(vers == FORMAT_VERSION)
        return transcripts
    }

    enum class AssociationStrategy {
        SINGLE, TWO, BASAL_PLUS_EXT, MULTIPLE
    }

    /**
     * Determines the transcripts that could be regulated by the location. By default,
     * the association strategy is GREAT's "single nearest gene", adapted for use with transcripts.
     * See corresponding internal methods for more information.
     */
    fun associatedTranscripts(
        location: Location,
        strategy: AssociationStrategy = AssociationStrategy.SINGLE,
        limit: Int = 1000000,
        codingOnly: Boolean = false
    ): List<Transcript> = when (strategy) {
        AssociationStrategy.SINGLE -> associatedTranscriptsSingle(
            location, limit, codingOnly
        )
        AssociationStrategy.TWO -> associatedTranscriptsTwo(
            location, limit, codingOnly
        )
        AssociationStrategy.BASAL_PLUS_EXT -> associatedTranscriptsPlus(
            location, codingOnly = codingOnly, distal = limit
        )
        AssociationStrategy.MULTIPLE -> associatedTranscriptsMultiple(
            location, codingOnly = codingOnly, limit = limit
        )
    }

    fun associatedTranscriptsSingle(
        location: Location,
        limit: Int = 1000000,
        codingOnly: Boolean = false
    ): List<Transcript> {
        // XXX GREAT: "single nearest gene" strategy: we use it because it simple, feel free to change it to
        // XXX "basal plus ext" or "two nearest genes"

        // Description from http://bejerano.stanford.edu/help/display/GREAT/Association+Rules :
        // "GREAT calculates statistics by associating genomic regions with nearby genes and applying the gene
        // annotations to the regions. Association is a two step process. First, every gene is assigned a
        // regulatory domain. Then, each genomic region is associated with all genes whose regulatory domain
        // it overlaps."

        // In reality, the GREAT web tool reduces the location to its midpoint, and that's what we use here as well.

        // "Single nearest gene" strategy:
        // Gene regulatory domain definition: Each gene is assigned a regulatory domain that extends
        // in both directions to the midpoint between the gene's TSS and the nearest gene's TSS
        // but no more than the maximum extension in one direction (1000 kb)

        val chr = location.chromosome
        val (transcripts, bounds5, bound5Lut) = bound5Index(
            chr.genome, onlyCodingGenes = codingOnly
        ).getValue(chr)

        val midpoint = (location.startOffset + location.endOffset) / 2

        val (low, high) = bound5Lut.nearestElemDist(bounds5, midpoint)

        if (low == -1) {
            // actually empty transcripts list
            check(transcripts.isEmpty())
            return emptyList()
        }
        if (abs(bounds5[low] - midpoint) >= limit) {
            return emptyList()
        }
        return (low..high).map { transcripts[it] }
    }

    fun associatedTranscriptsTwo(
        location: Location,
        limit: Int = 1000000,
        codingOnly: Boolean = false
    ): List<Transcript> {
        // XXX GREAT: "two nearest genes" strategy

        // Description from http://bejerano.stanford.edu/help/display/GREAT/Association+Rules :
        // "GREAT calculates statistics by associating genomic regions with nearby genes and applying the gene
        // annotations to the regions. Association is a two step process. First, every gene is assigned a
        // regulatory domain. Then, each genomic region is associated with all genes whose regulatory domain
        // it overlaps."

        // In reality, the GREAT web tool reduces the location to its midpoint, and that's what we use here as well.

        // "Two nearest gene" strategy:
        // Gene regulatory domain definition: Each gene is assigned a regulatory domain that extends
        // in both directions to the nearest gene's TSS but no more
        // than the maximum extension in one direction (1000 kb)

        val chr = location.chromosome
        val (transcripts, bounds5, bound5Lut) = bound5Index(
            chr.genome, onlyCodingGenes = codingOnly
        ).getValue(chr)

        val midpoint = (location.startOffset + location.endOffset) / 2

        val (low, high) = bound5Lut.nearestElemLR(bounds5, midpoint)

        if (low == -1) {
            // actually empty transcripts list
            check(transcripts.isEmpty())
            return emptyList()
        }
        return (low..high).filter { abs(bounds5[it] - midpoint) <= limit }.map { transcripts[it] }
    }

    fun associatedTranscriptsMultiple(
        location: Location,
        limit: Int = 1000000,
        closeGenesAreaThreshold: Int = 50000,
        codingOnly: Boolean = false
    ): List<Transcript> {
        val chr = location.chromosome
        val (transcripts, bounds5, bound5Lut) = bound5Index(
            chr.genome, onlyCodingGenes = codingOnly
        ).getValue(chr)

        val midpoint = (location.startOffset + location.endOffset) / 2
        val midpointLeft = maxOf(0, midpoint - closeGenesAreaThreshold)
        val midpointRight = minOf(chr.length, midpoint + closeGenesAreaThreshold)

        if (bounds5.isEmpty()) {
            check(transcripts.isEmpty())
            return emptyList()
        }
        val (lowL, _) = bound5Lut.nearestElemLR(bounds5, midpointLeft)
        val (_, highR) = bound5Lut.nearestElemLR(bounds5, midpointRight)
        val low = lowL
        val high = highR

        check(low != -1)
        check(high != -1)

        return (low..high).filter { abs(bounds5[it] - midpoint) <= limit }.map { transcripts[it] }
    }

    fun associatedTranscriptsPlus(
        location: Location,
        upstream: Int = 5000,
        downstream: Int = 1000,
        distal: Int = 1000000,
        codingOnly: Boolean = false
    ): List<Transcript> {
        // XXX GREAT: "basal plus extension" strategy

        // Description from http://bejerano.stanford.edu/help/display/GREAT/Association+Rules :
        // "GREAT calculates statistics by associating genomic regions with nearby genes and applying the gene
        // annotations to the regions. Association is a two step process. First, every gene is assigned a
        // regulatory domain. Then, each genomic region is associated with all genes whose regulatory domain
        // it overlaps."

        // In reality, the GREAT web tool reduces the location to its midpoint, and that's what we use here as well.

        // "Basal plus extension" strategy:
        // Each gene is assigned a basal regulatory domain of a minimum distance upstream and downstream of the TSS
        // (regardless of other nearby genes). The gene regulatory domain is extended in both directions
        // to the nearest gene's basal domain but no more than a maximum extension in one direction.
        // When extending the regulatory domain of gene G beyond its basal domain, the extension to the "left"
        // extends until it reaches the first basal domain of any gene whose transcription start site
        // is "left" of G's transcription start site (and analogously for extending "right").

        val chr = location.chromosome
        val (transcripts, bounds5, bound5Lut) = bound5Index(
            chr.genome, onlyCodingGenes = codingOnly
        ).getValue(chr)

        val midpoint = (location.startOffset + location.endOffset) / 2

        // check basal regulatory domain
        val (lowBas, highBas) = bound5Lut.elemWithinDist(bounds5, midpoint, max(upstream, downstream))
        if (lowBas != -1) {
            val candidates = (lowBas..highBas).filter {
                bounds5[it] - midpoint in transcripts[it].strand.choose(
                    -downstream..upstream, -upstream..downstream
                )
            }
                .map { transcripts[it] }
            if (candidates.isNotEmpty()) return candidates
        }
        // check extended regulatory domain -- this is very tricky; we do it mostly separately for the two strands
        val plus = bound5Lut.framingElem(bounds5, midpoint) { transcripts[it].strand == Strand.PLUS }
        val minus = bound5Lut.framingElem(bounds5, midpoint) { transcripts[it].strand == Strand.MINUS }
        val leftBasalBoundaryPlus = plus.map { bounds5[it] + downstream }
            .filter { it < midpoint }.maxOrNull()
        val rightBasalBoundaryPlus = plus.map { bounds5[it] - upstream }
            .filter { it > midpoint }.minOrNull()
        val leftBasalBoundaryMinus = minus.map { bounds5[it] + upstream }
            .filter { it < midpoint }.maxOrNull()
        val rightBasalBoundaryMinus = minus.map { bounds5[it] - downstream }
            .filter { it > midpoint }.minOrNull()
        val leftBasalBoundary = when {
            leftBasalBoundaryPlus == null -> leftBasalBoundaryMinus
            leftBasalBoundaryMinus == null -> leftBasalBoundaryPlus
            else -> max(leftBasalBoundaryPlus, leftBasalBoundaryMinus)
                .let { if (midpoint - it <= distal) it else null }
        }
        val rightBasalBoundary = when {
            rightBasalBoundaryPlus == null -> rightBasalBoundaryMinus
            rightBasalBoundaryMinus == null -> rightBasalBoundaryPlus
            else -> min(rightBasalBoundaryPlus, rightBasalBoundaryMinus)
                .let { if (it - midpoint <= distal) it else null }
        }
        val resultsPlus = plus.filter {
            bounds5[it] + downstream == leftBasalBoundary || bounds5[it] - upstream == rightBasalBoundary
        }
        val resultsMinus = minus.filter {
            bounds5[it] + upstream == leftBasalBoundary || bounds5[it] - downstream == rightBasalBoundary
        }
        return (resultsPlus + resultsMinus).sorted().map { transcripts[it] }
    }


    /**
     * Returns the distance implemented in "single nearest gene" strategy from GREAT
     * (see above), i.e. the distance between transcript's TSS and location's midpoint;
     * or -1 when any argument is null or the chromosomes differ
     */
    fun greatDistance(transcript: LocationAware?, location: Location?): Int {
        if (transcript == null || location == null || transcript.location.chromosome != location.chromosome) return -1
        val midpoint = (location.startOffset + location.endOffset) / 2
        return abs(transcript.location.get5Bound() - midpoint)
    }

    /**
     * Returns the signed distance implemented in "single nearest gene" strategy from GREAT
     * (see above), i.e. the distance between transcript's TSS and location's midpoint;
     * positive when the midpoint is downstream of the TSS and negative otherwise.
     */
    fun signedGreatDistance(transcript: LocationAware, location: Location): Int {
        val midpoint = (location.startOffset + location.endOffset) / 2
        return midpoint - transcript.location.get5Bound()
    }
}