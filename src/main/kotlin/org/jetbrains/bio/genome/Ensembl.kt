package org.jetbrains.bio.genome

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.log4j.Logger
import org.jetbrains.bio.util.bufferedWriter
import java.io.BufferedReader
import java.io.BufferedWriter
import java.nio.file.Path
import java.util.*
import kotlin.collections.ArrayList
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

object Ensembl {

    fun convertGTF(genome: Genome, inputStream: BufferedReader, outputPath: Path) {
        val mapping = genome.annotationsConfig?.chrAltName2CanonicalMapping ?: emptyMap()
        val writer = outputPath.bufferedWriter()
        for (line in inputStream.lines()) {
            if (line.startsWith("#")) {
                writer.write(line)
                writer.newLine()
            } else {
                val parts = line.split("\t").toMutableList()
                val name = parts[0]
                parts[0] = mapping[name] ?: when {
                    name.startsWith("chr") -> name
                    else -> "chr$name"
                }
                writer.write(parts.joinToString(separator = "\t"))
                writer.newLine()
            }

        }
        writer.close()
    }
}

typealias Attributes = Array<String?>
class GtfReader(val reader: BufferedReader, val genome: Genome) {

    fun readTranscripts(): List<Transcript> {
        var hasDetailedUTRInfo: Boolean? = null
        val transcriptsMap = HashMap<String, TranscriptInfo>()
        val genomeQuery = genome.toQuery()
        val chrNamesMap =  genome.chromosomeNamesMap
        for (line in reader.lineSequence()) {
            val featureType = parseLine(line, transcriptsMap, genomeQuery, chrNamesMap)

            if (hasDetailedUTRInfo == null) {
                if (featureType == "five_prime_utr" || featureType == "three_prime_utr") {
                    hasDetailedUTRInfo = true
                } else if (featureType == "UTR") {
                    hasDetailedUTRInfo = false
                }
            }
        }
        if (hasDetailedUTRInfo == null) {
            hasDetailedUTRInfo = false
        }

        val transcripts = ArrayList<Transcript>()
        for (tInfo in transcriptsMap.values) {
            //XXX: Transcript could hypothetically start with intron
            // "t" is optional field, not available in legacy gtf (e.g. for hg18) so if not specified
            // use just merged exons bounds. Otherwise use GTF annotation
            //
            // P.S: At the moment all existing GTF annotations doesn't contain transcripts which starts
            // with introns.
            val transcript = if (tInfo.transcript != null) {
                // assert(mergeLocations(tInfo.exons) == tInfo.t)
                tInfo.transcript!!
            } else {
                mergeLocations(tInfo.exons)
            }

            val sortedExonRanges = tInfo.exons.map { it.toRange() }.sorted()

            var utr3End5 = -1
            val cdsBounds: Location?
            if (tInfo.cds.isEmpty()) {
                // non-coding
                cdsBounds = null
            } else {
                cdsBounds = mergeLocations(tInfo.cds)

                //Validate that we parse 'start_codon' GTF correctly:
                if (tInfo.startCodon != -1) {
                    // 'start_codon' field my be absent, not clear why, may be there are not AUG at CDS 5' end
                    // But if it is specified, ensure that is the same as CDS start
                    val cdsStart = cdsBounds.get5Bound()
                    check(cdsStart == tInfo.startCodon) {
                        "[START_CODON] We assume 'start_codon' (${tInfo.startCodon}) starts in CDS 5'end ($cdsStart)" +
                                " but is false here. Gene: ${tInfo.transcriptId}"
                    }
                }

                //'stop_codon' and UTR3 story:
                //
                // In new GTF format we normally have 'three_prime_utr' field, in some previous versions only
                // 'UTR' both for UTR3, UTR5. And for hg18 no UTR info.
                // Mostly in GTF files UTR3' starts after 3-bp lenght stop_codon which could be split by exons
                // but some times 'stop_codon' is missed and UTR3 starts at CDS 3' end + 1 or even
                // 'stop_codon' my be included in CDS
                utr3End5 = if (hasDetailedUTRInfo) {
                    if (tInfo.utr3.isEmpty()) {
                        -1
                    } else {
                        mergeLocations(tInfo.utr3).get5Bound()
                    }
                } else {
                    // No 'three_prime_utr' => detect UTR3 5' end as next position after 3-bp stop codon.
                    // Stop codon starts after CDS 3'. All offsets are calculated in transcript's exome, not in genome.
                    Companion.determineUTR3End5(cdsBounds, sortedExonRanges, tInfo.transcriptId)
                }
            }

            transcripts.add(Transcript(tInfo.transcriptId, tInfo.geneId, tInfo.geneName,
                                       transcript, cdsBounds?.toRange(), utr3End5,
                                       sortedExonRanges))
        }
        return transcripts
    }

    private fun parseLine(
        line: String, transcriptsMap: HashMap<String, TranscriptInfo>,
        genomeQuery: GenomeQuery,
        chrsNamesMapping: Map<String, Chromosome>
    ): String? {

        if (line.startsWith("#")) {
            return null
        }

        val parts = line.split("\t")

        /*
         * Different genome variants has different sets of contigs.
         * 1. old: exon, CDS, start_codon, stop_codon
         * 2. 75: exon, CDS, start_codon, stop_codon, transcript, UTR
         * 3. 87: exon, CDS, stop_codon, gene, transcript, three_prime_utr, five_prime_utr
         */
        val chrStr = parts[0]
        val type = parts[2]

        val chr = chrsNamesMapping[chrStr]?.let { chr ->
            if (genomeQuery.accepts(chr)) chr else null
        } ?: return type
        
        when (type) {
            // just not to write long if condition:
            "transcript", "exon", "CDS", "start_codon", "three_prime_utr" -> { /* noop */ }
            else -> return type
        }

        val start = Integer.parseInt(parts[3])
        val end = Integer.parseInt(parts[4])
        val strand = parts[6].toStrand()
        val attributes = parseAttributes(parts[8])

        val transcriptId = attributes[GtfAttrTypes.TRANSCRIPT_ID.ordinal]!!

        val transcriptInfo = transcriptsMap.getOrPut(transcriptId) {
            TranscriptInfo(transcriptId,
                           attributes[GtfAttrTypes.GENE_NAME.ordinal]!!,
                           attributes[GtfAttrTypes.GENE_ID.ordinal]!!)
        }

        val location = Location(start - 1, end, chr, strand)

        when (type) {
            "exon" -> {  transcriptInfo.exons.add(location) }
            "CDS" -> { transcriptInfo.cds.add(location) }
            "transcript" -> { transcriptInfo.transcript = location }
            "three_prime_utr" -> { transcriptInfo.utr3.add(location) }
            "start_codon" -> {
                // Used only for data validation
                if (transcriptInfo.startCodon == -1) {
                    // no defined
                    transcriptInfo.startCodon = location.get5Bound()
                } else {
                    if (strand.isPlus()) {
                        transcriptInfo.startCodon = min(transcriptInfo.startCodon, location.get5Bound())
                    } else {
                        transcriptInfo.startCodon = max(transcriptInfo.startCodon, location.get5Bound())
                    }
                }
            }
        }
        return type
    }

    private fun mergeLocations(loci: List<Location>): Location {
        val start = (loci.map { it.startOffset }).min()!!
        val end = (loci.map { it.endOffset }).max()!!
        val chromosome = loci[0].chromosome
        val strand = loci[0].strand

        return Location(start, end, chromosome, strand)
    }

    private fun parseAttributes(rest: String): Attributes {
        val attrTypes = GtfAttrTypes.values()
        val attributes = Array<String?>(attrTypes.size) { null }

        for (chunk in rest.split(";")) {
            val trimmed = chunk.trimStart()

            for (i in 0 until attrTypes.size) {
                val key = attrTypes[i].key

                // if attr not set & is our key
                if (trimmed.startsWith(key)) {
                    val value = trimmed.substring(key.length).trimStart()
                    check(value[0] == '"' && value.last() == '"') { "Cannot parse: $key value, attrs list = $rest" }
                    attributes[i] = value.substring(1, value.length - 1)
                    break
                }
            }
        }

        if (attributes[GtfAttrTypes.GENE_NAME.ordinal] == null) {
            // not all gtf genes has 'gene_name' attr, i.e. defines gene symbol, e.g not info in ce11, rn5
            attributes[GtfAttrTypes.GENE_NAME.ordinal] = attributes[GtfAttrTypes.GENE_ID .ordinal]
        }

        check(attributes.all { it != null }) {
            val missed = attrTypes.zip(attributes).filter { it.second == null }.map { it.first.key }
            "Required attributes [${missed.joinToString()}] not all required attrs found in: $rest"
        }
        return attributes
    }

    enum class GtfAttrTypes(val key: String) {
        TRANSCRIPT_ID("transcript_id"),
        GENE_ID("gene_id"),
        GENE_NAME("gene_name");
    }

    class TranscriptInfo(val transcriptId: String,
                         val geneName: String,
                         val geneId: String) {
        val exons: MutableList<Location> = ArrayList()
        val cds: MutableList<Location> = ArrayList()
        var transcript: Location? = null
        var gene: Location? = null
        val utr3: MutableList<Location> = ArrayList(4)
        var startCodon: Int = -1
    }

    companion object {
        private val LOG = Logger.getLogger(GtfReader::class.java)
        fun determineUTR3End5(cdsBounds: Location, sortedExonRanges: List<Range>, transcriptId: String): Int {
            val strand = cdsBounds.strand
            val cdsEnd3 = cdsBounds.get3Bound(0)

            var cdsEnd3ExonIdx = -1
            val exonsNumber = sortedExonRanges.size
            for (i in 0 until exonsNumber) {
                val exon = sortedExonRanges[i]
                if (cdsEnd3 in exon) {
                    cdsEnd3ExonIdx = i
                    break
                }
            }
            check(cdsEnd3ExonIdx >= 0) { "Cannot find $cdsEnd3 in exons for $transcriptId" }
            val cdsEnd3Exon = sortedExonRanges[cdsEnd3ExonIdx]
            val cdsEnd3ExonEnd3 = Location.get3Bound(cdsEnd3Exon.startOffset, cdsEnd3Exon.endOffset, strand, 0)
            val cdsExonNonCodingLength = abs(cdsEnd3ExonEnd3 - cdsEnd3)
            if (cdsExonNonCodingLength > 3) {
                return cdsBounds.get3Bound(3 + 1)
            }

            var restCodonPart = 3 - cdsExonNonCodingLength
            if ((cdsExonNonCodingLength == 0 || restCodonPart == 0)
                    && strand.choose(cdsEnd3ExonIdx == exonsNumber - 1,
                                     cdsEnd3ExonIdx == 0)) {
                // no more exons
                return -1
            }

            // partly intersects:
            val progression = if (strand.isPlus()) {
                cdsEnd3ExonIdx + 1 until exonsNumber
            } else {
                cdsEnd3ExonIdx - 1 downTo 0
            }

            for (i in progression) {
                val exon = sortedExonRanges[i]
                val lastExon = strand.choose(i == exonsNumber - 1,i == 0)
                if (restCodonPart < exon.length()) {
                    return exon.let { (startOffset, endOffset) ->
                        Location.get5Bound(startOffset, endOffset, strand, restCodonPart)
                    }
                } else if (lastExon && restCodonPart == exon.length()) {
                    return -1
                }
                restCodonPart -= exon.length()
            }
            // partial exons
            LOG.warn("Cannot detect UTR3 5' end for $transcriptId. Not enough space for STOP codon only " +
                             "${3 - restCodonPart} bp available in exons after CDS 3' end.")
            // Unfortunately exons ENSMUST00000168714 in Mus_musculus.NCBIM37.67.gtf.gz have not enough space
            // for stop codon. So let's return -1 (no UTR3) or maybe first exonic bp after CDS 3' end.
            return -1
        }
    }
}

fun writeGtf(writer: BufferedWriter, transcripts: Collection<Transcript>) {
    val csvPrinter = CSVPrinter(writer, CSVFormat.TDF.withQuote(null))
    for (transcript in transcripts) {
        val baseAttributes = listOf(
                "gene_id" to transcript.ensemblGeneId,
                "transcript_id" to transcript.ensemblId,
                "gene_name" to transcript.geneSymbol)

        val location = transcript.location
        val feature = GtfFeature(
                location,
                "biomarkt",
                FeatureType.TRANSCRIPT,
                attributes = baseAttributes)
        feature.write(csvPrinter)

        var number = 1

        var frame = 0

        val exons = if (transcript.strand.isPlus()) transcript.exons else transcript.exons.reversed()
        for (exon in exons) {
            GtfFeature(
                    exon,
                    "biomarkt",
                    FeatureType.EXON,
                    attributes = baseAttributes + ("exon_number" to number.toString())
            ).write(csvPrinter)
            val cdsBounds = transcript.cdsRange
            if (cdsBounds != null && cdsBounds.intersects(exon.toRange())) {
                val cdsPart = cdsBounds.intersection(exon.toRange())
                GtfFeature(
                        cdsPart.on(location.chromosome, location.strand),
                        "biomarkt",
                        FeatureType.CDS,
                        frame = frame,
                        attributes = baseAttributes + ("exon_number" to number.toString())
                ).write(csvPrinter)

                frame = (frame + cdsPart.length()) % 3
            }
            number++
        }
    }
}

enum class FeatureType(val featureName: String) {
    GENE("gene"),
    TRANSCRIPT("transcript"),
    EXON("exon"),
    CDS("CDS")
}


data class GtfFeature(val location: Location,
                      val source: String,
                      val type: FeatureType,
                      val score: Double? = null,
                      val frame: Int? = null,
                      val attributes: List<Pair<String, String>>) {


    fun write(csvPrinter: CSVPrinter) {
        csvPrinter.printRecord(
                location.chromosome.name,
                source,
                type.featureName,
                location.startOffset + 1,
                location.endOffset,
                score?.toString() ?: ".",
                location.strand.char,
                frame?.toString() ?: ".",
                attributes.map { "${it.first} \"${it.second}\"" }.joinToString("; ")
        )
    }

}