package org.jetbrains.bio.genome

import org.jetbrains.bio.util.withResource
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * @author Roman.Chernyatchik
 */
class AnnotationsConfigLoaderTest {
    @Test
    fun checkVersionInYamlAgainstCode() {
        // XXX It's important to check that versions match, otherwise
        // AnnotationsConfig fails during static initialization and NoClassDefFound error is shown in logs
        withResource(AnnotationsConfigLoader::class.java, "annotations.yaml") { path ->
            val (version, _, _) = AnnotationsConfigLoader.parseYaml(path, AnnotationsConfigLoader.VERSION)
            assertEquals(AnnotationsConfigLoader.VERSION, version)
        }
    }

    @Test
    fun parseYamlOutdated() {
        withResource(AnnotationsConfigLoader::class.java, "test_annotations.yaml") { path ->
            val (version, mapping, yaml) = AnnotationsConfigLoader.parseYaml(path, 0)
            assertEquals(1, version)
            assertNull(mapping)
            assertEquals(6, yaml.genomes.size)
        }
    }

    @Test
    fun parseYaml() {
        withResource(AnnotationsConfigLoader::class.java, "test_annotations.yaml") { path ->
            val (version, mapping, _) = AnnotationsConfigLoader.parseYaml(path, 1)
            assertEquals(1, version)
            assertEquals(6, mapping!!.entries.size)
            assertEquals(listOf("ce6", "dm3", "hg19", "hg38", "hs1", "mm9"), mapping.keys.sorted())
            assertEquals("Drosophila melanogaster", mapping["dm3"]!!.species)
            assertEquals(listOf("mm9", "NCBIM37"), mapping["mm9"]!!.names)
            assertEquals("mm9", mapping["mm9"]!!.ucscAlias)
            assertEquals(
                Biomart(
                    "hsapiens_gene_ensembl",
                    "http://mar2016.archive.ensembl.org/biomart/martservice"
                ),
                mapping["hg38"]!!.mart
            )
            assertEquals(
                Biomart(
                    "dmelanogaster_gene_ensembl",
                    "http://dec2014.archive.ensembl.org/biomart/martservice"
                ),
                mapping["dm3"]!!.mart
            )
            assertEquals(
                "ftp://anonymous@ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz",
                mapping["mm9"]!!.gtfUrl
            )
            assertEquals("chrM", mapping["mm9"]!!.chrAltName2CanonicalMapping["MT"])
            assertEquals("chrM", mapping["dm3"]!!.chrAltName2CanonicalMapping["mitochondrion_genome"])
            assertTrue(mapping["mm9"]!!.ucscAnnLegacyFormat)
            assertFalse(mapping["hg19"]!!.ucscAnnLegacyFormat)

            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
                mapping["hg38"]!!.repeatsUrl
            )
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
                mapping["hg38"]!!.cytobandsUrl
            )
            assertNull(mapping["ce6"]!!.cytobandsUrl)
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz",
                mapping["hg38"]!!.gapsUrl
            )
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz",
                mapping["hg38"]!!.centromeresUrl
            )
            assertEquals(
                "http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v2.0/download/chm13v2.0_cytobands_allchrs.bed.gz",
                mapping["hs1"]!!.cytobandsUrl
            )
            assertNull(
                mapping["hs1"]!!.gapsUrl
            )
            assertEquals(
                "http://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz",
                mapping["hs1"]!!.repeatsUrl
            )
            assertEquals(
                "http://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/cpgIslandExt.bb",
                mapping["hs1"]!!.cpgIslandsUrl
            )
            assertNull(mapping["hg19"]!!.centromeresUrl)
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
                mapping["hg19"]!!.cpgIslandsUrl
            )
            assertNull(mapping["dm3"]!!.cpgIslandsUrl)
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.2bit",
                mapping["dm3"]!!.sequenceUrl
            )
            assertEquals(
                "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes",
                mapping["dm3"]!!.chromsizesUrl
            )
        }
    }

    @Test
    fun serialization() {
        val mapping =
            mapOf(
                "foo" to GenomeAnnotationsConfig(
                    species = "species",
                    ucscAlias = null,
                    names = listOf("foo"),
                    description = "description",
                    gtfUrl = "https://gtf",
                    chrAltName2CanonicalMapping = mapOf("MT" to "chrM"),
                    ucscAnnLegacyFormat = false,
                    sequenceUrl = "https://sequence",
                    chromsizesUrl = "https://chromsizes",
                    repeatsUrl = null,
                    cytobandsUrl = "https://cytobands",
                    gapsUrl = "https://gaps",
                    centromeresUrl = null,
                    cpgIslandsUrl = null,
                    mart = null
                )
            )
        withTempFile("foo", ".yaml") {
            AnnotationsConfigLoader.saveYaml(it, mapping)
            val (_, mapping2, _) = AnnotationsConfigLoader.parseYaml(it, 7)
            assertEquals(mapping, mapping2)
        }
    }

    @Test
    fun serializationPredefined() {
        withResource(AnnotationsConfigLoader::class.java, "test_annotations.yaml") { path ->
            val (_, mapping, _) = AnnotationsConfigLoader.parseYaml(path, 1)
            withTempFile("foo", ".yaml") {
                AnnotationsConfigLoader.saveYaml(it, mapping!!)
                val (_, mapping2, _) = AnnotationsConfigLoader.parseYaml(it, 7)
                assertEquals(mapping, mapping2)
            }
        }
    }

}