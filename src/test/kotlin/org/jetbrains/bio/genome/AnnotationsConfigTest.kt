package org.jetbrains.bio.genome

import org.jetbrains.bio.io.BedParserTest
import org.jetbrains.bio.util.withResource
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * @author Roman.Chernyatchik
 */
class AnnotationsConfigTest {
    @Test
    fun checkVersionInYamlAgainstCode() {
        // XXX It's important to check that versions match, otherwise
        // AnnotationsConfig fails during static initialization and NoClassDefFound error is shown in logs
        withResource(AnnotationsConfig::class.java, "annotations.yaml") { path ->
            val (version, _, _) = AnnotationsConfig.parseYaml(path, AnnotationsConfig.VERSION)
            assertEquals(AnnotationsConfig.VERSION, version)
        }
    }

    @Test
    fun parseYamlOutdated() {
        withResource(BedParserTest::class.java, "test_annotations.yaml") { path ->
            val (version, mapping, yaml) = AnnotationsConfig.parseYaml(path, 0)
            assertEquals(1, version)
            assertNull(mapping)
            assertEquals(5, yaml.genomes.size)
        }
    }

    @Test
    fun parseYaml() {
        withResource(BedParserTest::class.java, "test_annotations.yaml") { path ->
            val (version, mapping, _) = AnnotationsConfig.parseYaml(path, 1)
            assertEquals(1, version)
            assertEquals(5, mapping!!.entries.size)
            assertEquals(listOf("ce6", "dm3", "hg19", "hg38", "mm9"), mapping.keys.sorted())
            assertEquals("Drosophila melanogaster", mapping["dm3"]!!.species)
            assertEquals("NCBIM37", mapping["mm9"]!!.alias)
            assertEquals(Mart(
                    "hsapiens_gene_ensembl",
                    "http://mar2016.archive.ensembl.org/biomart/martservice"),
                    mapping["hg38"]!!.mart)
            assertEquals(Mart(
                    "dmelanogaster_gene_ensembl",
                    "http://dec2014.archive.ensembl.org/biomart/martservice"),
                    mapping["dm3"]!!.mart)
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
                    mapping["hg38"]!!.repeatsUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
                    mapping["hg38"]!!.cytobandsUrl)
            assertNull(mapping["ce6"]!!.cytobandsUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz",
                    mapping["hg38"]!!.gapsUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz",
                    mapping["hg38"]!!.centromeresUrl)
            assertNull(mapping["hg19"]!!.centromeresUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
                    mapping["hg19"]!!.cpgIslandsUrl)
            assertNull(mapping["dm3"]!!.cpgIslandsUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.2bit",
                    mapping["dm3"]!!.sequenceUrl)
            assertEquals(
                    "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes",
                    mapping["dm3"]!!.chromsizesUrl)
        }
    }
}