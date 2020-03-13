package org.jetbrains.bio.genome.data

import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.stem
import org.jetbrains.bio.util.withTempDirectory
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import java.io.StringWriter
import java.util.regex.Pattern
import kotlin.test.*

class DataConfigTest {
    @Test
    fun consistency() {
        withConfig { config ->
            assertTrue(DataConfig.Key("H3K4me3", Cell("test1")) in config.tracksMap)
            val h3k4me3 = config.tracks["H3K4me3"]
            val reps = arrayListOf<String>()
            for (e in (h3k4me3!!["test1"]!! as Map<String, String>).entries) {
                reps.add(e.key)
            }
            assertEquals("[rep1, rep2]", reps.toString())
            assertTrue(DataConfig.Key("H3K27ac", Cell("test0")) in config.tracksMap)
            assertTrue(DataConfig.Key("H3K9me3", Cell("test1")) in config.tracksMap)
            assertTrue(DataConfig.Key("H3K4me3", Cell("test2")) in config.tracksMap)
        }
    }

    @Test
    fun technicalReplicate() {
        withConfig { config ->
            config.tracksMap[DataConfig.Key("H3K9me3", Cell("test1"))]!!.forEach {
                assertTrue { Pattern.matches(DataConfig.TECHNICAL_REPLICATE_PATTERN, it.first) }
            }
        }
    }


    @Test
    fun metaConsistency() {
        withConfig { config ->
            val section = config.tracksMap[DataConfig.Key("H3K27ac", Cell("test0"))]
            assertEquals(setOf("rep1", "rep2"), section!!.map { it.first }.toSet())
            assertTrue { "meta" in section.first().second.meta }
            assertTrue { section.last().second.meta.isEmpty() }
        }

    }

    @Test
    fun auxConsistency() {
        withConfig { config ->
            val key = ReplicateDataKey("test5") { it == true }
            assertTrue { config["H3K27ac", key] }
            assertFalse { config["H3K36me3", key] }
        }
    }

    @Test
    fun labeledConsistency() {
        withConfig { config ->
            val section = config.tracksMap[DataConfig.Key("H3K4me3", Cell("test1"))]
            assertEquals(setOf("rep1", "rep2"), section!!.map { it.first }.toSet())
        }
    }

    @Test
    fun implicitConsistency() {
        withConfig { config ->
            val section = config.tracksMap[DataConfig.Key("H3K4me3", Cell("test2"))]
            assertEquals(1, section!!.size)
        }
    }

    @Test
    fun loadReloadDump() {
        withConfig { config ->
            val writer1 = StringWriter()
            config.genomeQuery
            config.tracksMap
            config.save(writer1)
            val save1 = writer1.toString()
            val reloaded = DataConfig.load(save1.reader(), "test")
            assertNotNull(reloaded)
            assertTrue('!' !in save1) // no YAML reads.
            assertTrue('"' !in save1) // no quotes
            val writer2 = StringWriter()
            config.save(writer2)
            val save2 = writer1.toString()
            assertEquals(save1, save2)
        }
    }

    @Test
    fun equalsHashCode() {
        withTempFile("test1", ".bed") { p ->
            val config1 = """genome: to1
tracks:
   H3K4me3:
      t:
      - $p
"""

            val config2 = """genome: to1
tracks:
   H3K4me3:
      t:
      - $p
"""
            val dc1 = DataConfig.load(config1.reader(), "test")
            val dc2 = DataConfig.load(config2.reader(), "test")
            assertEquals(dc1, dc2)
            assertEquals(dc1.hashCode(), dc2.hashCode())
        }
    }

    @Test
    fun testMethylation() {
        withTempFile("test1", ".bam") { p ->
            val config = """genome: to1
tracks:
   methylation:
      t:
      - $p
"""
            val dc = DataConfig.load(config.reader(), "test")
            assertEquals(1, dc.tracksMap.size)
            assertEquals(1, dc.tracksMap.keys.size)
            assertEquals(DataConfig.Key(METHYLATION.id, Cell("t")), dc.tracksMap.keys.first())
        }
    }

    @Test
    fun testNonExistentFile() {
        assertFailsWith<IllegalStateException> {
            val config = """genome: to1
tracks:
   methylation:
      t:
      - /foo/bar/baz
"""
            DataConfig.load(config.reader(), "test")
        }
    }


    @Test
    fun testTranscription() {
        withTempFile("test1", ".fastq") { p ->
            val config = """genome: to1
tracks:
   transcription:
      t:
      - $p
"""
            val dc = DataConfig.load(config.reader(), "test")
            assertEquals(1, dc.tracksMap.size)
            assertEquals(1, dc.tracksMap.keys.size)
            assertEquals(DataConfig.Key(TRANSCRIPTION.id, Cell("t")), dc.tracksMap.keys.first())
        }
    }


    @Test
    fun testSave() {
        withTempFile("test1", ".bam") { p ->
            val config = """genome: to1
tracks:
   methylation:
      t:
      - $p
"""
            val dc = DataConfig.load(config.reader(), "test")
            val writer = StringWriter()
            dc.save(writer)
            assert('#' !in writer.toString())
        }
    }

    @Test
    fun testId() {
        withTempFile("foo", ".yaml") { f ->
            withConfig { config ->
                assertEquals("test", config.id)
                f.bufferedWriter().use(config::save)
                val loaded = DataConfig.load(f)
                assertEquals(f.stem, loaded.id)
            }
        }
    }

    @Test
    fun testFailedTrack() {
        withConfig { config ->
            assertTrue(config.tracksMap[DataConfig.Key("H3K27ac", Cell("test0"))]!![0].second.failedTrack)
        }
    }


    companion object {
        private val METHYLATION = DataType("methylation")
        private val TRANSCRIPTION = DataType("transcription")

        fun withConfig(block: (DataConfig) -> Unit) {
            withTempFile("test1", ".bed") { p11 ->
                withTempFile("test1", ".bed") { p12 ->
                    withTempFile("test2", ".bed") { p2 ->
                        withTempFile("test3", ".bam") { p3 ->
                            withTempDirectory("test4") { p4 ->
                                withTempFile("test5", ".txt") { p5 ->
                                    val config = DataConfig.load("""genome: to1
tracks:
   H3K27ac:
      test0:
        rep1:
          path: $p11
          meta: true
          failed: true
        rep2: $p12
   H3K4me3:
      test1:
        rep1: $p11
        rep2: $p12
      test2:
      - $p2
   H3K9me3:
      test1:
      - $p2
   methylation:
      test1:
        rep1: $p3
   transcription:
      test1:
      - $p4
   mirna:
      test1:
        rep1: $p5
aux:
   H3K27ac:
      test5: true
""".reader(), "test")
                                    block(config)
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


class ReplicatedDataTest {
    @Test
    fun testList() {
        withTempFile("test1", ".bam") { p ->
            assertNotNull(ReplicateData.of(p.toString()))
        }
    }

    @Test
    fun testMeta() {
        withTempFile("test1", ".fastq") { p ->
            val data = ReplicateData.of(mapOf("path" to p.toString(), "failed" to true))
            assertNotNull(data)
            assertEquals(p, data.path)
            assertTrue(data.failedTrack)
        }
    }

    @Test
    fun testMetaOutlierNoTag() {
        withTempFile("test1", ".fastq") { p ->
            val data = ReplicateData.of(mapOf("path" to p.toString()))
            assertFalse(data.failedTrack)
        }
    }
}