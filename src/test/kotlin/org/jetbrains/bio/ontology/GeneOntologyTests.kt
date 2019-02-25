@file:Suppress("unused_variable")

import com.google.common.collect.HashMultimap
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.ontology.OboFile
import org.jetbrains.bio.ontology.Ontology
import org.jetbrains.bio.ontology.mapDepth
import org.jetbrains.bio.ontology.prune
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withResource
import org.jgrapht.EdgeFactory
import org.jgrapht.Graphs
import org.jgrapht.experimental.dag.DirectedAcyclicGraph
import org.junit.Test
import org.junit.runner.RunWith
import org.junit.runners.Parameterized
import org.junit.runners.Parameterized.Parameters
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class OboFileTest {
    @Test fun testParse() {
        withResource(OboFileTest::class.java, "go-basic-excerpt.obo") { path ->
            path.bufferedReader().use { reader ->
                val it = OboFile(Ontology.BIOLOGICAL_PROCESS, reader).iterator()
                val term1 = it.next()
                assertEquals("GO:0000001", term1.id)
                assertEquals("mitochondrion inheritance", term1.description)
                assertEquals(setOf("GO:0048308", "GO:0048311"), term1.children)

                val term2 = it.next()
                assertEquals("GO:0000002", term2.id)

                val term3 = it.next()
                assertEquals("GO:0000003", term3.id)
                assertFalse(term3.isObsolete)

                assertFalse(it.hasNext())
            }
        }
    }

    @Test fun testRead() {
        withResource(OboFileTest::class.java, "go-basic-excerpt.obo") { path ->
            val (g, _terms) = OboFile.read(Ontology.BIOLOGICAL_PROCESS, path)
            assertEquals(emptySet(),
                         setOf("GO:0000001", "GO:0000002", "GO:0000003") - g.vertexSet())
            assertFalse("GO:0000005" in g.vertexSet())
        }
    }
}

class MapDepthTest {
    @Test fun example() {
        val graph = listOf("1" to "2", "1" to "3",
                           "2" to "4", "2" to "3").toGraph()
        val depth = graph.mapDepth()

        assertEquals(0, depth["1"])
        assertEquals(1, depth["2"])
        assertEquals(1, depth["3"])
        assertEquals(2, depth["4"])
    }
}

class PruneTest {
    @Test fun example() {
        val graph = listOf("1" to "2", "1" to "3").toGraph()
        val associations = HashMultimap.create<String, String>().apply {
            putAll("1", "abc".map { it.toString() })
            putAll("2", listOf("b"))
            putAll("3", listOf("c"))
        }

        val pruned = associations.prune(graph)
        assertEquals(setOf("a"), pruned["1"])
        assertEquals(setOf("b"), pruned["2"])
        assertEquals(setOf("c"), pruned["3"])
    }
}

// The examples are from Figure 1 of
// Cao & Zhang "A Bayesian extension of the hypergeometric test for
//              functional enrichment analysis"
class PruneCaoZhangTest {
    @Test fun partial() {
        val graph = listOf("4" to "5", "4" to "6",
                           "5" to "7", "6" to "7",
                           "7" to "8").toGraph()

        val associations = HashMultimap.create<String, String>().apply {
            putAll("4", "abcdstuv".map { it.toString() })
            putAll("5", "abcdstu".map { it.toString() })
            putAll("6", "abdv".map { it.toString() })
            putAll("7", "abd".map { it.toString() })
            putAll("8", listOf("a"))
        }

        val pruned = associations.prune(graph)
        assertEquals(0, pruned["4"]!!.size)
        assertEquals(4, pruned["5"]!!.size)
        assertEquals(1, pruned["6"]!!.size)
        assertEquals(2, pruned["7"]!!.size)
        assertEquals(1, pruned["8"]!!.size)
    }

    @Test fun full() {
        val graph = listOf(
                "1" to "2", "1" to "3", "1" to "9",
                "2" to "4",
                "3" to "4", "3" to "6",
                "4" to "5",
                "5" to "7",
                "6" to "7",
                "7" to "8",
                "9" to "10", "9" to "11", "9" to "12",
                "10" to "13", "10" to "14",
                "11" to "15",
                "12" to "15", "12" to "17",
                "13" to "16",
                "14" to "17",
                "16" to "18",
                "17" to "19").toGraph()

        val associations = HashMultimap.create<String, String>().apply {
            putAll("1", "abcdefghijklmnopqrstuvwxyz".map { it.toString() })
            putAll("2", "abcdefstuvw".map { it.toString() })
            putAll("3", "abcdstuvxyz".map { it.toString() })
            putAll("4", "abcdstuv".map { it.toString() })
            putAll("5", "abcdstu".map { it.toString() })
            putAll("6", "abdtxy".map { it.toString() })
            putAll("7", "abdt".map { it.toString() })
            putAll("8", "a".map { it.toString() })
            putAll("9", "ghijklmnopqr".map { it.toString() })
            putAll("10", "ghjklopqr".map { it.toString() })
            putAll("11", "glmnopq".map { it.toString() })
            putAll("12", "hjlopqr".map { it.toString() })
            putAll("13", "gjklop".map { it.toString() })
            putAll("14", "hjklopq".map { it.toString() })
            putAll("15", "opq".map { it.toString() })
            putAll("16", "gjk".map { it.toString() })
            putAll("17", "hjlopq".map { it.toString() })
            putAll("18", listOf("g"))
            putAll("19", listOf("h"))
        }

        val pruned = associations.prune(graph)
        assertEquals(0, pruned["1"]!!.size)
        assertEquals(3, pruned["2"]!!.size)
        assertEquals(1, pruned["3"]!!.size)
        assertEquals(1, pruned["4"]!!.size)
        assertEquals(3, pruned["5"]!!.size)
        assertEquals(2, pruned["6"]!!.size)
        assertEquals(3, pruned["7"]!!.size)
        assertEquals(1, pruned["8"]!!.size)
        assertEquals(1, pruned["9"]!!.size)
        assertEquals(1, pruned["10"]!!.size)
        assertEquals(4, pruned["11"]!!.size)
        assertEquals(1, pruned["12"]!!.size)
        assertEquals(3, pruned["13"]!!.size)
        assertEquals(1, pruned["14"]!!.size)
        assertEquals(3, pruned["15"]!!.size)
        assertEquals(2, pruned["16"]!!.size)
        assertEquals(5, pruned["17"]!!.size)
        assertEquals(1, pruned["18"]!!.size)
        assertEquals(1, pruned["19"]!!.size)
    }
}

@RunWith(Parameterized::class)
class PruneRealDataTest(val ontology: Ontology) {
    @Test fun complete() {
        //XXX: this test forces mm9 annotations downloading if not available in genomes dir

        // TODO GO test data for test organism #1283
        val associations = ontology.associations(Genome["mm9"])
        val pruned = associations.prune(ontology.graph)
        // Not strictly equals, because some terms might loose all
        // associations as a result of pruning.
        assertTrue(associations.keySet().containsAll(pruned.keySet()))
        assertEquals(associations.values().toSet(), pruned.values().toSet())
    }

    @Test fun nonIncreasing() {
        val associations = ontology.associations(Genome["mm9"])
        val pruned = associations.prune(ontology.graph)
        for (term in associations.keySet()) {
            assertTrue(pruned[term].size <= associations[term].size,
                       message = term)
        }
    }

    companion object {
        @Parameters(name = "{0}}")
        @JvmStatic fun `data`() = listOf(Ontology.BIOLOGICAL_PROCESS)
    }
}

private fun List<Pair<String, String>>.toGraph(): DirectedAcyclicGraph<String, Pair<String, String>> {
    val g = DirectedAcyclicGraph(EdgeFactory { u: String, v -> u to v })
    for ((u, v) in this) {
        Graphs.addEdgeWithVertices(g, u, v)
    }

    return g
}
