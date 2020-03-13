@file:Suppress("unused_variable")

import org.jetbrains.bio.genome.OboFile
import org.jetbrains.bio.genome.Ontology
import org.jetbrains.bio.genome.mapDepth
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withResource
import org.jgrapht.EdgeFactory
import org.jgrapht.Graphs
import org.jgrapht.experimental.dag.DirectedAcyclicGraph
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse

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



private fun List<Pair<String, String>>.toGraph(): DirectedAcyclicGraph<String, Pair<String, String>> {
    val g = DirectedAcyclicGraph(EdgeFactory { u: String, v -> u to v })
    for ((u, v) in this) {
        Graphs.addEdgeWithVertices(g, u, v)
    }

    return g
}
