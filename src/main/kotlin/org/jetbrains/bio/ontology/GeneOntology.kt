package org.jetbrains.bio.ontology

import com.google.common.annotations.VisibleForTesting
import com.google.common.collect.*
import gnu.trove.impl.Constants
import gnu.trove.map.TObjectIntMap
import gnu.trove.map.hash.TObjectIntHashMap
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.io.CachingIterator
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.checkOrRecalculate
import org.jetbrains.bio.util.div
import org.jgrapht.DirectedGraph
import org.jgrapht.EdgeFactory
import org.jgrapht.event.TraversalListenerAdapter
import org.jgrapht.event.VertexTraversalEvent
import org.jgrapht.experimental.dag.DirectedAcyclicGraph
import org.jgrapht.graph.UnmodifiableDirectedGraph
import org.jgrapht.traverse.DepthFirstIterator
import java.io.BufferedReader
import java.io.IOException
import java.io.Reader
import java.nio.file.Path
import java.util.*
import java.util.concurrent.ConcurrentMap

enum class Ontology {
    BIOLOGICAL_PROCESS,
    CELLULAR_COMPONENT,
    MOLECULAR_FUNCTION;

    private val meta: Pair<DirectedGraph<String, Pair<String, String>>, List<Term>> by lazy {
        val path = Configuration.genomesPath / "go-basic.obo"
        path.checkOrRecalculate("GO") { output ->
            "http://geneontology.org/ontology/go-basic.obo"
                    .downloadTo(output.path)
        }

        OboFile.read(this, path)
    }

    /** A two-letter ontology identifier. */
    val id: String get() = when (this) {
        BIOLOGICAL_PROCESS -> "BP"
        CELLULAR_COMPONENT -> "CC"
        MOLECULAR_FUNCTION -> "MF"
    }

    /** Ontology DAG. */
    val graph: DirectedGraph<String, Pair<String, String>> get() = meta.first

    /** A list of terms in this ontology. */
    val terms: List<Term> get() = meta.second

    /**
     * Returns a mapping of GO terms to gene symbols.
     */
    fun associations(genome: Genome) = CACHE.computeIfAbsent(genome to this) {
        val speciesInformal = when (genome.species) {
            "hg" -> "human"
            "mm" -> "mouse"
            else -> error("unsupported species: ${genome.species}")
        }

        val name = "goa_$speciesInformal.gaf.gz"
        val path = Configuration.genomesPath / name
        path.checkOrRecalculate("GOA") { output ->
            ("http://ftp.ebi.ac.uk/pub/databases/GO/goa/" +
             "${speciesInformal.toUpperCase()}/$name").downloadTo(output.path)
        }

        GafFile.read(genome, path, this)
    }

    companion object {
        private var CACHE: ConcurrentMap<Pair<Genome, Ontology>, SetMultimap<String, String>> =
                Maps.newConcurrentMap()
    }
}

/**
 * Computes a mapping from GO terms to their corresponding depths in
 * the GO graph.
 */
fun DirectedGraph<String, Pair<String, String>>.mapDepth(): TObjectIntMap<String> {
    val depth = TObjectIntHashMap<String>(Constants.DEFAULT_CAPACITY,
                                          Constants.DEFAULT_LOAD_FACTOR,
                                          Int.MAX_VALUE)
    val root = vertexSet().filter { inDegreeOf(it) == 0 }.single()
    val it = DepthFirstIterator(this, root)
    it.addTraversalListener(object : TraversalListenerAdapter<String, Pair<String, String>>() {
        private var current = 0

        override fun vertexTraversed(e: VertexTraversalEvent<String>) {
            current++
        }

        override fun vertexFinished(e: VertexTraversalEvent<String>) {
            depth.put(e.vertex, Math.min(--current, depth[e.vertex]))
        }
    })

    it.forEach {}  // Drain the iterator.
    return depth
}

/**
 * For each term computes the number of genes annotated by it but
 * not any of its parent terms. For instance, here
 *
 *        {a, b, c}
 *            1
 *           / \
 *     {b}  2   3 {c}
 *
 * the gene b counts only towards the term 2.
 *
 * A more elaborate example from the Cao & Zhang paper can be found
 * in the tests.
 *
 * Note that two terms *at the same level* can be annotated by the
 * same set of genes.
 */
fun SetMultimap<String, String>.prune(
        graph: DirectedGraph<String, Pair<String, String>>): SetMultimap<String, String> {
    val pruned = HashMultimap.create<String, String>()

    val root = graph.vertexSet().filter { graph.inDegreeOf(it) == 0 }.single()
    val it = DepthFirstIterator(graph, root)
    it.addTraversalListener(object : TraversalListenerAdapter<String, Pair<String, String>>() {
        /**
         * A mapping from terms to sets of genes, defined as
         *
         *     U[t] = genes annotated by `t` and a union of U[t]
         *            over its children
         *
         * where `t` is a term and `U` is the mapping.
         */
        private val unions = HashMap<String, Set<String>>()

        override fun vertexFinished(e: VertexTraversalEvent<String>) {
            val term = e.vertex
            val seen = HashSet<String>()
            for ((_, child) in graph.outgoingEdgesOf(term)) {
                seen.addAll(unions[child]!!)
            }

            val associated = this@prune[term]!! - seen
            pruned.putAll(term, associated)
            seen.addAll(associated)
            unions[term] = seen
        }
    })

    for (term in it) {}  // Drain the iterator.
    return pruned
}

data class Term(val id: String, val description: String, val isObsolete: Boolean,
                val children: Set<String>) {
    companion object {
        internal fun ofProperties(properties: ListMultimap<String, String>): Term {
            @Suppress("unused_variable") val children = properties["is_a"].map {
                val (id, _) = it.split(" ! ", limit = 2)
                id
            }.toSet()

            return Term(properties.single("id"), properties.single("name"),
                        properties["is_obsolete"] == listOf("true"),
                        children)
        }
    }
}

/**
 * OBO is a generic format for representing ontologies.
 *
 * This class is only suitable for reading Gene Ontology OBO files
 * available from this page: http://geneontology.org/page/download-ontology.
 */
@VisibleForTesting
internal class OboFile constructor(private val ontology: Ontology,
                                   private val reader: Reader) : Iterable<Term> {
    override fun iterator(): Iterator<Term> = OboIterator(ontology, reader.buffered())

    companion object {
        /**
         * Returns a DAG for a given [ontology].
         *
         * The vertices are GO term IDs.
         */
        @Throws(IOException::class)
        fun read(ontology: Ontology, path: Path) = path.bufferedReader().use { reader ->
            val g = DirectedAcyclicGraph(EdgeFactory { u: String, v -> u to v })
            val terms = OboFile(ontology, reader).filterNot(Term::isObsolete)
                    .toList()

            for ((id, _, _, children) in terms) {
                g.addVertex(id)
                for (child in children) {
                    g.addVertex(child)
                    g.addEdge(child, id)
                }
            }

            UnmodifiableDirectedGraph(g) to terms
        }
    }
}

private class OboIterator(private val ontology: Ontology,
                          reader: BufferedReader)
:
        CachingIterator<String, Term>(reader.lines().iterator()) {

    override fun cache(): Term? {
        while (it.hasNext()) {
            if (it.next() == "[Term]") {
                val properties = ArrayListMultimap.create<String, String>()
                while (it.hasNext()) {
                    val line = it.next()
                    if (line.isBlank()) {
                        break
                    }

                    val (key, value) = line.split(": ", limit = 2)
                    properties.put(key, value)
                }

                if (properties.single("namespace").toUpperCase() == ontology.toString()) {
                    return Term.ofProperties(properties)
                }
            }
        }

        return null
    }
}

@Suppress("nothing_to_inline")
private inline fun <K, V> ListMultimap<K, V>.single(key: K): V = this[key].single()

/**
 * GAF is GO annotation format.
 *
 * See http://geneontology.org/page/go-annotation-file-gaf-format-21
 * for format description.
 */
object GafFile {
    private val LOG = Logger.getLogger(GafFile::class.java)

    // See
    // for field descriptions.
    private val FORMAT = CSVFormat.TDF
            .withCommentMarker('!')
            .withHeader("db", "db_object_id", "db_object_symbol", "qualifier",
                        "go_id", "db_reference", "evidence_code", "with_or_from",
                        "aspect", "db_object_name", "db_object_synonym",
                        "do_object_type", "taxon", "date", "assigned_by",
                        "annotation_extension", "gene_product_form_id")

    fun read(genome: Genome, path: Path, ontology: Ontology): SetMultimap<String, String> {
        val aspect = when (ontology) {
            Ontology.BIOLOGICAL_PROCESS -> "P"
            Ontology.CELLULAR_COMPONENT -> "C"
            Ontology.MOLECULAR_FUNCTION -> "F"
        }

        val expected = ontology.graph.vertexSet()
        val acc = HashMultimap.create<String, String>()
        FORMAT.parse(path.bufferedReader()).use {
            var total = 0L
            var unrecognized = 0L

            LOG.info("Reading GO->gene annotations for ${genome.build}")
            for (row in it) {
                if (row["aspect"] != aspect) {
                    continue  // Wrong ontology.
                } else if (row["go_id"] !in expected) {
                    continue  // Deprecated GO term.
                }

                val aliases = (sequenceOf(row["db_object_symbol"]) +
                        row["db_object_synonym"].splitToSequence('|'))
                        .iterator()

                // XXX we can only distinguish between genes with different
                // gene symbols. This is a limitation of the annotations.
                var transcript: Transcript?
                do {
                    transcript = GeneResolver.get(genome.build, aliases.next(),
                                                  GeneAliasType.GENE_SYMBOL)
                            .firstOrNull()
                } while (transcript == null && aliases.hasNext())

                total++
                if (transcript == null) {
                    unrecognized++
                } else {
                    acc.put(row["go_id"], transcript.geneSymbol)
                }
            }

            LOG.info("%.2f%% of records unrecognized".format(unrecognized * 100.0 / total))
        }

        return acc
    }
}
