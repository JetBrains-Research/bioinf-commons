package org.jetbrains.bio.genome.containers

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.util.await
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.ConcurrentHashMap

/**
 * Creates a new genome map from a given [init] function.
 *
 * If [parallel] is `false`, initialization is performed
 * sequentially otherwise for each chromosome [init] is
 * called in a separate thread.
 */
fun <T> genomeMap(
        genomeQuery: GenomeQuery,
        parallel: Boolean = false,
        init: (Chromosome) -> T
) = GenomeMap(genomeQuery, parallel, init)

fun <T> genomeStrandMap(
        genomeQuery: GenomeQuery,
        parallel: Boolean = false,
        init: (Chromosome, Strand) -> T
) = GenomeStrandMap(genomeQuery, parallel, init)

/**
 * A map with a fixed set of keys, defined by [GenomeQuery].
 *
 * All map operations are thread-safe.
 *
 * @author Sergei Lebedev
 */
class GenomeMap<T> internal constructor(
        val genomeQuery: GenomeQuery,
        parallel: Boolean,
        f: (Chromosome) -> T
) {

    private val data: ConcurrentHashMap<String, T> = ConcurrentHashMap()

    init {
        genomeQuery.get().map { chromosome ->
            Callable {
                data[chromosome.name] = f(chromosome)
            }
        }.await(parallel)
    }

    operator fun contains(chromosome: Chromosome) =
            chromosome.genome.build == genomeQuery.build && chromosome.name in data.keys


    operator fun get(chromosome: Chromosome): T {
        if (chromosome !in this) {
            throw NoSuchElementException(chromosome.toString())
        }
        return data[chromosome.name]!!
    }

    operator fun set(chromosome: Chromosome, value: T) {
        if (chromosome !in this) {
            throw NoSuchElementException(chromosome.toString())
        }
        data[chromosome.name] = value
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is GenomeMap<*>) return false

        if (genomeQuery != other.genomeQuery) return false
        if (data != other.data) return false

        return true
    }

    override fun hashCode(): Int {
        var result = genomeQuery.hashCode()
        result = 31 * result + data.hashCode()
        return result
    }

}

/**
 * A [GenomeMap] interface where each chromosome is allowed to have
 * per [Strand] values.
 *
 * All map operations are thread-safe.
 *
 * @author Roman Chernyatchik
 */
interface GenomeStrandMapLike<T> {
    val genomeQuery: GenomeQuery
    operator fun get(chromosome: Chromosome, strand: Strand): T

    fun asSequence(): Sequence<T> = genomeQuery.get().asSequence().flatMap { chromosome ->
        Strand.values().asSequence().map { strand -> this[chromosome, strand] }
    }
}

/**
 * A [GenomeMap] where each chromosome is allowed to have
 * per [Strand] values.
 *
 * All map operations are thread-safe.
 *
 * @author Sergei Lebedev
 */
class GenomeStrandMap<T> internal constructor(
        override val genomeQuery: GenomeQuery,
        parallel: Boolean,
        f: (Chromosome, Strand) -> T)
    : GenomeStrandMapLike<T> {

    private val data: ConcurrentHashMap<Pair<String, Strand>, T> = ConcurrentHashMap()

    init {
        genomeQuery.get().flatMap { chromosome ->
            Strand.values().map { strand ->
                Callable {
                    data[chromosome.name to strand] = f(chromosome, strand)
                }
            }
        }.await(parallel)
    }

    fun contains(chromosome: Chromosome, strand: Strand) =
            chromosome.genome.build == genomeQuery.build && chromosome.name to strand in data.keys

    override operator fun get(chromosome: Chromosome, strand: Strand): T {
        if (!contains(chromosome, strand)) {
            throw NoSuchElementException("${chromosome.name} $strand")
        }
        return data[chromosome.name to strand]!!
    }

    operator fun set(chromosome: Chromosome, strand: Strand, value: T) {
        if (!contains(chromosome, strand)) {
            throw NoSuchElementException("${chromosome.name} $strand")
        }
        data[chromosome.name to strand] = value
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is GenomeStrandMap<*>) return false

        if (genomeQuery != other.genomeQuery) return false
        if (data != other.data) return false

        return true
    }

    override fun hashCode(): Int {
        var result = genomeQuery.hashCode()
        result = 31 * result + data.hashCode()
        return result
    }
}