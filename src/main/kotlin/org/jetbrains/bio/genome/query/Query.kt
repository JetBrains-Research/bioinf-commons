package org.jetbrains.bio.genome.query

import com.github.benmanes.caffeine.cache.Cache
import com.github.benmanes.caffeine.cache.Caffeine
import java.util.function.Function

/**
 * Named [Function] implementation capable of caching the results.
 */
interface Query<I, O> : Function<I, O> {
    /**
     * Returns unique identifier suitable for use in a file name.
     */
    val id: String

    /**
     * Returns a human-readable description for this query.
     */
    val description: String get() = id
}

/** A query which caches the result of [apply]. */
abstract class CachingQuery<I : Any, O> : Query<I, O> {

    abstract fun getUncached(input: I): O

    private val cache: Cache<I, O> =  Caffeine.newBuilder().build()

    override fun apply(t: I): O = cache.get(t) { getUncached(t) }
}