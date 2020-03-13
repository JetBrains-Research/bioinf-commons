package org.jetbrains.bio.genome.query

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.jetbrains.bio.util.LockManager
import java.util.concurrent.locks.ReentrantLock
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

/** An query which caches the result of [apply]. */
abstract class CachingQuery<I : Any, O> : Query<I, O> {

    abstract fun getUncached(input: I): O

    private val cache: Cache<I, O> = CacheBuilder.newBuilder().build()
    private val lock = ReentrantLock()

    override fun apply(t: I): O = LockManager.synchronized(t) {
        check(lock.holdCount <= 0) { "Attempt to call CachingInputQuery#apply recursively" }
        return@synchronized cache.get(t) { getUncached(t) }
    }
}