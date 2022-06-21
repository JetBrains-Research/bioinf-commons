package org.jetbrains.bio.genome.query

import java.lang.ref.SoftReference
import java.util.function.Supplier

/**
 * Named [Supplier] implementation capable of caching the result.
 */
interface InputQuery<T> : Supplier<T> {
    fun getUncached(): T

    override fun get() = getUncached()  // no caching in default implementation.

    /**
     * A unique identifier suitable for use in a file name.
     */
    val id: String

    /**
     * A human-readable description for this query.
     */
    val description: String get() = id
}

/** An input query which caches the result of [get] in a [SoftReference]. */
abstract class CachingInputQuery<T> : InputQuery<T> {

    private val cachedValue: T by lazy { getUncached() }
    override fun get() = cachedValue
}