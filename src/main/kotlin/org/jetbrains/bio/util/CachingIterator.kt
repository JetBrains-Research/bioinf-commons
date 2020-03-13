package org.jetbrains.bio.util

import com.google.common.collect.Iterators
import com.google.common.collect.PeekingIterator

/** A base class for lazy I/O-based iterators. */
abstract class CachingIterator<S, out T>(wrapped: Iterator<S>) : Iterator<T> {
    protected var it: PeekingIterator<S> = Iterators.peekingIterator(wrapped)
    private var cached: T? = null

    override fun hasNext(): Boolean {
        if (cached == null && it.hasNext()) {
            cached = cache()  // Got some?
        }

        return cached != null
    }

    override fun next(): T {
        check(hasNext())
        val next = cached
        cached = null
        return next!!
    }

    protected abstract fun cache(): T?
}