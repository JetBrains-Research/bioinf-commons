package org.jetbrains.bio.query

import org.jetbrains.bio.util.LockManager
import java.lang.ref.SoftReference
import java.util.concurrent.locks.ReentrantLock
import java.util.function.Supplier

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
    private var cached = SoftReference<T>(null)
    private val lock = ReentrantLock()

    override fun get() = LockManager.synchronized(this) {
        check(lock.holdCount <= 0) { "Attempt to call CachingInputQuery#get recursively" }
        var value = cached.get()
        if (value == null) {
            value = getUncached()
            cached = SoftReference(value)
        }

        value!!    // Kotlin 1.1-rc workaround (KT-16368)
    }
}