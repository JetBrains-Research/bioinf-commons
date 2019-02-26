package org.jetbrains.bio.statistics.gson

import com.google.gson.Gson
import com.google.gson.TypeAdapter
import com.google.gson.TypeAdapterFactory
import com.google.gson.reflect.TypeToken
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter

/**
 * This interface describes a class that cannot be easily correctly
 * deserialized, typically because it contains useful transients which
 * are not filled during deserialization.
 *
 * The interface provides a method [.updateTransients] which should be
 * called after deserialization (i.e. after the non-transient fields are
 * filled) and a [com.google.gson.TypeAdapterFactory] to be registered
 * in the reader which ensures that the method will be called properly
 * (otherwise delegating reading and writing to default adapters).
 *
 * Take care to register this factory before any other factory that could
 * process an implementation of the interface.
 *
 * @author Alexey Dievsky
 * @since 16/04/15
 */
interface NotDirectlyDeserializable {
    /**
     * Fill the transient fields left blank by GSON deserializer
     * (and possibly do other actions to leave the object in a consistent state).
     */
    fun updateTransients()

    companion object {
        val ADAPTER_FACTORY: TypeAdapterFactory = object : TypeAdapterFactory {
            override fun <T> create(gson: Gson, type: TypeToken<T>): TypeAdapter<T>? {
                if (!NotDirectlyDeserializable::class.java.isAssignableFrom(type.rawType)) {
                    return null
                }
                val delegateAdapter = gson.getDelegateAdapter(this, type)
                return object : TypeAdapter<T>() {
                    override fun write(out: JsonWriter, value: T) = delegateAdapter.write(out, value)

                    override fun read(`in`: JsonReader): T {
                        val res = delegateAdapter.read(`in`)
                        (res as NotDirectlyDeserializable).updateTransients()
                        return res
                    }
                }
            }
        }
    }

}
