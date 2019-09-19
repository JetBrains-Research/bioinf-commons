package org.jetbrains.bio.statistics.gson

import com.google.common.base.CaseFormat
import com.google.gson.*
import com.google.gson.internal.Streams
import com.google.gson.internal.bind.JsonTreeReader
import com.google.gson.internal.bind.JsonTreeWriter
import com.google.gson.reflect.TypeToken
import com.google.gson.stream.JsonReader
import com.google.gson.stream.JsonWriter
import java.io.IOException

/**
 * @author Alexey Dievsky
 * @since 03/07/15
 */
object GSONUtil {

    /**
     * This naming strategy converts 'myFooBar' field into 'foo_bar'.
     */
    @JvmField
    val NO_MY_UNDESCORE_NAMING_STRATEGY = FieldNamingStrategy { f ->
        CaseFormat.LOWER_CAMEL.to(CaseFormat.LOWER_UNDERSCORE, f.name).replace("my_", "")
    }


    /**
     * This factory checks the incoming object type against the [aClass].
     * If it's a descendant, the factory creates and returns an adapter
     * via [typeAdapter] function, otherwise it returns null.
     */
    @JvmStatic fun <C> classSpecificFactory(
            aClass: Class<C>,
            typeAdapter: (Gson, TypeAdapterFactory) -> TypeAdapter<C>): TypeAdapterFactory {
        return object : TypeAdapterFactory {
            override fun <T> create(gson: Gson, type: TypeToken<T>): TypeAdapter<T>? {
                return if (aClass.isAssignableFrom(type.rawType)) {
                    (typeAdapter(gson, this)) as TypeAdapter<T>
                } else {
                    null
                }
            }
        }
    }

    /**
     * This adapter is used to work with polymorphic types. It saves the
     * fully qualified runtime class name under the [fqnField] entry and,
     * respectively, can read it from there and return an appropriate
     * type object.
     */
    @JvmStatic fun <T: Any> classAwareAdapter(
            gson: Gson, factory: TypeAdapterFactory, fqnField: String): TypeAdapter<T> {
        return object : TypeAdapter<T>() {
            @Throws(IOException::class)
            override fun write(out: JsonWriter, value: T) {

                val fqn = value.javaClass.name

                val delegate = gson.getDelegateAdapter(factory, TypeToken.get(value.javaClass))

                val tmpWriter = JsonTreeWriter()
                tmpWriter.isLenient = true
                delegate.write(tmpWriter, value)
                val element = tmpWriter.get().asJsonObject

                element.add(fqnField, JsonPrimitive(fqn))
                Streams.write(element, out)
            }

            @Throws(IOException::class)
            override fun read(reader: JsonReader): T? {
                val element = Streams.parse(reader)
                val fqnElement = element.asJsonObject.remove(fqnField)
                if (fqnElement == null) {
                    deserializationError("Class name (%s) field is missing." +
                                         " Please, recalculate the model.", fqnField)
                }
                val fqn = fqnElement!!.asString

                val aClass: Class<T>
                try {
                    aClass = Class.forName(fqn) as Class<T>
                } catch (e: ClassNotFoundException) {
                    deserializationError("Cannot load class %s", fqn, cause = e)
                    return null
                }


                val adapter = gson.getDelegateAdapter(factory, TypeToken.get(aClass))

                val tmpReader = JsonTreeReader(element)
                tmpReader.isLenient = true
                return adapter.read(tmpReader)
            }

        }
    }

    /**
     * Same as [classAwareAdapter], but also saves the contents of a static field named "VERSION" and later
     * checks it against the actual contents of this field when loading.
     * Helps to work around the changes in the class field layout;
     * a version change is more self-explanatory than a random field reading failure.
     */
    @JvmStatic fun <T : Any> classAndVersionAdapter(gson: Gson, factory: TypeAdapterFactory,
                                                    fqnField: String, versionField: String): TypeAdapter<T> {
        return object : TypeAdapter<T>() {
            override fun write(out: JsonWriter, value: T) {
                val modelFQN = value.javaClass.name
                val modelVersion = getSerializationFormatVersion(value.javaClass)

                val delegate = gson.getDelegateAdapter(factory, TypeToken.get(value.javaClass))

                val tmpWriter = JsonTreeWriter()
                tmpWriter.isLenient = true
                delegate.write(tmpWriter, value)
                val element = tmpWriter.get().asJsonObject

                element.add(fqnField, JsonPrimitive(modelFQN))
                element.add(versionField, JsonPrimitive(modelVersion))
                Streams.write(element, out)
            }

            override fun read(`in`: JsonReader): T? {
                val element = Streams.parse(`in`)
                val modelFQNElement = element.asJsonObject.remove(fqnField)
                val modelFormatVersElement = element.asJsonObject.remove(versionField)
                if (modelFQNElement == null) {
                    deserializationError("Class name ($fqnField) is missing.")
                }
                if (modelFormatVersElement == null) {
                    deserializationError("Version field ($versionField) is missing.")
                }
                val fqn = modelFQNElement!!.asString
                val vers = modelFormatVersElement!!.asString

                val aClass: Class<T>
                try {
                    aClass = Class.forName(fqn) as Class<T>
                } catch (e: ClassNotFoundException) {
                    deserializationError("Cannot load class %s", fqn, cause = e)
                    return null
                }


                val recentVersion = getSerializationFormatVersion(aClass)
                if (recentVersion != vers) {
                    deserializationError("Format has changed, '%s' expects '%s' version, but got '%s'",
                                         fqn, recentVersion, vers)
                }

                val adapter = gson.getDelegateAdapter(factory, TypeToken.get(aClass))

                val tmpReader = JsonTreeReader(element)
                tmpReader.isLenient = true
                return adapter.read(tmpReader)
            }

            private fun getSerializationFormatVersion(aClass: Class<T>): String {
                try {
                    return (aClass.getDeclaredField("VERSION").get(null)).toString()
                } catch (e: Exception) {
                    deserializationError("Cannot get serialization format version." +
                                         " Probably VERSION field is missing in %s. Exception message: %s",
                                         aClass.name, e.message, cause = e)
                    return "" // this statement can never be reached
                }
            }
        }
    }

    private fun deserializationError(format: String, vararg args: Any?, cause: Exception? = null) {
        throw if (cause == null) {
            JsonParseException("Deserialization error: " + format.format(*args))
        } else {
            JsonParseException("Deserialization error: " + format.format(*args), cause)
        }
    }
}
