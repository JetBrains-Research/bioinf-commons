package org.jetbrains.bio.statistics.gson

import com.google.gson.*
import com.google.gson.reflect.TypeToken
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.jetbrains.bio.viktor.toF64Array
import java.lang.reflect.Type

/**
 * Historically strided containers were serialized *as-is*, i.e. each
 * container was a JSON object. This function is to support this legacy
 * format.
 *
 * Creating [com.google.gson.Gson] on each call is obviously a hack,
 * but the API for calling default deserializer seems to be missing.
 */
private inline fun <reified T> fromJsonFallback(json: JsonElement): T {
    return GsonBuilder()
        .setFieldNamingStrategy(GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY)
        .create()
        .fromJson(json, object : TypeToken<T>() {}.type)
}

object F64ArrayTypeAdapter : JsonSerializer<F64Array>, JsonDeserializer<F64Array> {
    override fun serialize(
        src: F64Array, typeOfSrc: Type,
        context: JsonSerializationContext
    ): JsonElement {
        return context.serialize(src.toArray())
    }

    override fun deserialize(
        json: JsonElement, typeOfT: Type,
        context: JsonDeserializationContext
    ): F64Array {
        if (json.isJsonObject) {
            val a = fromJsonFallback<F64Array>(json)
            return a
        }

        // XXX there is no way to deserialize an array of arbitrary
        //     depth in Gson, therefore we special case.
        return when (json.guessDepth()) {
            1 -> context.deserialize<DoubleArray>(
                json, object : TypeToken<DoubleArray>() {}.type
            ).asF64Array()
            2 -> context.deserialize<Array<DoubleArray>>(
                json, object : TypeToken<Array<DoubleArray>>() {}.type
            ).toF64Array()
            3 -> context.deserialize<Array<Array<DoubleArray>>>(
                json, object : TypeToken<Array<Array<DoubleArray>>>() {}.type
            ).toF64Array()
            else -> throw IllegalStateException()
        }
    }
}

private fun JsonElement.guessDepth(): Int {
    val item = asJsonArray.first()
    return when {
        item.isJsonPrimitive -> 1
        item.isJsonArray -> 1 + item.guessDepth()
        else -> throw IllegalStateException()
    }
}