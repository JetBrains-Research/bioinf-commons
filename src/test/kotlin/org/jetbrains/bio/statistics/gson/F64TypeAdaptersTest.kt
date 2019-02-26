package org.jetbrains.bio.statistics.gson

import com.google.common.reflect.TypeToken
import com.google.gson.GsonBuilder
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.test.assertEquals

private val GSON_BUILDER = GsonBuilder()
        .setFieldNamingStrategy(GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY)

class F64ArrayTypeAdapterTest {
    @Test fun toJsonFromJson() {
        val gson = GSON_BUILDER
                .registerTypeAdapter(F64Array::class.java, F64ArrayTypeAdapter)
                .create()

        val token = object : TypeToken<F64Array>() {}.type
        val m = F64Array(3, 4) { i, j -> Math.pow(i.toDouble(), j.toDouble()) }

        // as array
        assertEquals(m, gson.fromJson(gson.toJson(m.toArray()), token))
        assertEquals(m, gson.fromJson(gson.toJson(m), token))
        // as object
        assertEquals(m, gson.fromJson(GSON_BUILDER.create().toJson(m), token))
    }
}