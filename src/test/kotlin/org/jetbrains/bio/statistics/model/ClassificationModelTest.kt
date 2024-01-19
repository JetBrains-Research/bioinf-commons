package org.jetbrains.bio.statistics.model

import com.google.gson.GsonBuilder
import com.google.gson.JsonParseException
import org.jetbrains.bio.Tests
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.util.*
import org.jetbrains.bio.viktor.F64Array
import org.junit.Before
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals

class ClassificationModelTest {
    @Before
    fun setUp() {
        Boo.VERSION = 222
    }

    @Test
    fun testSave() = withTempFile("model", ".json") { path ->
        Boo(1).save(path)
        assertEquals(
            """{
  "value": 1,
  "b": "bbb",
  "name": "bboo",
  "model.class.fqn": "org.jetbrains.bio.statistics.model.Boo",
  "model.class.format": 222
}""", path.read()
        )
    }

    @Test
    fun testLoad() = withTempFile("model", ".json") { path ->
        val obj = Boo(1)
        obj.save(path)
        val boo = ClassificationModel.load<Boo>(path)
        assertEquals(obj, boo)
    }

    @Test
    fun testLoadNoVersion() {
        withTempFile("model", ".json") { path ->
            val gson = GsonBuilder().setPrettyPrinting()
                .serializeSpecialFloatingPointValues()
                .create()

            val obj = Boo(1)
            path.bufferedWriter().use { gson.toJson(obj, it) }

            Tests.assertThrowsWithMessage(
                "Deserialization error: Class name (model.class.fqn) is missing.",
                JsonParseException::class.java,
            ) {
                ClassificationModel.load<Boo>(path)
            }
        }
    }

    @Test
    fun testSaveNaN() = withTempFile("model", ".json") { path ->
        val nanBoo = NanBoo()
        nanBoo.save(path)
        assertEquals(
            """{
  "value": NaN,
  "name": "NaN-boo",
  "model.class.fqn": "org.jetbrains.bio.statistics.model.NanBoo",
  "model.class.format": 0
}""",
            path.read()
        )

        assertEquals(nanBoo, ClassificationModel.load<NanBoo>(path))
    }

    @Test
    fun testLoadWrongVersion1() {
        withTempFile("model", ".json") { path ->
            val obj = Boo(1)
            obj.save(path)

            // Change model
            Boo.VERSION = 123

            Tests.assertThrowsWithMessage(
                "Deserialization error: Format has changed, " +
                        "'org.jetbrains.bio.statistics.model.Boo' expects '123' version, but got '222'",
                JsonParseException::class.java,
            ) {
                ClassificationModel.load<Boo>(path)
            }
        }
    }

    @Test
    fun testLoadWrongVersion2() {
        withTempFile("model", ".json") { path ->
            val obj = Boo(1)

            obj.save(path)

            // Change model
            Boo.VERSION = 666

            Tests.assertThrowsWithMessage(
                "Deserialization error: Format has changed, " +
                        "'org.jetbrains.bio.statistics.model.Boo' expects '666' version, but got '222'",
                JsonParseException::class.java,
            ) {
                ClassificationModel.load<Boo>(path)
            }
        }
    }
}

abstract class AbstractBoo(private val name: String) : ClassificationModel {
    override fun degreesOfFreedom() = 0

    override fun fit(
        preprocessed: Preprocessed<DataFrame>,
        title: String, threshold: Double, maxIterations: Int
    ) {
    }

    override fun predict(preprocessed: Preprocessed<DataFrame>) = IntArray(0)

    override fun logLikelihood(preprocessed: Preprocessed<DataFrame>) = 0.0

    override fun toString() = name
}

class Boo(private val value: Int) : AbstractBoo("bboo") {
    override fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array = TODO()

    private val b = "bbb"

    override fun toString(): String {
        return "%s, b = %s, value = %d".format(super.toString(), b, value)
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is Boo) return false

        if (value != other.value) return false
        return b == other.b

    }

    override fun hashCode() = Objects.hash(b, value)

    companion object {
        @Transient
        @JvmField
        var VERSION = 222
    }
}

class NanBoo : AbstractBoo("NaN-boo") {
    override fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array = TODO()

    private val value = java.lang.Double.NaN

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is NanBoo) return false

        return other.value.compareTo(value) == 0

    }

    override fun hashCode() = value.hashCode()

    companion object {
        @Suppress("unused")
        @Transient
        const val VERSION = 0
    }
}
