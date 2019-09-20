package org.jetbrains.bio.statistics

import com.google.gson.GsonBuilder
import com.google.gson.JsonParseException
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.read
import org.jetbrains.bio.util.withTempFile
import org.jetbrains.bio.viktor.F64Array
import org.junit.Before
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
import java.util.*
import kotlin.test.assertEquals

class ClassificationModelTest {
    @Rule @JvmField val thrown = ExpectedException.none()

    @Before fun setUp() {
        Boo.VERSION = 222
    }

    @Test fun testSave() = withTempFile("model", ".json") { path ->
        Boo(1).save(path)
        assertEquals("{\n  \"b\": \"bbb\",\n" +
                     "  \"value\": 1,\n" +
                     "  \"name\": \"bboo\",\n" +
                     "  \"model.class.fqn\": \"org.jetbrains.bio.statistics.Boo\",\n" +
                     "  \"model.class.format\": \"222\"\n}",
                path.read())
    }

    @Test fun testLoad() = withTempFile("model", ".json") { path ->
        val obj = Boo(1)
        obj.save(path)
        val boo = ClassificationModel.load<Boo>(path)
        assertEquals(obj, boo)
    }

    @Test fun testLoadNoVersion() {
        withTempFile("model", ".json") { path ->
            val gson = GsonBuilder().setPrettyPrinting()
                    .serializeSpecialFloatingPointValues()
                    .create()

            val obj = Boo(1)
            path.bufferedWriter().use { gson.toJson(obj, it) }

            thrown.expect(JsonParseException::class.java)
            thrown.expectMessage("Deserialization error: Class name (model.class.fqn) is missing.")
            ClassificationModel.load<Boo>(path)
        }
    }

    @Test fun testSaveNaN() = withTempFile("model", ".json") { path ->
        val nanBoo = NanBoo()
        nanBoo.save(path)
        assertEquals("{\n  \"value\": NaN,\n  " +
                     "\"name\": \"NaN-boo\",\n  " +
                     "\"model.class.fqn\": \"org.jetbrains.bio.statistics.NanBoo\",\n  " +
                     "\"model.class.format\": \"0\"\n}",
                path.read())

        assertEquals(nanBoo, ClassificationModel.load<NanBoo>(path))
    }

    @Test fun testLoadWrongVersion1() {
        withTempFile("model", ".json") { path ->
            val obj = Boo(1)
            obj.save(path)

            // Change model
            Boo.VERSION = 123

            thrown.expect(JsonParseException::class.java)
            thrown.expectMessage("Deserialization error: Format has changed, " +
                                 "'org.jetbrains.bio.statistics.Boo' expects '123' version, but got '222'")
            ClassificationModel.load<Boo>(path)
        }
    }

    @Test fun testLoadWrongVersion2() {
        withTempFile("model", ".json") { path ->
            val obj = Boo(1)

            obj.save(path)

            // Change model
            Boo.VERSION = 666

            thrown.expect(JsonParseException::class.java)
            thrown.expectMessage("Deserialization error: Format has changed, " +
                                 "'org.jetbrains.bio.statistics.Boo' expects '666' version, but got '222'")

            ClassificationModel.load<Boo>(path)
        }
    }
}

abstract class AbstractBoo(private val name: String) : ClassificationModel {
    override fun degreesOfFreedom() = 0

    override fun fit(preprocessed: Preprocessed<DataFrame>,
                     title: String, threshold: Double, maxIter: Int) {
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
        @Transient @JvmField var VERSION = 222
    }
}

class NanBoo : AbstractBoo("NaN-boo") {
    override fun evaluate(preprocessed: Preprocessed<DataFrame>): F64Array = TODO()

    private val value = java.lang.Double.NaN

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is NanBoo) return false

        return java.lang.Double.compare(other.value, value) == 0

    }

    override fun hashCode() = value.hashCode()

    companion object {
        @Suppress("unused")
        @Transient @JvmField val VERSION = 0
    }
}
