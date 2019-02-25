package org.jetbrains.bio.util

import org.junit.Test
import java.util.function.Consumer
import kotlin.test.assertNotNull

class TestClass1

class TestClass2 {
    companion object{
        @Suppress("unused")
        @JvmStatic val INSTANCE : TestClass2 = TestClass2()
    }
}

class ClassesTest {

    @Test fun initializeTest() {
        assertNotNull(ClassProcessor.tryToInstantiate(TestClass1::class.java))
        assertNotNull(ClassProcessor.tryToInstantiate(TestClass2::class.java))
    }

    @Test fun classProcessor() {
        var junitSeen = false
        var classesTestSeen = false
        ClassProcessor.processClasses(Consumer<String> { t ->
            if (t.startsWith("org.junit")) {
                junitSeen = true
            }
            if (t == ClassesTest::class.qualifiedName) {
                classesTestSeen = true
            }
        })
        assert(junitSeen)
        assert(classesTestSeen)
    }
}