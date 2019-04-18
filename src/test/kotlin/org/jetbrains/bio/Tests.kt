package org.jetbrains.bio

import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * In my sincerest opinion, [assertTrue] without a custom message should not be used
 * under any circumstances. There's nothing more annoying than seeing a test fail
 * with an oh-so-informative message of "Expected the value to be true." Yeah, and I
 * expected the message to help me identify the problem. Guess we're both disappointed.
 *
 * So, I've collected a few assert helper methods which cover popular use cases of [assertTrue].
 */
object Tests {

    fun assertIn(substring: String, fullString: String) {
        // Process Windows with different line separators correctly.
        substring.lines().forEach { s ->
            assertTrue(s in fullString, "Expected <$s> to be in <$fullString>.")
        }
    }

    fun assertIn(regex: Regex, fullString: String) {
        // Process Windows with different line separators correctly.
        assertTrue(regex.matches(fullString), "Expected <$regex> to match <$fullString>.")
    }

    fun <T> assertIn(actual: T, expected: List<T>) = assertTrue(
        actual in expected,
        "Expected <$actual> to be in $expected."
    )

    fun <T> assertNotIn(actual: T, expected: List<T>) = assertFalse(
        actual in expected,
        "Expected <$actual> not to be in $expected."
    )

    fun assertIs(actual: Any, expected: Class<out Any>) {
        assertTrue(expected.isInstance(actual),
            "Expected ${expected.simpleName}, got ${actual::class.java.simpleName}.")
    }

    fun assertMatches(output: String, regex: Regex) = assertTrue(
        regex.matches(output),
        "Regex ${regex.pattern} doesn't match content:\n<$output>"
    )

}