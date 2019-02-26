package org.jetbrains.bio.statistics

import org.apache.commons.math3.distribution.AbstractRealDistribution
import org.apache.commons.math3.stat.inference.ChiSquareTest
import org.junit.Assert.assertTrue

/**
 * [org.junit.Assert] in a probabilistic way!
 *
 * @author Alexey Dievsky
 * @author Dmitry Groshev
 * @since 24/07/14
 */

/**
 * Asserts that state transitions estimated from given `states`
 * are consistent with given transition probabilities using a
 * chi^2-GOF-test.
 */
fun assertIsConsistent(states: IntArray,
                       transitionProbabilities: Array<*>,
                       threshold: Double) {
    val numStates = transitionProbabilities.size
    val transitions = Array(numStates) { LongArray(numStates) }
    for (t in 0 until states.size - 1) {
        val i = states[t]
        val j = states[t + 1]
        transitions[i][j]++
    }

    for (i in 0 until numStates) {
        val pValue = ChiSquareTest()
                .chiSquareTest(transitionProbabilities[i] as DoubleArray, transitions[i])
        val message = "chi^2-test failed for transitions from the state %d: %f < %f"
        assertTrue(String.format(message, i, pValue, threshold), pValue > threshold)
    }
}

/**
 * Asserts tha the probability P(|X - EX| > |value - EX|) is less
 * than `threshold`, where X follows a given `distribution`.
 */
@JvmOverloads
fun assertTwoSidedEquals(value: Double,
                         distribution: AbstractRealDistribution,
                         threshold: Double = 1e-3) {
    val mean = distribution.numericalMean
    val delta = Math.abs(value - mean)
    val testStatistic = (1 - distribution.cumulativeProbability(mean + delta)
            + distribution.cumulativeProbability(mean - delta))
    val message = "The probability of seeing X = %f (EX = %f) is too low: %f < %f"
    assertTrue(String.format(message, value, mean, testStatistic, threshold),
            testStatistic > threshold)
}