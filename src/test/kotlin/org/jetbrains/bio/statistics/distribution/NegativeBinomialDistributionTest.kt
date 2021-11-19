package org.jetbrains.bio.statistics.distribution

import gnu.trove.list.array.TIntArrayList
import junit.framework.AssertionFailedError
import junit.framework.TestCase.assertEquals
import org.jetbrains.bio.Retry
import org.jetbrains.bio.RetryRule
import org.jetbrains.bio.viktor.F64Array
import org.junit.Before
import org.junit.Rule
import org.junit.Test
import kotlin.test.assertTrue

class NegativeBinomialDistributionTest {
    @get:Rule
    var rule = RetryRule(3)

    @Before
    fun setUp() {
        Sampling.RANDOM_DATA_GENERATOR.reSeed(42)
    }

    @Test
    fun testMeanGreaterThanVariance() {
        assertTrue(NegativeBinomialDistribution.estimateFailuresUsingMoments(1.0, 0.5).isInfinite())
    }

    @Test
    fun testMeanVarianceZeros() {
        assertTrue(NegativeBinomialDistribution.estimateFailuresUsingMoments(0.0, 0.0).isInfinite())
    }

    @Retry
    @Test
    fun testMLE() {
        var numFailed = 0
        val numRuns = 5
        val sampleSize = 65536

        for (run in 0 until numRuns) {
            val numFailures = Sampling.sampleUniform(0.0, 42.0).toInt()
            val successProbability = Sampling.sampleUniform(0.0, 1.0)
            val nbSampling =
                NegativeBinomialDistribution.usingSuccessProbability(successProbability, numFailures.toDouble())
            val observations = TIntArrayList(nbSampling.sample(sampleSize))

            val data = observations.toArray()
            try {
                val nb = NegativeBinomialDistribution.of(data)
                assertEquals(numFailures.toDouble(), nb.failures, 1.0)
                assertEquals(successProbability, nb.probabilityOfSuccess, 0.05)
            } catch (ignored: AssertionFailedError) {
                numFailed++
            }

        }

        // Note(lebedev): we're happy if at most 1/5 runs fail.
        assertTrue(numFailed <= 1)
    }

    @Retry
    @Test
    fun testMLEF64Array() {
        var numFailed = 0
        val numRuns = 5
        val sampleSize = 65536

        for (run in 0 until numRuns) {
            val numFailures = Sampling.sampleUniform(1.0, 42.0).toInt()
            val successProbability = Sampling.sampleUniform(0.0, 1.0)
            val nbSampling =
                NegativeBinomialDistribution.usingSuccessProbability(successProbability, numFailures.toDouble())
            val samples = TIntArrayList(nbSampling.sample(sampleSize))
            val weights = F64Array(sampleSize) { 1.0 }
            try {
                val samplesArray = samples.toArray()
                val mean = samplesArray.average()
                val failuresGuess = Sampling.sampleUniform(1.0, 42.0)
                val failures =
                    NegativeBinomialDistribution.fitNumberOfFailures(samplesArray, weights, mean, failuresGuess)
                assertEquals(numFailures.toDouble(), failures, 1.0)
                assertEquals(successProbability, mean / (mean + failures), 0.05)
            } catch (ignored: AssertionFailedError) {
                numFailed++
            }

        }

        // Note(lebedev): we're happy if at most 1/5 runs fail.
        assertTrue(numFailed <= 1)
    }

    @Retry
    @Test
    fun testFitGamma() {
        val sampleSize = 65536
        val numRuns = 5
        var numFailed = 0
        for (i in 0.until(numRuns)) {
            val shape = Sampling.sampleUniform(0.0, 20.0)
            val rate = Sampling.sampleUniform(0.0, 2.0)
            val samples = DoubleArray(sampleSize) { Sampling.sampleGamma(shape, rate) }
            val logSamples = DoubleArray(sampleSize) { Math.log(samples[it]) }
            val shapeGuess = Sampling.sampleUniform(0.0, 20.0)
            val fittedShape = NegativeBinomialDistribution.fitGamma(logSamples.average(), samples.average(), shapeGuess)
            try {
                assertEquals(shape, fittedShape, 0.1)
            } catch (ignored: AssertionFailedError) {
                numFailed++
            }
        }
        assertTrue(numFailed <= 1)
    }

    /**
     * All the values generated with R.
     */
    @Test
    fun testFixedValues() {
        val failures = doubleArrayOf(
            10.344394, 6.636840, 8.272800, 5.790573, 15.066386, 8.219735, 7.611640,
            5.793857, 8.699888, 14.319623
        )
        val successProbabilities = doubleArrayOf(
            0.60400253, 0.65120170, 0.03357039, 0.20654823, 0.60348266,
            0.54437009, 0.57465349, 0.60595111, 0.39196646, 0.29655913
        )
        val means = doubleArrayOf(
            15.7779795, 12.3908904, 0.2873682, 1.5073790, 22.9304040, 9.8206412,
            10.2835114, 8.9095387, 5.6083489, 6.0369179
        )
        val sample = intArrayOf(1, 6, 5, 0, 4, 2, 4, 2, 4, 6)
        val probabilities = doubleArrayOf(
            4.306413e-04, 5.092278e-02, 2.920415e-05, 2.619186e-01,
            3.653011e-04, 1.754760e-02, 4.536087e-02, 3.277967e-02,
            1.370312e-01, 1.345970e-01
        )
        val logProbabilities = doubleArrayOf(
            -7.750235, -2.977445, -10.441200, -1.339722, -7.914789,
            -4.042838, -3.093105, -3.417947, -1.987547, -2.005470
        )
        val eps = 1E-6
        for (i in 0..9) {
            val fromMean = NegativeBinomialDistribution.usingMean(means[i], failures[i])
            val fromProbability = NegativeBinomialDistribution.usingSuccessProbability(
                successProbabilities[i], failures[i]
            )
            assertEquals(fromMean.numericalMean, fromProbability.numericalMean, eps)
            assertEquals(fromMean.probabilityOfSuccess, fromProbability.probabilityOfSuccess, eps)
            assertEquals(fromMean.failures, fromProbability.failures, eps)
            assertEquals(probabilities[i], fromMean.probability(sample[i]), eps)
            assertEquals(probabilities[i], fromProbability.probability(sample[i]), eps)
            assertEquals(logProbabilities[i], fromMean.logProbability(sample[i]), eps)
            assertEquals(logProbabilities[i], fromProbability.logProbability(sample[i]), eps)
            assertEquals(
                logProbabilities[i],
                Densities.logNegativeBinomialDensity(sample[i], means[i], failures[i]), eps
            )
        }
    }
}