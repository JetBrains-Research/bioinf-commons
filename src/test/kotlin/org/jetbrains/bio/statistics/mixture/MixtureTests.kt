package org.jetbrains.bio.statistics.mixture

import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.viktor.F64Array
import org.junit.Assert.assertArrayEquals
import org.junit.Before

fun assertEqualsUnordered(values1: DoubleArray,
                                  values2: DoubleArray,
                                  eps: Double) {
    val tmp1 = values1.clone()
    val tmp2 = values2.clone()
    tmp1.sort()
    tmp2.sort()
    assertArrayEquals(tmp1, tmp2, eps)
}

/**
 * This is an _unstable_ randomized test. It might fail. Now you KNOW.
 *
 * @author Sergei Lebedev
 * @since 07/11/13
 */
abstract class AbstractMixtureTest {
    protected val weights = F64Array.of(.5, .3, .2)
    protected val numStates = weights.size

    @Before
    fun setUp() {
        Sampling.RANDOM_DATA_GENERATOR.reSeed(42)
    }
}









