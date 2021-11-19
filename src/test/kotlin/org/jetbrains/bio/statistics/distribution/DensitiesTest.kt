package org.jetbrains.bio.statistics.distribution

import junit.framework.TestCase.assertEquals
import org.apache.commons.math3.util.FastMath
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.test.assertTrue

class DensitiesTest {
    @Test
    fun testLogDirichletDensity() {
        // Validated against 'ddensity' from 'gtools' CRAN package.
        assertEquals(
            0.08106397,
            FastMath.exp(
                Densities.logDirichletDensity(
                    F64Array.of(.1, .6, .3),
                    F64Array.of(.5, 42.0, 13.0)
                )
            ),
            1e-8
        )
        assertEquals(
            2.569954,
            FastMath.exp(
                Densities.logDirichletDensity(
                    F64Array.of(0.2059613, 0.2574219, 0.5366168),
                    F64Array.of(7.0, 22.0, 27.0)
                )
            ),
            1e-6
        )
    }

    @Test
    fun testOverflow() {
        assertTrue(
            Densities.logPoissonDensity(1000, 10.0).isFinite(),
            "log Poisson density overflowed"
        )
    }
}
