package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.viktor.F64Array
import org.junit.Assert
import org.junit.Test

class BenjaminiHochbergTest {
    @Test
    fun againstR() {
        Assert.assertEquals(
            F64Array.of(0.4, 0.56, 0.80, 0.04),
            BenjaminiHochberg.adjust(F64Array.of(0.2, 0.42, 0.8, 0.01))
        )

        Assert.assertArrayEquals(
            doubleArrayOf(
                0.00949222403622329, 0.00949222403622329, 0.0610324793680683,
                0.078829589901044, 0.0963387314699263, 0.176146462440991,
                0.670529787479448, 0.670529787479448, 0.670529787479448,
                0.955690764788587
            ),
            BenjaminiHochberg.adjust(
                F64Array.of(
                    0.000962882346117542, 0.00189844480724466,
                    0.0183097438104205, 0.0315318359604176,
                    0.0481693657349631, 0.105687877464594,
                    0.543211136961355, 0.565056666152251,
                    0.603476808731503, 0.955690764788587
                )
            ).toDoubleArray(),
            1e-6
        )
    }
}