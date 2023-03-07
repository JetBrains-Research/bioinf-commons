package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.statistics.hypothesis.StofferLiptakTest.Companion.NORMAL
import org.junit.Assert
import org.junit.Test

class StofferLiptakTestTest {

    @Test
    fun testNormalInverseCumulative() {
        Assert.assertTrue(NORMAL.inverseCumulativeProbability(1.0 - StofferLiptakTest.EPSILON).isFinite())
        Assert.assertFalse(NORMAL.inverseCumulativeProbability(1.0 - StofferLiptakTest.EPSILON * 0.5).isFinite())
    }

    @Test
    fun testCorrelations() {
        val correlations = StofferLiptakTest.computeCorrelations(doubleArrayOf(0.0, 0.1, 0.2, 0.3, 0.4), 20)
        Assert.assertEquals(3, correlations.size)
        Assert.assertTrue(doubleArrayOf(0.0, 1.0, 1.0) contentEquals correlations)
    }

    @Test
    fun testZeroesCorrelations() {
        val correlations = StofferLiptakTest.computeCorrelations(doubleArrayOf(0.0, 0.0, 0.0, 0.0), 2)
        Assert.assertEquals(3, correlations.size)
        Assert.assertTrue(doubleArrayOf(0.0, 0.0, 0.0) contentEquals correlations)
    }

    @Test
    fun testZScore() {
        val stofferLiptakTest = StofferLiptakTest(doubleArrayOf(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.2, 0.3))
        Assert.assertEquals(8.209536151601387, stofferLiptakTest.zscore(0.0), 1e-6)
        Assert.assertEquals(-8.209536151601387, stofferLiptakTest.zscore(1.0), 1e-6)
        Assert.assertEquals(1.6448536269514724, stofferLiptakTest.zscore(0.05), 1e-6)
    }

    @Test
    fun testCombine() {
        val stofferLiptakTest = StofferLiptakTest(doubleArrayOf(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.2, 0.3))
        // Check that results are monotonous by pvalues and size
        Assert.assertEquals(1e-6, stofferLiptakTest.combine(doubleArrayOf(1e-6)), 1e-6)
        Assert.assertEquals(0.004272487313437656, stofferLiptakTest.combine(doubleArrayOf(1e-6, 1e-4)), 1e-6)
        Assert.assertEquals(0.1824166413929258, stofferLiptakTest.combine(doubleArrayOf(1e-2, 0.1)), 1e-6)
        // Improving relaxed pvalue 1e-4
        Assert.assertEquals(5.961866446146935E-5, stofferLiptakTest.combine(doubleArrayOf(1e-8, 1e-6, 1e-4)), 1e-6)
        Assert.assertEquals(
            2.351597618743817E-6, stofferLiptakTest.combine(doubleArrayOf(1e-10, 1e-8, 1e-6, 1e-4)), 1e-6
        )
        Assert.assertEquals(
            1.6964616000869626E-7, stofferLiptakTest.combine(doubleArrayOf(1e-10, 1e-8, 1e-10, 1e-4)), 1e-12
        )
        // Case when combining is much worse
        Assert.assertEquals(
            0.010885145722788425,
            stofferLiptakTest.combine(doubleArrayOf(9.61429557064954E-8, 5.878646086208946E-4)), 1e-6
        )
    }

    @Test
    fun testNegativePCorrelations() {
        // See https://github.com/JetBrains-Research/span/issues/32
        val pValues = doubleArrayOf(
            4.536832820715596E-24, 1.8393263841715094E-27, 5.994971392624784E-21, 2.06835737356540E-42,
            1.8393263841715094E-27, 6.085984087566133E-11, 7.803572981068073E-13, 3.0593896259456266E-34,
            5.994971392624784E-21, 4.2782295607092716E-33, 2.06835737356540E-42, 1.3611161918897815E-65,
            4.2782295607092716E-33, 2.068357373565406E-42, 1.504096577367313E-36, 6.066738108307561E-77,
            3.2379502029748828E-78, 3.62203538305709E-24, 5.534060003552777E-46, 7.022847588883002E-12,
            8.022790630864596E-31, 6.003874037766238E-21, 0.007552237026621174, 0.007552237026621174,
            8.034673975762777E-31, 5.994971392624784E-21, 2.2119764512649164E-7, 2.9106729037243954E-26,
            1.0831485550464778E-15, 3.2206438254678727E-4, 2.942794634147554E-25, 8.117383475641987E-18,
            1.686895628392727E-15, 6.085984087566133E-11, 1.8393263841715094E-27, 4.034142241199322E-9,
            4.034142241199322E-9, 2.1589050666765635E-35, 3.622035383057091E-24, 8.64047594152862E-45,
            3.2206438254678727E-4, 4.731380466089632E-40, 5.176489595662274E-22, 0.0016398706491104705,
            6.794465920204226E-20, 8.370554457775231E-14, 3.06498754606286E-8, 5.994971392624784E-21,
            4.251326388864969E-68, 1.2167142708556864E-93, 1.2948467168860364E-60, 7.755118251507132E-124,
            2.2698495981315425E-105, 3.62235383057091E-24, 4.3732795626831537E-23, 9.090850871918053E-81,
            6.536893123315908E-57, 3.0593896259456266E-34, 2.0979958890372006E-74, 3.8375216365123104E-58,
            6.06738108307561E-77, 5.412172008363167E-108, 7.755118251507132E-124, 5.105010789765497E-52,
            3.4976834298945663E-47, 3.4976834298945663E-47, 3.622035383057091E-24, 5.90734275309314E-32,
            6.536893123315908E-57, 4.2524239202042556E-63, 6.792263176744292E-20, 8.539658036579827E-17,
            3.622035383057091E-24, 8.117383475641987E-18, 9.25062722468583E-90, 8.022790630864596E-31,
            3.089724755968893E-53, 4.722908695981398E-91, 1.7324958657514746E-130, 2.2368774805061296E-59,
            8.539658036579827E-17, 1.773316710600536E-28
        )
        val stofferLiptakTest = StofferLiptakTest(
            pValues
        )
        // Check that result is not Nan
        Assert.assertFalse(stofferLiptakTest.combine(pValues).isNaN())
    }


}