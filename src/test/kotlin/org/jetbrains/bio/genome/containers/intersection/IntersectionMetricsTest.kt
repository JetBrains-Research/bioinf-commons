package org.jetbrains.bio.genome.containers.intersection

import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.LocationsSortedList
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.toQuery
import org.junit.Test
import java.io.StringReader
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class IntersectionMetricsTest {
    @Test
    fun overlapTest_NoIntersection() {
        doCheckOverlapAB(
            """
                |chr1,5,10
            """.trimMargin(),
            """
                |chr1,10,12
            """.trimMargin(),
            0, 0.0
        )
    }

    @Test
    fun overlapTest_NoIntersection2() {
        doCheckOverlapAB(
            """
                |chr1,9,10
            """.trimMargin(),
            """
                |chr1,8,9
            """.trimMargin(),
            0, 0.0
        )
    }

    @Test
    fun overlapTest_1bpSingle() {
        doCheckOverlapAB(
            """
                |chr1,5,10
            """.trimMargin(),
            """
                |chr1,9,12
            """.trimMargin(),
            1, 1.0
        )
    }

    @Test
    fun overlapTest_1bpSingle2() {
        doCheckOverlapAB(
            """
                |chr1,9,10
            """.trimMargin(),
            """
                |chr1,9,10
            """.trimMargin(),
            1, 1.0
        )
    }

    @Test
    fun overlapTest_1bpMultiple_Merged() {
        doCheckOverlapAB(
            """
                |chr1,5,10
                |chr1,6,11
                |chr1,20,21
            """.trimMargin(),
            """
                |chr1,9,12
                |chr1,20,21
            """.trimMargin(),
            2, 1.0, "merged"
        )
    }

    @Test
    fun overlapTest_1bpMultiple_Sorted() {
        doCheckOverlapAB(
            """
                |chr1,5,10
                |chr1,6,11
                |chr1,20,21
            """.trimMargin(),
            """
                |chr1,9,12
                |chr1,20,21
            """.trimMargin(),
            3, 1.0, "sorted"
        )
    }

    @Test
    fun overlapTest_Other_Merged() {
        doCheckOverlapAB(
            """
                |chr1,5,15
                |chr1,6,11
                |chr1,18,21
                |chr1,22,23
                |chr1,23,24
                |chr1,26,28
            """.trimMargin(),
            """
                |chr1,6,7
                |chr1,9,12
                |chr1,20,26
                |chr1,30,36
            """.trimMargin(),
            3, 0.75, "merged"
        )
    }

    @Test
    fun overlapTest_Other_Sorted() {
        doCheckOverlapAB(
            """
                |chr1,5,15
                |chr1,6,11
                |chr1,18,21
                |chr1,22,23
                |chr1,23,24
                |chr1,26,28
            """.trimMargin(),
            """
                |chr1,6,7
                |chr1,9,12
                |chr1,20,26
                |chr1,30,36
            """.trimMargin(),
            5, 0.8333, "sorted"
        )
    }

    @Test
    fun jaccardTest_NoIntersection() {
        doCheckJaccardAB(
            """
                |chr1,5,10
            """.trimMargin(),
            """
                |chr1,10,12
            """.trimMargin(),
            0.0
        )
    }

    @Test
    fun jaccardTest_NoIntersection2() {
        doCheckJaccardAB(
            """
                |chr1,9,10
            """.trimMargin(),
            """
                |chr1,8,9
            """.trimMargin(),
            0.0
        )
    }

    @Test
    fun jaccardTest_1bpSingle() {
        doCheckJaccardAB(
            """
                |chr1,5,10
            """.trimMargin(),
            """
                |chr1,9,12
            """.trimMargin(),
            1.0 * 1 / (5 + 3 - 1)
        )
    }

    @Test
    fun jaccardTest_1bpSingle2() {
        doCheckJaccardAB(
            """
                |chr1,9,10
            """.trimMargin(),
            """
                |chr1,9,10
            """.trimMargin(),
            1.0
        )
    }

    @Test
    fun jaccardTest_1bpMultiple_Merged() {
        doCheckJaccardAB(
            """
                |chr1,5,10
                |chr1,6,11
                |chr1,20,21
            """.trimMargin(),
            """ 
                |chr1,9,12
                |chr1,20,21
            """.trimMargin(),
            // {[5,11), [20, 21)}, {[9,12), [20, 21)}
            1.0 * (2 + 1) / (6 + 1 + 3 + 1 - (2 + 1)),
            "merged"
        )
    }

    @Test
    fun jaccardTest_1bpMultiple_Sorted() {
        doCheckJaccardAB(
            """
                   |chr1,5,10
                   |chr1,6,11
                   |chr1,20,21
               """.trimMargin(),
            """ 
                   |chr1,9,12
                   |chr1,20,21
               """.trimMargin(),
            // {[9,10)}, {[9,11)}, {[20,21)},
            1.0 * (1 + 2 + 1) / (5 + 5 + 1 + 3 + 1 - (1 + 2 + 1)),
            "sorted"
        )
    }

    @Test
    fun jaccardTest_Other_Merged() {
        doCheckJaccardAB(
            """
                |chr1,5,15
                |chr1,6,11
                |chr1,18,21
                |chr1,22,23
                |chr1,23,24
                |chr1,26,28
            """.trimMargin(),
            """
                |chr1,6,7
                |chr1,9,12
                |chr1,20,26
                |chr1,30,36
            """.trimMargin(),
            // {[5,15), [18, 21), [22,24), [26,28)}, {[6, 7), [9,12), [20, 26), [30, 36)}
            1.0 * (1 + 3 + 1 + 2) / ((5 + 5 + 3 + 2 + 2) + (1 + 3 + 6 + 6) - (1 + 3 + 1 + 2)),
            "merged"
        )
    }

    @Test
    fun jaccardTest_Other_Sorted() {
        doCheckJaccardAB(
            """
                |chr1,5,15
                |chr1,6,11
                |chr1,18,21
                |chr1,22,23
                |chr1,23,24
                |chr1,26,28
            """.trimMargin(),
            """
                |chr1,6,7
                |chr1,9,12
                |chr1,20,26
                |chr1,30,36
            """.trimMargin(),
            // {[6,7), [9,12)}, {[6,7), [9,11)}, {[20,21)}, {[22,23)}, {[23,24)}, {}}
            1.0 * (1 + 3 + 1 + 2 + 1 + 1 + 1) / (10 + 5 + 3 + 1 + 1 + 2 + 1 + 3 + 6 + 6 - (1 + 3 + 1 + 2 + 1 + 1 + 1)),
            "sorted"
        )
    }


    private fun doCheckOverlapAB(
        aContent: String, bContent: String, overlapNum: Long, overlapFract: Double, mode: String = "both"
    ) {
        val gq = Genome["to1"].toQuery()
        val bedFormat = BedFormat.from("bed3", ',')

        require(mode in listOf("both", "sorted", "merged"))
        if (mode == "both" || mode == "merged") {
            val a = LocationsMergingList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsMergingList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(overlapNum, IntersectionMetrics.OVERLAP.calcMetric(a, b).toLong())
            assertEquals(overlapFract, IntersectionMetrics.OVERLAP_FRACTION.calcMetric(a, b))
        }

        if (mode == "both" || mode == "sorted") {
            val a = LocationsSortedList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsSortedList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(overlapNum, IntersectionMetrics.OVERLAP.calcMetric(a, b).toLong())
            assertEquals(overlapFract, IntersectionMetrics.OVERLAP_FRACTION.calcMetric(a, b))
        }
    }

    @Test
    fun intersectionNumberTest_NoIntersection() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,10
              """.trimMargin(),
            """
                  |chr1,10,12
              """.trimMargin(),
            0
        )
    }

    @Test
    fun intersectionNumberTest_NoIntersection2() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,9,10
              """.trimMargin(),
            """
                  |chr1,8,9
              """.trimMargin(),
            0
        )
    }

    @Test
    fun intersectionNumberTest_1bpSingle() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,10
              """.trimMargin(),
            """
                  |chr1,9,12
              """.trimMargin(),
            1
        )
    }

    @Test
    fun intersectionNumberTest_1bpSingle2() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,9,10
              """.trimMargin(),
            """
                  |chr1,9,10
              """.trimMargin(),
            1
        )
    }

    @Test
    fun intersectionNumberTest_1bpMultiple_Merged() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,10
                  |chr1,6,11
                  |chr1,20,21
              """.trimMargin(),
            """
                  |chr1,9,12
                  |chr1,20,21
              """.trimMargin(),
            2, "merged"
        )
    }

    @Test
    fun intersectionNumberTest_1bpMultiple_Sorted() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,10
                  |chr1,6,11
                  |chr1,20,21
              """.trimMargin(),
            """
                  |chr1,9,12
                  |chr1,20,21
              """.trimMargin(),
            3, "sorted"
        )
    }

    @Test
    fun intersectionNumberTest_Other_Merged() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,15
                  |chr1,6,11
                  |chr1,18,21
                  |chr1,22,23
                  |chr1,23,24
                  |chr1,26,28
              """.trimMargin(),
            """
                  |chr1,6,7
                  |chr1,9,12
                  |chr1,20,26
                  |chr1,30,36
              """.trimMargin(),
            4, "merged"
        )
    }

    @Test
    fun intersectionNumberTest_Other_Sorted() {
        doCheckIntersectionNumberAB(
            """
                  |chr1,5,15
                  |chr1,6,11
                  |chr1,18,21
                  |chr1,22,23
                  |chr1,23,24
                  |chr1,26,28
              """.trimMargin(),
            """
                  |chr1,6,7
                  |chr1,9,12
                  |chr1,20,26
                  |chr1,30,36
              """.trimMargin(),
            7, "sorted"
        )
    }

    private fun doCheckIntersectionNumberAB(
        aContent: String, bContent: String, intersectionNum: Long, mode: String = "both"
    ) {
        val gq = Genome["to1"].toQuery()
        val bedFormat = BedFormat.from("bed3", ',')

        require(mode in listOf("both", "sorted", "merged"))
        if (mode == "both" || mode == "merged") {
            val a = LocationsMergingList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsMergingList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(intersectionNum, IntersectionMetrics.INTERSECTION_NUMBER.calcMetric(a, b).toLong())
        }

        if (mode == "both" || mode == "sorted") {
            val a = LocationsSortedList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsSortedList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(intersectionNum, IntersectionMetrics.INTERSECTION_NUMBER.calcMetric(a, b).toLong())
        }
    }

    private fun doCheckJaccardAB(
        aContent: String, bContent: String, jaccardIndex: Double, mode: String = "both"
    ) {
        val gq = Genome["to1"].toQuery()
        val bedFormat = BedFormat.from("bed3", ',')

        require(mode in listOf("both", "sorted", "merged"))
        if (mode == "both" || mode == "merged") {
            val a = LocationsMergingList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsMergingList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(jaccardIndex, IntersectionMetrics.JACCARD.calcMetric(a, b))
            assertEquals(jaccardIndex, IntersectionMetrics.JACCARD.calcMetric(b, a))
        }

        if (mode == "both" || mode == "sorted") {
            val a = LocationsSortedList.load(gq, StringReader(aContent), "a.csv", bedFormat)
            val b = LocationsSortedList.load(gq, StringReader(bContent), "b.csv", bedFormat)
            assertEquals(jaccardIndex, IntersectionMetrics.JACCARD.calcMetric(a, b))
            assertEquals(jaccardIndex, IntersectionMetrics.JACCARD.calcMetric(b, a))
        }
    }

    fun assertEquals(expected: Double, actual: Double, eps: Double = 0.0001) {
        assertTrue(
            Precision.equals(expected, actual, eps),
            "Expected: $expected, but was $actual (precision: $eps)"
        )
    }
}
