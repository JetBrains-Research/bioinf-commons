package org.jetbrains.bio.genome.containers.intersection

import org.jetbrains.bio.util.read
import org.jetbrains.bio.viktor.F64Array
import org.junit.Assert
import org.junit.Test
import java.io.File
import kotlin.math.round
import kotlin.test.assertEquals

class IntersectionInfoTest {
    @Test fun saveToCsv() {
        val info = generateTestData(3, 2)
        val outputPath = File.createTempFile("intersection_df", ".csv").toPath()
        info.saveToCSV(outputPath)
        val outputContent = outputPath.read().replace("\r", "")
        assertEquals(
            """
                |a_names,a_merged_loci_count,col_0,col_1
                |b_merged_loci_count,.,100,101
                |row_0,0,0.0,0.17
                |row_1,1,0.33,0.5
                |row_2,2,0.67,0.83
                |
            """.trimMargin(),
            outputContent)
    }

    @Test fun dumpToString_2x3() {
        assertEquals(
            """
                |a_names,a_merged_loci_count,col_0,col_1,col_2
                |b_merged_loci_count,.,100,101,102
                |row_0,0,0.0,0.17,0.33
                |row_1,1,0.5,0.67,0.83
                |
            """.trimMargin(),
            generateTestData(2, 3).toString()
        )
    }

    @Test fun dumpToString_2x2() {
        Assert.assertEquals(
            """
                |a_names,a_merged_loci_count,col_0,col_1
                |b_merged_loci_count,.,100,101
                |row_0,0,0.0,0.25
                |row_1,1,0.5,0.75
                |
            """.trimMargin(),
            generateTestData(2, 2).toString()
        )
    }

    @Test fun dumpToString_3x2() {
        Assert.assertEquals(
            """
                |a_names,a_merged_loci_count,col_0,col_1
                |b_merged_loci_count,.,100,101
                |row_0,0,0.0,0.17
                |row_1,1,0.33,0.5
                |row_2,2,0.67,0.83
                |
            """.trimMargin(),
            generateTestData(3, 2).toString()
        )
    }

    @Test fun loadFromToCsv() {
        val expectedInfo = generateTestData(3, 2)

        val outputPath = File.createTempFile("intersection_df", ".csv").toPath()
        expectedInfo.saveToCSV(outputPath)
        val loadedInfo = IntersectionInfo.loadFromCSV(outputPath)

        Assert.assertEquals(expectedInfo.toString(), loadedInfo.toString())
        assertEquals(expectedInfo, loadedInfo)
    }

    private fun generateTestData(n: Int, m: Int): IntersectionInfo {
        val nm = n.toDouble() * m
        val data = F64Array(n, m) { i, j ->
            round(100 * (i * m + j) / nm) / 100.0
        }
        val info = IntersectionInfo(
            data,
            (0 until n).map { it }.toIntArray(),
            (0 until m).map { 100 + it }.toIntArray(),
            (0 until n).map { "row_${it}" },
            (0 until m).map { "col_${it}" }
        )
        return info
    }
}