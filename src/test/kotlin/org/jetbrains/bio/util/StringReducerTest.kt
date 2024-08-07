package org.jetbrains.bio.util

import org.junit.Test
import kotlin.test.assertEquals

class StringReducerTest {
    @Test
    fun reducedIdSingle() {
        assertEquals("GSM42", reduceIds(listOf("GSM42")))
        assertEquals("GSM42", reduceIds(listOf("GSM42_something_else"), 12))
        assertEquals("non_gsm", reduceIds(listOf("non_gsm")))
        assertEquals("YD_ac_BC1_R1_hg19_bam", reduceIds(listOf("YD_ac_BC1_R1_hg19.bam")))
    }

    @Test
    fun realNames() {
        val ids = listOf(
            "GSM1102782_UW.CD14_Primary_Cells.H3K27ac.RO_01721.Histone.DS22926.bed.gz",
            "GSM772872_BI.CD8_Naive_Primary_Cells.H3K36me3.Donor_100_7_pooled_leukopaks_Jan_7_2011.bed.gz",
            "GSM613881_UCSF-UBC.CD8_Naive_Primary_Cells.H3K4me3.TC010.bed.gz"
        )
        assertEquals("GSM1102782_GSM772872_GSM613881", reduceIds(ids))
    }

    @Test
    fun methylomes10vs10() {
        val path = "/mnt/stripe/bio/experiments/cache/fit/vbdpcsmbhmm/methylationenrichment_hg19_3_".length
        assertEquals(
            "OD1_sorted_OD10_sorted_OD2_sorted_OD3_sorted_OD4_sorted_YD1_sorte#f6361",

            reduceIds(
                listOf(
                    "OD1.sorted#00000",
                    "OD10.sorted#11111",
                    "OD2.sorted#22222",
                    "OD3.sorted#33333",
                    "OD4.sorted#44444",
                    "YD1.sorted#55555",
                    "YD10.sorted#66666",
                    "YD2.sorted#77777",
                    "YD3.sorted#88888"
                ), maxLength = 150 - path
            )
        )
        assertEquals(
            "OD1_sorted_OD10_sorted_OD2_sorted_OD3_sorted_OD4_sorted_YD1_sorted",
            reduceIds(
                listOf(
                    "OD1.sorted#00000",
                    "OD10.sorted#11111",
                    "OD2.sorted#22222",
                    "OD3.sorted#33333",
                    "OD4.sorted#44444",
                    "YD1.sorted#55555"
                ), maxLength = 150 - path
            )
        )
    }

    @Test
    fun sha() {
        val ids = listOf("0123456789012345678901234567890")
        assertEquals("01234567890123#a4a44", reduceIds(ids, maxLength = 20))
        assertEquals("0123456789#a4a44", reduceIds(ids, maxLength = 16))
    }
}