package org.jetbrains.bio.query

import org.junit.Test
import kotlin.test.assertEquals

class IdReducerTest {
    @Test
    fun reducedIdSingle() {
        assertEquals("GSM42", reduceIds(listOf("GSM42")))
        assertEquals("GSM42", reduceIds(listOf("GSM42_something_else"), 5))
        assertEquals("non_gsm", reduceIds(listOf("non_gsm")))
        assertEquals("YD_ac_BC1_R1_hg19.bam", reduceIds(listOf("YD_ac_BC1_R1_hg19.bam")))
    }

    @Test
    fun reducedIdMultipleNoCommon() {
        assertEquals("foo_bar", reduceIds(listOf("foo", "bar")))
    }

    @Test
    fun reducedIdMultipleCommonPrefix() {
        assertEquals("GSM42_foo_GSM43_bar", reduceIds(listOf("GSM42_foo", "GSM43_bar")))
        assertEquals("GSM42", reduceIds(listOf("GSM42_foo", "GSM42_bar"), 5))
        assertEquals("GSM42", reduceIds(listOf("foo_GSM42", "bar_GSM42"), 5))
        assertEquals("GSM423_GSM424", reduceIds(listOf("GSM423", "GSM424")))
        assertEquals("foo_bar", reduceIds(listOf("foo", "bar")))
    }

    @Test
    fun reducedIdMultipleCommonSuffix() {
        assertEquals("foo_boo_bar_boo", reduceIds(listOf("foo_boo", "bar_boo")))
        assertEquals("boo_foo_bar", reduceIds(listOf("foo_boo", "bar_boo"), 11))
    }

    @Test
    fun reduceInclusion() {
        assertEquals("span_track1584577470751151686_unique_200",
                reduceIds(listOf("span_track1584577470751151686_unique#sha1.bam", "200")))
    }


    @Test
    fun realNames() {
        val ids = listOf(
                "GSM1102782_UW.CD14_Primary_Cells.H3K27ac.RO_01721.Histone.DS22926.bed.gz",
                "GSM772872_BI.CD8_Naive_Primary_Cells.H3K36me3.Donor_100_7_pooled_leukopaks_Jan_7_2011.bed.gz",
                "GSM613881_UCSF-UBC.CD8_Naive_Primary_Cells.H3K4me3.TC010.bed.gz")
        assertEquals("GSM1102782_GSM772872_GSM613881", reduceIds(ids))
    }

    @Test
    fun realNamesMethylome() {
        val ids = listOf(
                "GSM4326_Bisulfite-Seq.methylC-seq_87_IMR90_r1a",
                "Bisulfite-Seq.methylC-seq_88_IMR90_r1b",
                "Bisulfite-Seq.methylC-seq_89_IMR90_r1c",
                "Bisulfite-Seq.methylC-seq_90_IMR90_r2a",
                "Bisulfite-Seq.methylC-seq_91_IMR90_r2b",
                "Bisulfite-Seq.methylC-seq_92_IMR90_r2c")
        assertEquals("methylC-seq_IMR90_Bisulfite-Seq_GSM4326_88_r1b_89_r1c_90_r2a_91_r2b_92_r2c",
                reduceIds(ids))
    }

    @Test
    fun methylomes10vs10() {
        val path = "/mnt/stripe/bio/experiments/cache/fit/vbdpcsmbhmm/methylationenrichment_hg19_3_".length
        val ids = listOf(
                "OD1.sorted#00000",
                "OD10.sorted#11111",
                "OD2.sorted#22222",
                "OD3.sorted#33333",
                "OD4.sorted#44444",
                "YD1.sorted#55555",
                "YD10.sorted#66666",
                "YD2.sorted#77777",
                "YD3.sorted#88888")
        assertEquals("sorted_OD1_OD10_OD2_OD3_OD4_YD1_YD10_YD2_YD3", reduceIds(ids, maxLength = 150 - path))
        assertEquals("sorted_OD1_OD10#f6361", reduceIds(ids, maxLength = 100 - path))
    }

    @Test
    fun sha() {
        val ids = listOf("0123456789012345678901234567890")
        assertEquals("01234567890123#a4a44", reduceIds(ids, maxLength = 20))
        assertEquals("0123456789#a4a44", reduceIds(ids, maxLength = 16))
    }
}