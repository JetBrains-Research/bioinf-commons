package org.jetbrains.bio.genome

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.assertNull

class RepeatsHs1ParsingTest {

    private val chromosomes = mapOf("chr1" to Chromosome(Genome["to1"], "chr1"))

    @Test
    fun shouldIgnoreFirstLineOfHeader() {
        val header1 = "SW  perc perc perc  query      position in query           matching       repeat              position in  repeat"
        assertNull(RepeatsHs1.parseHs1RepeatsLine(header1, chromosomes))
    }

    @Test
    fun shouldIgnoreSecondLineOfHeader() {
        val header2 = "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID"
        assertNull(RepeatsHs1.parseHs1RepeatsLine(header2, chromosomes))
    }

    @Test
    fun shouldIgnoreEmptyStrings() {
        assertNull(RepeatsHs1.parseHs1RepeatsLine("", chromosomes))
    }

    @Test
    fun shouldParseSimpleRepeat() {
        val simpleRepeat = " 2590   2.9  0.6  1.0  chr1            2    2705 (248384623) +  (ACCCTA)n      Simple_repeat            1 2693    (0) 4166944"
        val result = RepeatsHs1.parseHs1RepeatsLine(simpleRepeat, chromosomes)
        val expectedRepeat = Repeat(name = "(ACCCTA)n", location = Location(2, 2705, chromosomes["chr1"]!!, Strand.PLUS), repeatClass = "Simple_repeat", family = "Simple_repeat")
        assertNotNull(result)
        assertEquals(expectedRepeat, result.second)
    }

    @Test
    fun shouldParseRepeatOnMinusStrand() {
        val minus = " 1134  25.5  9.2  7.6  chr1         4083    4660 (248382668) C  LTR60B         LTR/ERV1               (0)  765    178 4166946"
        val result = RepeatsHs1.parseHs1RepeatsLine(minus, chromosomes)
        val expectedRepeat = Repeat(name = "LTR60B", location = Location(4083, 4660, chromosomes["chr1"]!!, Strand.MINUS), repeatClass = "ltr", family = "erv1")
        assertNotNull(result)
        assertEquals(expectedRepeat, result.second)
    }

    @Test
    fun shouldParseRepeatOnPlusStrand() {
        val plus = " 1293  21.3  7.2  1.9  chr1         4664    5263 (248382065) +  L1MC3          LINE/L1               4913 5546 (2239) 4166947"
        val result = RepeatsHs1.parseHs1RepeatsLine(plus, chromosomes)
        val expectedRepeat = Repeat(name = "L1MC3", location = Location(4664, 5263, chromosomes["chr1"]!!, Strand.PLUS), repeatClass = "line", family = "l1")
        assertNotNull(result)
        assertEquals(expectedRepeat, result.second)
    }

}