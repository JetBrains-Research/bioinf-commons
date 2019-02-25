package org.jetbrains.bio.util

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * @author Oleg Shpynov
 * @since 01/04/16.
 */
class ParseTest {

    @Test fun testTokens() {
        val tokenizer = Tokenizer(" { 1 , 10 , 1 } ", setOf(Lexeme("{"), Lexeme(","), Lexeme("}")))
        val lexemes = arrayListOf<String>()
        while (!tokenizer.atEnd()) {
            lexemes.add(tokenizer.fetch().toString())
            tokenizer.next()
        }
        assertEquals("{;1;,;10;,;1;}", lexemes.joinToString(";"))
    }

    @Test fun testLookahead() {
        val tokenizer = Tokenizer("{ 1, 10, 1 }", setOf(Lexeme("{"), Lexeme(","), Lexeme("}")))
        tokenizer.next()
        tokenizer.lookahead(Match(Lexeme("XXX"), 1, 3))
        assertEquals("XXX", tokenizer.fetch()!!.token)
        tokenizer.next()
        assertEquals(",", tokenizer.fetch()!!.token)
    }


    @Test fun testParseInt() {
        val tokenizer = Tokenizer("123456", emptySet())
        assertEquals(123456, tokenizer.parseInt())
        assertTrue(tokenizer.atEnd())
    }

    @Test fun testParseIntWhitespaces() {
        val tokenizer = Tokenizer("\t123456 ", emptySet())
        assertEquals(123456, tokenizer.parseInt())
        assertTrue(tokenizer.atEnd())
    }

    @Test fun testParseDouble() {
        val tokenizer = Tokenizer("123456.5", emptySet())
        assertEquals(123456.5, tokenizer.parseDouble())
    }

    @Test fun testRegex() {
        val delimiter = RegexLexeme("[_;,]|\\.\\.")
        val tokenizer = Tokenizer("1 _ 2 ; 3 , 4 .. 5", setOf(delimiter))
        assertEquals(1, tokenizer.parseInt())
        val list = mutableListOf(1)
        while (tokenizer.fetch() == delimiter) {
            tokenizer.next()
            list.add(tokenizer.parseInt())
        }
        assertTrue(tokenizer.atEnd())
        assertEquals("[1, 2, 3, 4, 5]", list.toString())
    }

}