package org.jetbrains.bio.genome.search.sa

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.search.Searcher
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream
import kotlin.math.min

/**
 * @author Oleg Shpynov
 */
class SuffixArraySearcher(val sa: SuffixArray) : Searcher {

    constructor(chromosome: Chromosome) : this(SuffixArray.load(chromosome))

    override fun find(text: NucleotideSequence): IntStream {
        var startOffset = 0
        // Important: SA length may be technically bigger than length of sequence
        var endOffset = sa.sequence.length - 1
        startOffset = searchMatchAtLeftOrRightBound(sa, text, startOffset, endOffset, true)
        endOffset = searchMatchAtLeftOrRightBound(sa, text, startOffset, endOffset, false)
        return if (startOffset <= endOffset) {
            IntStream.rangeClosed(startOffset, endOffset).map { sa.index(it) }
        } else
            IntStream.empty()
    }

    companion object {

        fun searchMatchAtLeftOrRightBound(
            suffixArray: SuffixArray,
            pattern: NucleotideSequence,
            startOffset: Int,
            endOffset: Int,
            searchLeft: Boolean
        ): Int {
            val sequence = suffixArray.sequence
            val patternLength = pattern.length
            val textLength = sequence.length

            // matched symbols quantity
            var matchedSymbolsAtLeftBound = 0
            var matchedSymbolsAtRightBound = 0

            var left = startOffset
            var right = endOffset
            while (left <= right) {
                val m = (left + right + 1) / 2

                // guaranteed matched symbols number in current text[start, end] region
                // on each iteration we choosing left:text[start, median] or right[median, end] part
                // let's track amount of matched symbols in left/right parts separately, thus on each
                // iteration we can guarantee that Min(matchedSymbolsAtLeftBound, matchedSymbolsAtRightBound)
                // symbols are matched in whole text[start, end] region
                val patternMatchedSymbols = min(matchedSymbolsAtLeftBound, matchedSymbolsAtRightBound)

                // let's continue from guaranteed matched symbol
                var posInPattern = patternMatchedSymbols

                // current position in text = suffix start offset + guaranteed matched symbols
                val suffixStartPosInText = suffixArray.index(m)
                var posInText = suffixStartPosInText + patternMatchedSymbols

                // if pattern matches text suffix - match as long prefix as possible
                while (0 <= posInText && posInPattern < patternLength && posInText < textLength &&
                    sequence.charAt(posInText) == pattern.charAt(posInPattern)
                ) {
                    posInPattern++
                    posInText++
                }
                // * whole pattern already matched ---> find start position of suffixes which matches pattern
                if (posInPattern >= patternLength) {
                    // 2 policies: look up/down
                    if (searchLeft) {
                        // * look left ("up") -->  [left, m]
                        matchedSymbolsAtRightBound = posInPattern
                        right = m - 1
                    } else {
                        // * look right ("down") --> [m, right]
                        matchedSymbolsAtLeftBound = posInPattern
                        left = m + 1
                    }
                } else {
                    // suffix is shorter - look to next one or due to mismatch type
                    if (posInText < 0 || posInText >= textLength ||
                        sequence.charAt(posInText) < pattern.charAt(posInPattern)
                    ) {
                        // [m, right]
                        matchedSymbolsAtLeftBound = posInPattern
                        left = m + 1
                    } else {
                        // [left, m]
                        matchedSymbolsAtRightBound = posInPattern
                        right = m - 1
                    }
                }
            }
            return if (searchLeft) right + 1 else left - 1
        }
    }
}
