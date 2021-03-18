package org.jetbrains.bio.genome.search

import gnu.trove.list.array.TIntArrayList
import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream

class SimpleSearcher(private val sequence: String) : Searcher {
    override fun find(text: NucleotideSequence): IntStream {
        val string = text.substring(0, text.length)
        val result = TIntArrayList()
        var index = 0
        while (true) {
            index = sequence.indexOf(string, index)
            if (index == -1) {
                break
            }

            result.add(index++)
        }

        return IntStream.range(0, result.size()).map { result[it] }.parallel()
    }
}
