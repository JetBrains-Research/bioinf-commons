package org.jetbrains.bio.genome.search

import org.jetbrains.bio.genome.sequence.NucleotideSequence
import java.util.stream.IntStream

interface Searcher {
    fun find(text: NucleotideSequence): IntStream
}
