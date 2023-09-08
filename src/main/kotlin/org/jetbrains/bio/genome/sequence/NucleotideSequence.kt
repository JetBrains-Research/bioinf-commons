package org.jetbrains.bio.genome.sequence

import com.google.common.base.Preconditions.checkPositionIndexes
import org.jetbrains.bio.genome.Strand

/**
 * An immutable container for nucleotide sequence.
 *
 * @author Sergei Lebedev
 */
interface NucleotideSequence {
    /**
     * Returns the lowercase nucleotide at specific position.
     */
    fun charAt(pos: Int): Char

    fun nucleotideAt(pos: Int) = Nucleotide.fromChar(charAt(pos))

    /**
     * Returns the lowercase nucleotide at specific position and strand.
     */
    fun charAt(pos: Int, strand: Strand): Char {
        val ch = charAt(pos)
        return if (strand.isPlus()) {
            ch
        } else {
            Nucleotide.complement(ch)
        }
    }

    fun substring(from: Int, to: Int, strand: Strand = Strand.PLUS): String {
        checkPositionIndexes(from, to, length)
        val acc = CharArray(to - from)
        if (strand.isPlus()) {
            for (offset in acc.indices) {
                acc[offset] = charAt(from + offset)
            }
        } else {
            for (offset in acc.indices) {
                acc[acc.size - 1 - offset] = charAt(from + offset, strand)
            }
        }

        return String(acc)
    }

    val length: Int
}

fun CharSequence.asNucleotideSequence(): NucleotideSequence {
    return WrappedCharSequence(this)
}

/**
 * A [NucleotideSequence] adapter for [java.lang.CharSequence].
 *
 * @author Sergei Lebedev
 */
internal data class WrappedCharSequence(val wrapped: CharSequence) : NucleotideSequence {
    override fun charAt(pos: Int) = wrapped[pos].lowercaseChar()

    override val length: Int get() = wrapped.length

    override fun toString() = wrapped.toString()
}