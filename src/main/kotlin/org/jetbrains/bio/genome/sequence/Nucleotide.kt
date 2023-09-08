package org.jetbrains.bio.genome.sequence


import com.google.common.base.Preconditions.checkElementIndex

/**
 * DNA nucleotide.
 */
enum class Nucleotide constructor(
    /**
     * Two bit coding as defined by the 2bit format.
     *
     * See http://genome.ucsc.edu/FAQ/FAQformat.html#format7.
     */
    val byte: Byte
) {
    // Please keep 'em in bytes order
    T(0b00),
    C(0b01),
    A(0b10),
    G(0b11);

    val char: Char get() = getChar(byte)

    init {
        // is important for Nucleotide.fromChar method
        require(byte.toInt() == this.ordinal)
    }

    companion object {
        /**
         * Lower case letters
         */
        val ALPHABET = CharArray(values().size)

        init {
            for (nucleotide in values()) {
                ALPHABET[nucleotide.byte.toInt()] = nucleotide.toString().first().lowercaseChar()
            }
        }

        internal const val ANY_NUCLEOTIDE_BYTE: Byte = 0b101
        const val N = 'n'

        fun fromChar(c: Char): Nucleotide? {
            val byte = getByte(c)
            if (byte == ANY_NUCLEOTIDE_BYTE) {
                return null
            }
            return Nucleotide.values()[byte.toInt()]
        }

        internal fun getByte(c: Char) = when (c) {
            'T', 't' -> T.byte
            'C', 'c' -> C.byte
            'A', 'a' -> A.byte
            'G', 'g' -> G.byte
            else -> ANY_NUCLEOTIDE_BYTE
        }

        internal fun getChar(b: Byte) = when (b) {
            ANY_NUCLEOTIDE_BYTE -> N
            else -> {
                checkElementIndex(b.toInt(), 4, "invalid nucleotide byte")
                ALPHABET[b.toInt()]
            }
        }

        fun complement(ch: Char) = when (ch) {
            'a', 'A' -> 't'
            't', 'T' -> 'a'
            'c', 'C' -> 'g'
            'g', 'G' -> 'c'
            'n' -> 'n'
            else -> error("invalid nucleotide: $ch")
        }
    }
}
