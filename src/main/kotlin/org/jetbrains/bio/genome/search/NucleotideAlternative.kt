package org.jetbrains.bio.genome.sequence

/**
 * Extended DNA / RNA alphabet.
 *
 * Symbol Meaning           Nucleic Acid
 * A     A                  Adenine
 * C     C                  Cytosine
 * G     G                  Guanine
 * T     T                  Thymine
 * U     U                  Uracil
 * M     A or C             aMino
 * R     A or G             puRine
 * W     A or T             Weak
 * S     C or G             Strong
 * Y     C or T             pYrimidine
 * K     G or T             Keto
 * V     A or C or G        not T (V)
 * H     A or C or T        not G (H)
 * D     A or G or T        not C (D)
 * B     C or G or T        not A (B)
 * X     G or A or T or C   any (not recommended)
 * N     G or A or T or C   aNy
 */
enum class NucleotideAlternative(vararg nucleotides: Nucleotide) {
    U(Nucleotide.T),
    A(Nucleotide.A),
    C(Nucleotide.C),
    G(Nucleotide.G),
    T(Nucleotide.T),

    M(Nucleotide.A, Nucleotide.C),
    R(Nucleotide.A, Nucleotide.G),
    W(Nucleotide.A, Nucleotide.T),
    S(Nucleotide.C, Nucleotide.G),
    Y(Nucleotide.C, Nucleotide.T),
    K(Nucleotide.G, Nucleotide.T),

    V(Nucleotide.A, Nucleotide.C, Nucleotide.G),
    H(Nucleotide.A, Nucleotide.C, Nucleotide.T),
    D(Nucleotide.A, Nucleotide.T, Nucleotide.G),
    B(Nucleotide.C, Nucleotide.T, Nucleotide.G),
    N(*Nucleotide.values());

    val alternatives: ByteArray

    init {
        alternatives = ByteArray(nucleotides.size)
        for (i in nucleotides.indices) {
            alternatives[i] = nucleotides[i].byte
        }
    }

    val alternativesCount: Int get() = alternatives.size

    fun match(b: Byte): Boolean {
        for (alternative in alternatives) {
            if (alternative == b) {
                return true
            }
        }

        return false
    }

    fun match(c: Char): Boolean {
        for (alternative in alternatives) {
            if (alternative == Nucleotide.getByte(c)) {
                return true
            }
        }

        return false
    }

    companion object {
        @JvmStatic
        fun fromChar(c: Char) = when (c) {
            'a', 'A' -> A
            't', 'T' -> T
            'c', 'C' -> C
            'g', 'G' -> G
            'u', 'U' -> U
            'm', 'M' -> M
            'r', 'R' -> R
            'w', 'W' -> W
            's', 'S' -> S
            'y', 'Y' -> Y
            'k', 'K' -> K
            'v', 'V' -> V
            'h', 'H' -> H
            'd', 'D' -> D
            'b', 'B' -> B
            else -> N
        }
    }
}
