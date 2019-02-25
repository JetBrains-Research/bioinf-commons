package org.jetbrains.bio.genome.sequence

import java.util.*
import java.util.stream.Collectors

fun Random.nextString(alphabet: String, length: Int): String {
    return ints(length.toLong(), 0, alphabet.length)
            .mapToObj { alphabet[it].toString() }
            .collect(Collectors.joining())
}