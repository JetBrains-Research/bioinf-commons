package org.jetbrains.bio.query

import com.google.common.collect.Maps
import com.google.common.collect.Multimap
import com.google.common.collect.Multimaps
import com.google.common.collect.Sets.newHashSet
import gnu.trove.map.hash.TObjectIntHashMap
import org.jetbrains.bio.util.HASH_PREFIX
import org.jetbrains.bio.util.sha
import java.util.*

private const val ID_SEPARATOR = '_'
private val ID_SPLITTERS = charArrayOf(ID_SEPARATOR, ',', '.')

private fun getSubChunkOrder(chunks: List<Set<String>>): Multimap<Int, String> {
    // Lets collect map predicate -> summary number of usages
    val numbersMap = TObjectIntHashMap<String>()
    chunks.flatten().forEach { p -> numbersMap.adjustOrPutValue(p, 1, 1) }
    val result = Multimaps.newSetMultimap(Maps.newHashMap<Int, Collection<String>>()) { newHashSet() }
    numbersMap.forEachEntry { a, b ->
        result.put(b, a)
        true
    }
    return result
}

/**
 * Collect all the sub chunks separated by underscore
 */
private fun collectSubChunks(p: String): Set<String> {
    val cumulativeLengths = p.split(*ID_SPLITTERS).stream().mapToInt { it.length }.toArray()
            .let { Arrays.parallelPrefix(it) { a, b -> a + b }; it }
    for (i in 0 until cumulativeLengths.size) {
        cumulativeLengths[i] = cumulativeLengths[i] + i
    }
    val result = hashSetOf<String>()
    for (i in -1 until cumulativeLengths.size) {
        for (j in i + 1 until cumulativeLengths.size) {
            val start = if (i == -1) 0 else cumulativeLengths[i] + 1
            val end = cumulativeLengths[j]
            result.add(p.substring(start, end))
        }
    }
    return result
}

/**
 * Constructs a short and human-readable identifier for a list of ids.
 *
 * The strategy is to
 *   0. Try as is
 *   1. Try to reduce IDs to GSM IDs
 *   2. Find common subsequence limited by _
 *   3. Select subsequence applicable to max number of ids > 2 and remove it
 *   4. Iterate if necessary
 *   5. Concatenate the resulting strings.
 *   6. Truncate and append hash if still too long
 */
fun reduceIds(files: List<String>, maxLength: Int = 100): String {
    val sha = files.joinToString(ID_SEPARATOR.toString()).sha
    val names = files.map { it.substringBefore(HASH_PREFIX) }.distinct()
    // Try to join as is
    val fullId = names.joinToString(ID_SEPARATOR.toString())
    if (fullId.length < maxLength) {
        return fullId
    }
    // Try to filter our GSM ids
    val gsmRegex = "GSM\\d+".toRegex()
    val gsmIds = names.map {  gsmRegex.find(it)?.value ?: it }
    if (gsmIds != names) {
        return reduceIds(gsmIds, maxLength)
    }

    var ids = names
    // Filter case of inclusion
    val filtered = ids.filter { s -> ids.count { s in it } == 1 }
    ids = filtered
    while (true) {
        // Iterative algorithm, start with greatest number of common sub chunk, remove it etc
        val map = getSubChunkOrder(ids.map(::collectSubChunks))
        val max = map.keySet().max()!!
        val common = map.get(max).sorted().first()
        if (max < 2 || ids.any { it.length == common.length }) {
            break
        }
        ids = listOf(common) + ids.map {
            var result = it.replace(common, "")
            for (splitter in ID_SPLITTERS) {
                result = result.replace("$splitter$splitter", "$splitter").trim(splitter)
            }
            result
        }
    }
    val joinedId = ids.joinToString(ID_SEPARATOR.toString())
    return if (joinedId.length <= maxLength) {
        joinedId
    } else {
        "${joinedId.substring(0, Math.min(joinedId.length, maxLength - sha.length))}$sha"
    }
}