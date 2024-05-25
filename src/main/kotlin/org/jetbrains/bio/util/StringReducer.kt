package org.jetbrains.bio.util

import kotlin.math.min

/**
 * Reduces list of string into truncated id with hash
 */
fun reduceIds(strings: List<String>, maxLength: Int = 100, predefinedSha: String? = null): String {
    val sha = predefinedSha ?: strings.joinToString("_").sha
    require(maxLength > sha.length) { "Too short length required" }
    // Reduce possible sha in names and fix illegal symbols
    val ids = strings.map {
        it
            .substringAfterLast(if (isWindows()) "\\" else "/")
            .substringBefore(HASH_PREFIX)
            .replace(Regex("[^a-zA-Z0-9_]+"), "_")
            .replace(Regex("__+"), "_")
            .replace(Regex("(^_+)|(_+$)"), "")
    }.distinct()
    var joined = ids.joinToString("_")
    if (joined.length <= maxLength) {
        return joined
    }
    joined += sha
    // Try to filter our GSM ids
    val gsmRegex = "GSM\\d+".toRegex()
    val gsmIds = ids.map { gsmRegex.find(it)?.value ?: it }
    if (gsmIds != ids) {
        return reduceIds(gsmIds, maxLength, sha)
    }
    return "${joined.substring(0, min(joined.length, maxLength - sha.length))}$sha"
}