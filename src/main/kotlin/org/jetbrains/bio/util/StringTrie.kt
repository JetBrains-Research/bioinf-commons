package org.jetbrains.bio.util

class StringTrieNode(
        val value: Char = 0.toChar(),
        val parent: StringTrieNode? = null,
        var child: StringTrieNode? = null,
        var sibling: StringTrieNode? = null
) {

    operator fun get(c: Char): StringTrieNode {
        if (child == null) {
            val newNode = StringTrieNode(c, this)
            child = newNode
            return newNode
        }
        var node = child!!
        while (true) {
            if (node.value == c) {
                return node
            }
            if (node.sibling == null) {
                val newNode = StringTrieNode(c, this)
                node.sibling = newNode
                return newNode
            }
            node = node.sibling!!
        }
    }

    fun value(): String = buildString {
        var node = this@StringTrieNode
        while (node.parent != null) {
            insert(0, node.value)
            node = node.parent!!
        }
    }

    operator fun get(s: String): StringTrieNode = getRecursive(s, 0, this)

    companion object {
        tailrec fun getRecursive(s: String, offset: Int, node: StringTrieNode): StringTrieNode {
            if (offset == s.length) {
                return node
            }
            val nextNode = node[s[offset]]
            return getRecursive(s, offset + 1, nextNode)
        }
    }

}
