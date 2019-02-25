package org.jetbrains.bio.util

import org.junit.Test
import kotlin.test.assertEquals

class TarTest {
    @Test
    fun checkPackUnpack() {
        withTempDirectory("foo") { dir ->
            val a = dir / "a.txt"
            a.bufferedWriter().use {
                it.write("Technology is not science!")
            }
            val b = dir / "b.txt"
            b.bufferedWriter().use {
                it.write("Goodbye cruel world!")
            }
            val p = dir / "package.tar"
            Tar.compress(p, a.toFile(), b.toFile())
            val uncompressed = dir / "uncompressed"
            Tar.decompress(p, uncompressed.toFile())
            assertEquals("[a.txt, b.txt]", uncompressed.list().map { it.fileName }.sortedBy { it.name }.toString())
            assertEquals("Technology is not science!", (uncompressed / "a.txt").useLines { it.firstOrNull()!! })
        }
    }

    @Test
    fun checkSkipExistingTrue() {
        withTempDirectory("foo") { dir ->
            val a = dir / "a.txt"
            a.bufferedWriter().use {
                it.write("Technology is not science!")
            }
            val b = dir / "b.txt"
            b.bufferedWriter().use {
                it.write("Goodbye cruel world!")
            }
            val p = dir / "package.tar"
            Tar.compress(p, a.toFile(), b.toFile())
            val uncompressed = (dir / "uncompressed").apply { createDirectories() }

            (uncompressed / "b.txt").bufferedWriter().use {
                it.write("Hello world!")
            }
            Tar.decompress(p.inputStream(), uncompressed.toFile(), skipExisting = true)
            assertEquals("[a.txt, b.txt]", uncompressed.list().map { it.fileName }.sortedBy { it.name }.toString())
            assertEquals("Hello world!", (uncompressed / "b.txt").useLines { it.firstOrNull()!! })
        }
    }

    @Test
    fun checkSkipExistingFalse() {
        withTempDirectory("foo") { dir ->
            val a = dir / "a.txt"
            a.bufferedWriter().use {
                it.write("Technology is not science!")
            }
            val b = dir / "b.txt"
            b.bufferedWriter().use {
                it.write("Goodbye cruel world!")
            }
            val p = dir / "package.tar"
            Tar.compress(p, a.toFile(), b.toFile())
            val uncompressed = (dir / "uncompressed").apply { createDirectories() }

            (uncompressed / "b.txt").bufferedWriter().use {
                it.write("Hello world!")
            }
            check("Hello world!" == (uncompressed / "b.txt").useLines { it.firstOrNull()!! })
            Tar.decompress(p.inputStream(), uncompressed.toFile(), skipExisting = false)
            assertEquals("Goodbye cruel world!", (uncompressed / "b.txt").useLines { it.firstOrNull()!! })
        }
    }

}
