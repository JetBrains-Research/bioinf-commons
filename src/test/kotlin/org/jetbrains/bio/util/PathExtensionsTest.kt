@file:Suppress("UsePropertyAccessSyntax")

package org.jetbrains.bio.util

import org.jetbrains.bio.genome.format.BedParserTest
import org.junit.Test
import java.io.File
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.CountDownLatch
import java.util.concurrent.Executors
import java.util.concurrent.atomic.AtomicInteger
import java.util.zip.GZIPInputStream
import java.util.zip.ZipInputStream
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue
import kotlin.test.fail

class PathExtensionsTest {
    @Test
    fun name() {
        assertEquals("to", "/path/to/".toPath().name)
        assertEquals("foo", "/path/to/foo".toPath().name)
        assertEquals("foo.bar", "/path/to/foo.bar".toPath().name)
    }

    @Test
    fun escape() {
        assertEquals("h3k27me3_log2_od_yd", "[H3K27me3 log2(OD/YD)]".escape())
        assertEquals("h3k27me3_log2_od_yd.bw", "[H3K27me3 log2(OD/YD)]().bw".escape())
        assertEquals("h3k27me3_od_mean.foo.bw", "H3K27me3 OD mean.foo.bw".escape())
    }


    @Test
    fun stem() {
        assertEquals("to", "/path/to/".toPath().stem)
        assertEquals("foo", "/path/to/foo".toPath().stem)
        assertEquals("foo", "/path/to/foo.bar".toPath().stem)
    }

    @Test
    fun stemMultipleExtensions() {
        assertEquals("foo.bar", "/path/to/foo.bar.baz".toPath().stem)

        // GZ and ZIP extensions are treated as all other extensions:
        assertEquals("foo.OD7", "/path/to/foo.OD7.zip".toPath().stem)
        assertEquals("foo.OD7", "/path/to/foo.OD7.gz".toPath().stem)
        assertEquals("foo.input.tar", "/path/to/foo.input.tar.gz".toPath().stem)
        assertEquals("foo.input.tar", "/path/to/foo.input.tar.zip".toPath().stem)
    }

    @Test
    fun stemAndExtConsistency() {
        val assertConsistent: (String) -> Unit = { pathStr ->
            val path = pathStr.toPath()
            assertEquals(path.name, "${path.stem}.${path.extension}")
        }

        assertConsistent("/path/to/foo.bar")
        assertConsistent("/path/to/foo.bar.baz")

        // GZ & Zip case
        assertConsistent("/path/to/foo.OD7.zip")
        assertConsistent("/path/to/foo.OD7.gz")
        assertConsistent("/path/to/foo.bar.tar.gz")
        assertConsistent("/path/to/foo.bar.tar.zip")
        assertConsistent("/path/to/input.gz.OD17.model.tar.zip")
        assertConsistent("/path/to/input.gz.OD17.model.tar.gz")
    }

    @Test
    fun hasChunk() {
        assertTrue("/path/to/foo".toPath().hasChunk("foo"))
        assertTrue("/path/to/bar_foo_baz".toPath().hasChunk("foo"))
        assertTrue("/path/to/bar.foo.baz".toPath().hasChunk("foo"))
        assertFalse("/path/to/foo12".toPath().hasChunk("foo"))
    }


    @Test
    fun withName() {
        assertEquals(
            "/path/to/boo".toPath(),
            "/path/to/foo.bar".toPath().withName("boo")
        )
    }

    @Test
    fun withExtension() {
        assertEquals(
            "/path/to/foo.baz".toPath(),
            "/path/to/foo.bar".toPath().withExtension("baz")
        )
        assertEquals(
            "/path/to/foo.e.baz".toPath(),
            "/path/to/foo.e.bar".toPath().withExtension("baz")
        )
        assertEquals(
            "/path/to/foo.baz".toPath(),
            "/path/to/foo".toPath().withExtension("baz")
        )
    }

    @Test
    fun withStem() {
        assertEquals(
            "/path/to/bar.baz".toPath(),
            "/path/to/foo.baz".toPath().withStem("bar")
        )
        assertEquals(
            "/path/to/bar.baz".toPath(),
            "/path/to/foo.baz".toPath().withStem("bar")
        )
        assertEquals(
            "/path/to/bar".toPath(),
            "/path/to/foo".toPath().withStem("bar")
        )
    }

    @Test
    fun glob() {
        withTempDirectory("glob") { path ->
            assertEquals(emptyList(), path.glob("*"))

            val dummy = byteArrayOf()
            (path / "foo.txt").write(dummy)
            (path / "bar.zip").write(dummy)

            assertEquals(path.list().toSet(), path.glob("*").toSet())
            assertEquals(listOf(path / "foo.txt"), path.glob("*.txt"))
        }
    }

    @Test
    fun globRecursive() {
        withTempDirectory("glob") { path ->
            assertEquals(emptyList(), path.glob("*"))
            val subdir = path / "foo" / "bar"
            subdir.createDirectories()

            val dummy = byteArrayOf()
            (subdir / "boo.txt").write(dummy)
            (path / "baz.zip").write(dummy)

            assertEquals(setOf(subdir / "boo.txt", path / "baz.zip"), path.glob("**").toSet())
            assertEquals(setOf(subdir / "boo.txt"), path.glob("**/*.txt").toSet())
        }
    }

    @Test
    fun divStringString() {
        assertEquals("foo${File.separatorChar}boo", ("foo" / "boo").toString())
    }

    @Test
    fun divStringPath() {
        assertEquals("foo${File.separatorChar}boo", ("foo" / "boo".toPath()).toString())
    }

    @Test
    fun divPathString() {
        assertEquals("foo${File.separatorChar}boo", ("foo".toPath() / "boo").toString())
    }

    @Test
    fun divPathPath() {
        assertEquals("foo${File.separatorChar}boo", ("foo".toPath() / "boo".toPath()).toString())
    }

    @Test
    fun listFiles_Dir() {
        val actual = withTempDirectory("foo") { dir ->
            (dir / "file1.txt").touch()
            (dir / "file2.txt").touch()
            (dir / "boo" / "boo1.txt")
            (dir / "doo")

            dir.list().map { it.name }
        }

        assertEquals(listOf("file1.txt", "file2.txt"), actual.sorted())
    }

    @Test
    fun gzipInputOutput() {
        withTempFile("test", ".gz") { path ->
            path.bufferedWriter().use { it.write("Hello, world!") }
            path.bufferedReader().use {
                assertEquals("Hello, world!", it.readText())
            }
        }
    }

    @Test
    fun zipInputOutput() {
        withTempFile("test", ".zip") { path ->
            path.bufferedWriter().use { it.write("Hello, world!") }
            path.bufferedReader().use {
                assertEquals("Hello, world!", it.readText())
            }
        }
    }
}

class FileSizeTest {
    @Test
    fun formatting() {
        assertEquals("<not accessible>", FileSize(-1).toString())
        assertEquals("0 b", FileSize(0L).toString())
        assertEquals("10 b", FileSize(10L).toString())
        assertEquals("100 b", FileSize(100L).toString())
        assertEquals("1020 b", FileSize(1020L).toString())
        assertEquals("1 kb", FileSize(1050L).toString())
        assertEquals("10,3 kb", FileSize(10500L).toString())
        assertEquals("102,5 kb", FileSize(105000L).toString())
        assertEquals("1 mb", FileSize(1050000L).toString())
        assertEquals("10 mb", FileSize(10500000L).toString())
        assertEquals("100,1 mb", FileSize(105000000L).toString())
        assertEquals("1001,4 mb", FileSize(1050000000L).toString())
        assertEquals("1,9 gb", FileSize(2050000000L).toString())
    }
}

class CheckOrRecalculateTest {
    @Test
    fun checkOrRecalculateMissing() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked.txt"
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id") { output ->
                    Thread.sleep(10)
                    output.path.write(byteArrayOf(0, 1, 2))
                    recalcFlag.andIncrement
                }
            }
            assertEquals(1, recalcFlag.get(), "A single block should've been executed!")
        }
    }

    @Test
    fun checkOrRecalculateMissingDir() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked_dir"
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id", isDirectory = true) { output ->
                    Thread.sleep(10)
                    (output.path / "result.txt").write(byteArrayOf(0, 1, 2))
                    recalcFlag.andIncrement
                }
            }
            assertEquals(1, recalcFlag.get(), "A single block should've been executed!")
        }
    }

    @Test
    fun checkOrRecalculateExisting() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked.txt"
            lockedPath.write(byteArrayOf(0, 1, 2))
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id") { _ ->
                    Thread.sleep(10)

                    // No write leave file empty
                    recalcFlag.andIncrement
                }
            }
            assertEquals(0, recalcFlag.get())
        }
    }

    @Test
    fun checkOrRecalculateExistingIgnoreEmptyFile() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val path = it / "path.txt"
            path.touch()
            assertTrue(path.exists)
            val recalculate: (PathWrapper) -> Unit = { (p) ->
                Thread.sleep(10)
                // No write leave file empty
                recalcFlag.andIncrement
                p.touch()
            }
            path.checkOrRecalculate(ignoreEmptyFile = true, recalculate = recalculate)
            assertEquals(0, recalcFlag.get())
            path.delete()
            assertFalse(path.exists)
            path.checkOrRecalculate(ignoreEmptyFile = true, recalculate = recalculate)
            assertEquals(1, recalcFlag.get())
            assertTrue(path.exists)
        }
    }


    @Test
    fun checkOrRecalculateExistingDir() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked_dir"
            lockedPath.createDirectories()
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id", isDirectory = true) { _ ->
                    Thread.sleep(10)

                    // No write leave file empty
                    recalcFlag.andIncrement
                }
            }
            assertEquals(0, recalcFlag.get())
        }
    }

    @Test
    fun checkOrRecalculateEmptyDir() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked_dir"
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id", isDirectory = true) { _ ->
                    Thread.sleep(10)

                    // No write leave file empty
                    recalcFlag.andIncrement
                }
            }
            assertEquals(1, recalcFlag.get())
        }
    }

    @Test
    fun stemExtension() {
        withTempDirectory("annotations") {
            val lockedPath = it / "locked.txt"
            lockedPath.checkOrRecalculate("test") { output ->
                assertTrue(output.path.stem.startsWith(lockedPath.stem))
                assertEquals("txt", output.path.extension)

                output.path.write(byteArrayOf(1, 2, 3))
            }
        }
    }

    @Test(expected = IllegalStateException::class)
    fun ignoredOutputPath() {
        withTempDirectory("annotations") {
            (it / "foo.txt").checkOrRecalculate("") {
                // ignore path wrapper!
            }
        }
    }

    @Test
    fun accessedOutputPath() {
        withTempDirectory("annotations") {
            (it / "foo.txt").checkOrRecalculate("") {
                it.path.write(byteArrayOf(1, 2, 3))
            }
        }
    }

    @Test(expected = IllegalStateException::class)
    fun zeroSizedOutputPath() {
        val recalcFlag = AtomicInteger(0)

        withTempDirectory("annotations") {
            val lockedPath = it / "locked.txt"
            execConcurrently { id ->
                lockedPath.checkOrRecalculate("thread #$id") { output ->
                    Thread.sleep(10)

                    output.path  // access, but don't write anything.
                    recalcFlag.getAndIncrement()
                }

            }
            assertEquals(2, recalcFlag.get())
        }
    }

    @Test
    fun missingParent() {
        withTempDirectory("annotations") {
            val lockedPath = it / "foo" / "bar" / "locked.txt"
            lockedPath.checkOrRecalculate("test") { output ->
                output.path.write(byteArrayOf(1, 2, 3))
            }

            assertEquals(3, lockedPath.size.bytes)
        }
    }


    @Test
    fun checkOrRecalculateName() {
        withTempDirectory("foo") { dir ->
            (dir / "bar.txt").checkOrRecalculate(ignoreEmptyFile = true) { (p) ->
                val fileName = p.fileName.toString()
                assertTrue(fileName.startsWith("bar_tmp"))
                assertTrue(fileName.endsWith(".txt"))
            }
        }
    }

    @Test
    fun checkOrRecalculateNameDir() {
        withTempDirectory("foo") { dir ->
            (dir / "bar").checkOrRecalculate(isDirectory = true) { (p) ->
                assertTrue(p.fileName.toString().startsWith("bar_tmp"))
            }
        }
    }

    @Test
    fun checkFailedDirDeleted() {
        var path: Path? = null
        try {
            withTempDirectory("fail") { dir ->
                path = dir
                fail("FAILED")
            }
        } catch (t: Throwable) {
            assertFalse(path!!.exists)
            return
        }
    }

    @Test
    fun streamFor() {
        withResource(BedParserTest::class.java, "bed12.bed.gz") { path ->
            assertTrue(Files.newInputStream(path).streamFor(path.name) is GZIPInputStream)
            assertTrue(Files.newInputStream(path).streamFor(path.toString()) is GZIPInputStream)
        }

        withResource(BedParserTest::class.java, "bed12.bed.zip") { path ->
            assertTrue(Files.newInputStream(path).streamFor(path.name) is ZipInputStream)
            assertTrue(Files.newInputStream(path).streamFor(path.toString()) is ZipInputStream)
        }

        withResource(BedParserTest::class.java, "bed12.bed") { path ->
            assertFalse(Files.newInputStream(path).streamFor(path.name) is GZIPInputStream)
            assertFalse(Files.newInputStream(path).streamFor(path.name) is ZipInputStream)
            assertFalse(Files.newInputStream(path).streamFor(path.toString()) is GZIPInputStream)
            assertFalse(Files.newInputStream(path).streamFor(path.toString()) is ZipInputStream)
        }
    }

    private inline fun execConcurrently(crossinline task: (Int) -> Unit) {
        val availableProcessors = Runtime.getRuntime().availableProcessors()
        if (availableProcessors == 1) {
            Logs.getRootLogger().warn("Cannot do honest test if no parallelism. Only on processor is available")
        }

        val executor = Executors.newFixedThreadPool(2)
        val latch = CountDownLatch(2)  // ensures concurrency.
        val tasks = (1..2).map {
            Callable {
                latch.countDown()
                task(it)
            }
        }

        executor.awaitAll(tasks)
        check(executor.shutdownNow().isEmpty())
    }

    @Test
    fun checkWithTempDirectoryMissingDir() {
        withTempDirectory("foo") { dir ->
            withTempDirectory("bar", dir / "missing") {
                assertTrue(it.isDirectory)
            }
        }
    }

    @Test
    fun checkWithTempFileMissingDir() {
        withTempDirectory("foo") { dir ->
            withTempFile("bar", ".txt", dir / "missing") {
                assertTrue(it.isRegularFile)
            }
        }
    }

    @Test
    fun checkPathToURIAndBack() {
        val homeDir = System.getProperty("user.home")
        assertEquals(homeDir, homeDir.toUri().path.toPath().toString())
    }

}