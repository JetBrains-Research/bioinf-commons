package org.jetbrains.bio.util

import gnu.trove.set.hash.TIntHashSet
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.fail

class THashSetRegression {
    @Test
    fun retainAllSideEffect() {
        val data = intArrayOf(5, 4, 3)

        val seen = TIntHashSet(intArrayOf(1, 2, 3))

        // Precondition
        assertEquals("[5, 4, 3]", data.contentToString())

        seen.retainAll(data)

        // Test
        //assertEquals("[5, 4, 3]", Arrays.toString(data));
        if (data.contentToString() == "[5, 4, 3]") {
            // NOTE[shpynov]
            // This issue is marked as closed in trove4j 3.2, which is not available as of Feb 28th 2019
            // Bump trove4j version once it is available
            // https://bitbucket.org/trove4j/trove/issues/53/_e_hashsettemplate-retainall-hidden-side
            fail("Issue #53 already fixed, time to remove explicit array.clone in DataFrame Column.intersect()")
        } else {
            // current side effect demo:
            assertEquals("[3, 4, 5]", data.contentToString())
        }
    }
}
