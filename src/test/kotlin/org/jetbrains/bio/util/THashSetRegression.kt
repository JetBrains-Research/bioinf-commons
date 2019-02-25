package org.jetbrains.bio.util

import gnu.trove.set.hash.TIntHashSet
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.fail

class THashSetRegression {
    @Test fun retainAllSideEffect() {
        val data = intArrayOf(5, 4, 3)

        val seen = TIntHashSet(intArrayOf(1, 2, 3))

        // Precondition
        assertEquals("[5, 4, 3]", Arrays.toString(data))

        seen.retainAll(data)

        // Test
        //assertEquals("[5, 4, 3]", Arrays.toString(data));
        if (Arrays.toString(data) == "[5, 4, 3]") {
            // https://bitbucket.org/robeden/trove/issue/53/_e_hashsettemplate-retainall-hidden-side
            fail("Issue #53 already fixed, time to remove explicit array.clone in DataFrame Column.intersect()")
        } else {
            // current side effect demo:
            assertEquals("[3, 4, 5]", Arrays.toString(data))
        }
    }
}
