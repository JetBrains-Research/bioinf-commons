package org.jetbrains.bio.genome.coverage

sealed class Fragment {

    /**
     * At some point, fragment was stored as Int?, with null corresponding to what now is [AutoFragment].
     * Since some parts of code (like ID generation) depend on this convention, this is a useful conversion property.
     */
    abstract val nullableInt: Int?

    companion object {
        @Throws(NumberFormatException::class)
        fun fromString(value: String): Fragment {
            if (value == "auto") {
                return AutoFragment
            }
            return FixedFragment(Integer.parseInt(value))
        }
    }
}

data class FixedFragment(val size: Int) : Fragment() {
    override fun toString(): String = size.toString()
    override val nullableInt: Int = size
}

object AutoFragment : Fragment() {
    override fun toString(): String = "auto"
    override val nullableInt: Int? = null
}

