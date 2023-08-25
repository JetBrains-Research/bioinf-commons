package org.jetbrains.bio.gsea

import joptsimple.ValueConversionException
import joptsimple.ValueConverter

enum class PermutationAltHypothesis(private val presentableString: String) {
    TWO_SIDED("two-sided"),
    GREATER("greater"),
    LESS("less");

    override fun toString() = presentableString

    companion object {
        fun converter() = object : ValueConverter<PermutationAltHypothesis> {
            override fun convert(value: String): PermutationAltHypothesis {
                for (h in values()) {
                    if (h.presentableString == value) {
                        return h
                    }
                }
                throw ValueConversionException("Unsupported hypothesis: $value")
            }

            override fun valueType() = PermutationAltHypothesis::class.java

            override fun valuePattern() = null
        }
    }
}