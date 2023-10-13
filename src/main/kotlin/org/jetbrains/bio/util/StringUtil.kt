package org.jetbrains.bio.util

import org.apache.commons.math3.util.Precision
import java.io.File
import java.util.*

fun String.separatorsToSystem(): String {
    return if (File.separatorChar == '\\') {
        // From Windows to Linux/Mac
        this.replace('/', File.separatorChar)
    } else {
        // From Linux/Mac to Windows
        this.replace('\\', File.separatorChar)
    }
}
fun Long.formatLongNumber(universalThousandsSeparator:Boolean=false): String {
    // XXX: Depending on country the thousands separator could be '.'
    // (e.g. Germany) or ',' (e.g. US).
    //
    // Let's make a solution that will not confuse people how get used to ','
    // as a decimal separator.
    // E.g., when making a screenshot of overlap than `1,034` will look like
    // a floating point number and not like an integer
    if (universalThousandsSeparator) {
        return String
            .format(locale = Locale.US, "%,d", this)
            .replace(',', ' ')
    } else {
        return String.format("%,d", this)
    }
}
fun Int.formatLongNumber(universalThousandsSeparator:Boolean=false) =
    this.toLong().formatLongNumber(universalThousandsSeparator=universalThousandsSeparator)

fun Int.asPercentOf(total: Int, digitsAfterDot: Int = 2) = Precision.round(100.0 * this / total, digitsAfterDot)
fun Int.asPercentStrOf(total: Int, digitsAfterDot: Int = 2) = when (digitsAfterDot) {
    0 -> "${this.asPercentOf(total, digitsAfterDot).toInt()}%"
    else -> String.format("%.${digitsAfterDot}f%%", this.asPercentOf(total, digitsAfterDot))
}
