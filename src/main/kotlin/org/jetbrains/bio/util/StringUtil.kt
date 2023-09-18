package org.jetbrains.bio.util

import java.util.*

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
