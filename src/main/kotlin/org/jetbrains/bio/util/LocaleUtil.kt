package org.jetbrains.bio.util

import java.util.*

/**
 * Is used for writing tests when numeric values are converted to strings with respect to user locale
 */
fun <T> withLocale(locale: Locale, action: () -> T): T {
    val originLocale = Locale.getDefault()

    return try {
        Locale.setDefault(locale)
        action()
    } finally {
        Locale.setDefault(originLocale)
    }
}

fun withLocaleEU(action: () -> Unit) {
    // Locale where decimal separator is '.' and a thousand separator is ','
    // Here we expect `1234.56789` may to look like `1.234,56789` (DE locale)
    // But also could be `1 234,56789` in `fr_FR` locale
    withLocale(locale = Locale.GERMANY, action)
}

fun withLocaleUS(action: () -> Unit) {
    // Locale where a decimal separator is ',' and a thousand separator is '.'
    // Here we expect `1234.56789` may to look like `1,234.56789`
    withLocale(locale = Locale.US, action)
}
