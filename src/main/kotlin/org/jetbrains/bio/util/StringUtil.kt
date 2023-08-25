package org.jetbrains.bio.util

import org.apache.commons.math3.util.Precision

fun Long.formatLongNumber() = String.format("%,d", this).replace(',', ' ')
fun Int.formatLongNumber() = this.toLong().formatLongNumber()

fun Int.asPercentOf(total: Int, digitsAfterDot: Int = 2) = Precision.round(100.0 * this / total, digitsAfterDot)
