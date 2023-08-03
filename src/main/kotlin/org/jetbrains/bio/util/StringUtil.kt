package org.jetbrains.bio.util

fun Long.formatLongNumber() = String.format("%,d", this).replace(',', ' ')
fun Int.formatLongNumber() = this.toLong().formatLongNumber()
