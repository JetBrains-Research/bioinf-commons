package org.jetbrains.bio.statistics

import org.jetbrains.bio.viktor.F64Array

fun F64Array.Companion.ones(vararg shape: Int): F64Array {
    return full(*shape, init = 1.0)
}

fun F64Array.Companion.stochastic(vararg shape: Int): F64Array {
    return full(*shape, init = 1.0 / shape.last())
}