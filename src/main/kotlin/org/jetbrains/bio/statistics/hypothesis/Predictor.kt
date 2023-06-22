package org.jetbrains.bio.statistics.hypothesis

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.viktor.F64Array

interface Predictor {
    fun predict(logNullMemberships: F64Array): BitList
}