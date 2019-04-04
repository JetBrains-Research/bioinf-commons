package org.jetbrains.bio.methylome

import com.google.common.base.MoreObjects
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame

/**
 * A view of a single chromosome/strand data.
 *
 * @author Sergei Lebedev
 */
data class StrandMethylomeView internal constructor(
        private val frame: MethylomeFrame) : MethylomeView {

    override fun peel() = frame.peel()

    override val size = frame.size
}

/**
 * A view of strand-combined data for a single chromosome.
 *
 * @author Sergei Lebedev
 */
data class ChromosomeMethylomeView internal constructor(
        private val framePlus: MethylomeFrame,
        private val frameMinus: MethylomeFrame) : MethylomeView {

    override fun peel() = DataFrame.rowBind(
            framePlus.peel().with("strand", BitterSet(framePlus.size) { true }),
            frameMinus.peel().with("strand", BitterSet(frameMinus.size) { false })
    ).reorder("offset")

    override val size = framePlus.size + frameMinus.size

    override fun toString() = MoreObjects.toStringHelper(this)
            .add("+", framePlus.size)
            .add("-", frameMinus.size)
            .toString()
}

/**
 * An immutable view of [Methylome] data.
 *
 * @author Sergei Lebedev
 */
interface MethylomeView {
    /**
     * Extracts a data frame with view contents.
     *
     * @see Methylome.MethylomeFrame for concrete implementations.
     * @see MethylomeToDataFrame for a cytosine context-specific way
     *                           of extracting a data frame.
     */
    fun peel(): DataFrame

    val size: Int
}
