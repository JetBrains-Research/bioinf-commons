package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.exception.NotStrictlyPositiveException
import org.apache.commons.math3.exception.OutOfRangeException
import org.apache.commons.math3.linear.AbstractRealMatrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealMatrix

/**
 * Matrix class for WLSMultipleLinearRegession that allows to use non transposed data matrix without intercept.
 * Save us from unnecessary copying.
 *
 * @author Elena Kartysheva
 * @date 7/8/2019
 */
class DesignMatrix(
        /** non-transposed data matrix without intercept */
        private val data: Array<DoubleArray>) : AbstractRealMatrix() {

    /** @inheritDoc */
    override fun getData(): Array<DoubleArray> {
        return this.data
    }

    /**
     * Get number of columns of transposed data with intercept
     */
    override fun getColumnDimension(): Int = data.size + 1

    /**
     * Get number of rows of transposed data
     */
    override fun getRowDimension(): Int = data[0].size

    /** @inheritDoc */
    override fun copy(): RealMatrix = DesignMatrix(copyOut())

    /**
     * Get a fresh copy of the underlying data array.
     *
     * @return a copy of the underlying data array.
     */
    private fun copyOut(): Array<DoubleArray> {
        val nRows = this.rowDimension
        val out = Array(nRows) { DoubleArray(this.columnDimension) }
        // can't copy 2-d array in one shot, otherwise get row references
        for (i in 0 until nRows) {
            System.arraycopy(data[i], 0, out[i], 0, data[i].size)
        }
        return out
    }

    /** @inheritDoc */
    @Throws(NotStrictlyPositiveException::class)
    override fun createMatrix(rowDimension: Int,
                              columnDimension: Int): RealMatrix {
        return Array2DRowRealMatrix(rowDimension, columnDimension)
    }

    /**
     * Get element of transposed data matrix with intercept.
     */
    @Throws(OutOfRangeException::class)
    override fun getEntry(row: Int, column: Int): Double {
        if (column == 0) return 1.0
        MatrixUtils.checkMatrixIndex(this, row, column - 1)
        return data[column - 1][row]
    }

    /**
     * Set element of transposed data matrix with intercept.
     */
    @Throws(OutOfRangeException::class)
    override fun setEntry(row: Int, column: Int, value: Double) {
        if (column > 0) {
            MatrixUtils.checkMatrixIndex(this, column - 1, row)
            data[column - 1][row] = value
        }
    }


}