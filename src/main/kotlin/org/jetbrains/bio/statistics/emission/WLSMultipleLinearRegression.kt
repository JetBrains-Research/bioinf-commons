package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.exception.DimensionMismatchException
import org.apache.commons.math3.exception.MathIllegalArgumentException
import org.apache.commons.math3.exception.NoDataException
import org.apache.commons.math3.exception.util.LocalizedFormats
import org.apache.commons.math3.linear.*
import org.apache.commons.math3.stat.regression.AbstractMultipleLinearRegression

class WLSMultipleLinearRegression : AbstractMultipleLinearRegression() {
    /** Entries of the matrix.  */
    private var xMatrix: DesignMatrix = DesignMatrix(Array(0){ DoubleArray(0) })
    private var yVector: RealVector = ArrayRealVector()
    //use DiagonalMatrix instead of RealMatrix
    private var Omega: DiagonalMatrix = DiagonalMatrix(DoubleArray(0))

    public override fun getX(): RealMatrix {
        return xMatrix
    }

    fun newSampleData(y: DoubleArray, x: Array<DoubleArray>, weights: DoubleArray) { //weight now is double[] instead of double[][]
        this.newXSampleData(x)
        this.validateSampleData(y)
        this.newYSampleData(y)
        this.validateCovarianceData(weights)
        this.newCovarianceData(weights)
    }

    fun newSampleData(y: DoubleArray, weights: DoubleArray) {
        this.validateSampleData(y)
        this.newYSampleData(y)
        this.validateCovarianceData(weights)
        this.newCovarianceData(weights)
    }

    //override
    override fun newXSampleData(x: Array<DoubleArray>) {
        if (x.size == 0) {
            throw NoDataException()
        }

        xMatrix = DesignMatrix(x)
    }

    //override
    override fun newYSampleData(y: DoubleArray) {
        if (y.size == 0) {
            throw NoDataException()
        }
        this.yVector = ArrayRealVector(y)
    }

    protected fun validateCovarianceData(weights: DoubleArray) {
        if (this.xMatrix.rowDimension != weights.size) {
            throw DimensionMismatchException(this.xMatrix.rowDimension, weights.size)
        }
    }

    protected fun validateSampleData(y: DoubleArray) {
        if (this.xMatrix.rowDimension != y.size) {
            throw DimensionMismatchException(y.size, this.xMatrix.rowDimension)
        }
        if (this.xMatrix.rowDimension == 0) {  // Must be no y data either
            throw NoDataException()
        }
        if (this.xMatrix.columnDimension > this.xMatrix.rowDimension) {
            throw MathIllegalArgumentException(
                    LocalizedFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS,
                    this.xMatrix.rowDimension, this.xMatrix.columnDimension)
        }
    }

    //use DiagonalMatrix
    protected fun newCovarianceData(omega: DoubleArray) {
        this.Omega = DiagonalMatrix(omega)
    }

    public override fun calculateBeta(): RealVector {
        val XT = xMatrix.transpose()
        val XTOTX = XT.multiply(Omega.multiply(xMatrix))
        val inverse = LUDecomposition(XTOTX).solver.inverse
        return inverse.multiply(XT).operate(Omega.operate(yVector))
    }

    public override fun calculateBetaVariance(): RealMatrix {
        val XTOX = xMatrix.transpose().multiply(Omega.multiply(xMatrix))
        return LUDecomposition(XTOX).solver.inverse
    }
}