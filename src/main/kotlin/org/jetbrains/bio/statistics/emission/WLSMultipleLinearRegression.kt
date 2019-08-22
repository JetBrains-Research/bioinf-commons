package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.exception.DimensionMismatchException
import org.apache.commons.math3.exception.MathIllegalArgumentException
import org.apache.commons.math3.exception.NoDataException
import org.apache.commons.math3.exception.util.LocalizedFormats
import org.apache.commons.math3.linear.*
import org.apache.commons.math3.stat.regression.AbstractMultipleLinearRegression
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array

class WLSMultipleLinearRegression : AbstractMultipleLinearRegression() {
    /** Entries of the matrix.  */
    private var xMatrix: DesignMatrix = DesignMatrix(Array(0){ DoubleArray(0) })
    private var yVector: F64Array = DoubleArray(0).asF64Array()
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
        this.yVector = y.asF64Array()
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
        var XTOTX = Array (xMatrix.columnDimension) { i ->
            DoubleArray(xMatrix.columnDimension) {j ->
                xMatrix.getColumn(i).asF64Array().apply { timesAssign(Omega.dataRef.asF64Array()) }.dot(xMatrix.getColumn(j))
            }}
        val inverse = LUDecomposition(Array2DRowRealMatrix(XTOTX)).solver.inverse
        val Oy = Omega.dataRef.asF64Array().apply { timesAssign(yVector) }
        val XTOy = F64Array (xMatrix.columnDimension) {xMatrix.getColumn(it).asF64Array().dot(Oy)}
        return ArrayRealVector(DoubleArray (xMatrix.columnDimension) {inverse.getRow(it).asF64Array().dot(XTOy)}, false)
    }

    public override fun calculateBetaVariance(): RealMatrix {
        val XTOX = xMatrix.transpose().multiply(Omega.multiply(xMatrix))
        return LUDecomposition(XTOX).solver.inverse
    }
}
