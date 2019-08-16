package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.exception.DimensionMismatchException
import org.apache.commons.math3.exception.MathIllegalArgumentException
import org.apache.commons.math3.exception.NoDataException
import org.apache.commons.math3.exception.util.LocalizedFormats
import org.apache.commons.math3.linear.*
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor._I
import org.jetbrains.bio.viktor.asF64Array
import org.jetbrains.bio.viktor.toF64Array

class WLSMultipleLinearRegression {
    /** Entries of the matrix.  */
    private var xMatrix: F64Array = DoubleArray(0).asF64Array()
    private var yVector: F64Array = DoubleArray(0).asF64Array()
    //use DiagonalMatrix instead of RealMatrix
    private var Omega: F64Array = DoubleArray(0).asF64Array()

    fun getX(): RealMatrix {
       return Array2DRowRealMatrix(Array (xMatrix.shape[0]) { xMatrix.V[it].toDoubleArray()})
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
    fun newXSampleData(x: Array<DoubleArray>) {
        if (x.size == 0) {
            throw NoDataException()
        }

        xMatrix = DoubleArray(x[0].size){1.0}.asF64Array().T.append(x.toF64Array().T, 1)
        val l = 0
    }

    //override
    fun newYSampleData(y: DoubleArray) {
        if (y.size == 0) {
            throw NoDataException()
        }
        this.yVector = y.asF64Array()
    }

    protected fun validateCovarianceData(weights: DoubleArray) {
        if (this.xMatrix.shape[0] != weights.size) {
            throw DimensionMismatchException(this.xMatrix.shape[0], weights.size)
        }
    }

    protected fun validateSampleData(y: DoubleArray) {
        if (this.xMatrix.shape[0] != y.size) {
            throw DimensionMismatchException(y.size, this.xMatrix.shape[0])
        }
        if (this.xMatrix.shape[0] == 0) {  // Must be no y data either
            throw NoDataException()
        }
        if (this.xMatrix.shape[1] > this.xMatrix.shape[0]) {
            throw MathIllegalArgumentException(
                    LocalizedFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS,
                    this.xMatrix.shape[0], this.xMatrix.shape[1])
        }
    }

    //use DiagonalMatrix
    protected fun newCovarianceData(omega: DoubleArray) {
        this.Omega = omega.asF64Array()
    }

    fun calculateBeta(): RealVector {
        var XTOTX = Array (xMatrix.shape[1]) { i ->
            DoubleArray(xMatrix.shape[1]) {j ->
                xMatrix.V[_I, i].times(Omega).dot(xMatrix.V[_I, j])
            }
        }
        val inverse = LUDecomposition(Array2DRowRealMatrix(XTOTX)).solver.inverse
        val Oy = Omega.times(yVector)
        val XTOy = F64Array (xMatrix.shape[1]) {xMatrix.V[_I, it].dot(Oy)}
        return ArrayRealVector(DoubleArray (xMatrix.shape[1]) {inverse.getRow(it).asF64Array().dot(XTOy)}, false)
    }
}