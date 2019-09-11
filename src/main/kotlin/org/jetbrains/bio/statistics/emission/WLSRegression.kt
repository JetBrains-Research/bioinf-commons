package org.jetbrains.bio.statistics.emission

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.RealMatrix
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array

/**
 * A utility object with various methods useful for WLS and IRLS.
 */
object WLSRegression {

    /**
     * Prepend an appropriate-sized array filled with 1.0.
     */
    fun designMatrix(x: Array<DoubleArray>): Array<DoubleArray> {
        check(x.drop(1).all { it.size == x[0].size }) { "different column sizes" }
        return Array(x.size + 1) { if (it == 0) DoubleArray(x[0].size) { 1.0 } else x[it - 1] }
    }

    /**
     * eta = X beta
     */
    fun calculateEta(x: Array<DoubleArray>, beta: DoubleArray): F64Array {
        check(x.all { it.size == beta.size }) { "X columns and beta have different size" }
        val xColumnsTimesBetaElements = Array(x.size) { i -> x[i].asF64Array() * beta[i] }
        xColumnsTimesBetaElements.drop(1).forEach { xColumnsTimesBetaElements[0].plusAssign(it) }
        return xColumnsTimesBetaElements[0]
    }

    /**
     * beta = (X^T W X)^{-1} X^T W y
     * where X is the design matrix, y is the observation vector, and W = diag(w), where w is the weight vector.
     */
    fun calculateBeta(x: Array<DoubleArray>, y: F64Array, weights: F64Array): DoubleArray {
        check(weights.size == y.size) { "weights and y have different size" }
        check(x.all { it.size == weights.size }) { "X columns and weights vector have different size" }

        val inverse = calculateBetaVariance(x, weights)
        val Wy = weights * y
        val XTWy = DoubleArray(x.size) { i -> x[i].asF64Array().dot(Wy) }
        return DoubleArray(x.size) { inverse.getRow(it).asF64Array().dot(XTWy) }
    }

    /**
     * Calculates (X^T W X)^{-1}, the variance multiplier of Var beta.
     *
     * Var beta = sigma^2 (X^T W X)^{-1}
     */
    fun calculateBetaVariance(x: Array<DoubleArray>, weights: F64Array): RealMatrix {
        // the resulting matrix is symmetric by construction, let's not waste time computing the upper half
        val XTOIX = Array(x.size) { i ->
            DoubleArray(x.size) { j ->
                if (i < j) 0.0 else x[i].asF64Array().times(weights).dot(x[j].asF64Array())
            }
        }
        // fill the upper half
        x.indices.forEach { i -> x.indices.forEach { j -> if (i < j) x[i][j] = x[j][i] }}
        return LUDecomposition(Array2DRowRealMatrix(XTOIX)).solver.inverse
    }
}