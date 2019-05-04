package org.jetbrains.bio.statistics.emission;

import org.apache.commons.math3.stat.regression.AbstractMultipleLinearRegression;
import org.apache.commons.math3.exception.DimensionMismatchException;

import org.apache.commons.math3.linear.*;

public class WLSMultipleLinearRegression extends AbstractMultipleLinearRegression {
    //use DiagonalMatrix instead of RealMatrix
    private DiagonalMatrix Omega;

    public WLSMultipleLinearRegression() {
    }

    public RealMatrix getx() {
        return getX();
    }

    public void newSampleData(double[] y, double[][] x, double[] weights) { //weight now is double[] instead of double[][]
        this.validateSampleData(x, y);
        this.newYSampleData(y);
        this.newXSampleData(x);
        this.validateCovarianceData(x, weights);
        this.newCovarianceData(weights);
    }
    //override
    protected void validateCovarianceData(double[][] x, double[] weights) {
        if (x.length != weights.length) {
            throw new DimensionMismatchException(x.length, weights.length);
        }
    }

    //use DiagonalMatrix
    protected void newCovarianceData(double[] omega) {
        this.Omega = new DiagonalMatrix(omega);
    }

    public RealVector calculateBeta() {
        RealMatrix XT = this.getX().transpose();
        RealMatrix XTOTX = XT.multiply(Omega.multiply(this.getX()));
        RealMatrix inverse = (new LUDecomposition(XTOTX)).getSolver().getInverse();
        return inverse.multiply(XT).operate((Omega).operate(this.getY()));
    }

    protected RealMatrix calculateBetaVariance() {
        RealMatrix XTOX = this.getX().transpose().multiply(Omega.multiply(this.getX()));
        return (new LUDecomposition(XTOX)).getSolver().getInverse();
    }

    protected double calculateErrorVariance() {
        RealVector residuals = this.calculateResiduals();
        double t = residuals.dotProduct(this.Omega.operate(residuals));
        return t / (double)(this.getX().getRowDimension() - this.getX().getColumnDimension());
    }
}


