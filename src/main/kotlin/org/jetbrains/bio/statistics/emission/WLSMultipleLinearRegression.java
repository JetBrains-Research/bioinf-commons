package org.jetbrains.bio.statistics.emission;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.stat.regression.AbstractMultipleLinearRegression;
import org.apache.commons.math3.exception.DimensionMismatchException;

import org.apache.commons.math3.linear.*;

public class WLSMultipleLinearRegression extends AbstractMultipleLinearRegression {
    /** Entries of the matrix. */
    private DesignMatrix X;
    //use DiagonalMatrix instead of RealMatrix
    private DiagonalMatrix Omega;

    public WLSMultipleLinearRegression() {
    }

    public RealMatrix getX() {
        return X;
    }

    public void newSampleData(double[] y, double[][] x, double[] weights) { //weight now is double[] instead of double[][]
        this.newXSampleData(x);
        this.validateSampleData(y);
        this.newYSampleData(y);
        this.validateCovarianceData(weights);
        this.newCovarianceData(weights);
    }

    public void newSampleData(double[] y, double[] weights) {
        if (this.X == null) {
            throw new NullArgumentException();
        }
        this.validateSampleData(y);
        this.newYSampleData(y);
        this.validateCovarianceData(weights);
        this.newCovarianceData(weights);
    }
    //override
    protected void newXSampleData(double[][] x) {
        if (x == null) {
            throw new NullArgumentException();
        }
        if (x.length == 0) {
            throw new NoDataException();
        }
        X = new DesignMatrix(x);
    }

    protected void validateCovarianceData(double[] weights) {
        if (this.X.getRowDimension() != weights.length) {
            throw new DimensionMismatchException(this.X.getRowDimension(), weights.length);
        }
    }

    protected void validateSampleData(double[] y){
        if ((this.X == null) || (y == null)) {
            throw new NullArgumentException();
        }
        if (this.X.getRowDimension() != y.length) {
            throw new DimensionMismatchException(y.length, this.X.getRowDimension());
        }
        if (this.X.getRowDimension() == 0) {  // Must be no y data either
            throw new NoDataException();
        }
        if (this.X.getColumnDimension() > this.X.getRowDimension()) {
            throw new MathIllegalArgumentException(
                    LocalizedFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS,
                    this.X.getRowDimension(), this.X.getColumnDimension());
        }
    }
    //use DiagonalMatrix
    protected void newCovarianceData(double[] omega) {
        this.Omega = new DiagonalMatrix(omega);
    }

    public RealVector calculateBeta() {
        RealMatrix XT = this.X.transpose();
        RealMatrix XTOTX = XT.multiply(Omega.multiply(this.X));
        RealMatrix inverse = (new LUDecomposition(XTOTX)).getSolver().getInverse();
        return inverse.multiply(XT).operate((Omega).operate(this.getY()));
    }

    protected RealMatrix calculateBetaVariance() {
        RealMatrix XTOX = this.X.transpose().multiply(Omega.multiply(this.X));
        return (new LUDecomposition(XTOX)).getSolver().getInverse();
    }

    protected double calculateErrorVariance() {
        RealVector residuals = this.calculateResiduals();
        double t = residuals.dotProduct(this.Omega.operate(residuals));
        return t / (double)(this.X.getRowDimension() - this.X.getColumnDimension());
    }
}


