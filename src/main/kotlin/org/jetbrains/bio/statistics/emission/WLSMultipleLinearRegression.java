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
    private DesignMatrix xMatrix;
    private RealVector yVector;
    //use DiagonalMatrix instead of RealMatrix
    private DiagonalMatrix Omega;

    public WLSMultipleLinearRegression() {
    }

    public RealMatrix getX() {
        return xMatrix;
    }

    public void newSampleData(double[] y, double[][] x, double[] weights) { //weight now is double[] instead of double[][]
        this.newXSampleData(x);
        this.validateSampleData(y);
        this.newYSampleData(y);
        this.validateCovarianceData(weights);
        this.newCovarianceData(weights);
    }

    public void newSampleData(double[] y, double[] weights) {
        if (this.xMatrix == null) {
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

        xMatrix = new DesignMatrix(x);
    }
    //override
    protected void newYSampleData(double[] y) {
        if (y == null) {
            throw new NullArgumentException();
        }
        if (y.length == 0) {
            throw new NoDataException();
        }
        this.yVector = new ArrayRealVector(y);
    }

    protected void validateCovarianceData(double[] weights) {
        if (this.xMatrix.getRowDimension() != weights.length) {
            throw new DimensionMismatchException(this.xMatrix.getRowDimension(), weights.length);
        }
    }

    protected void validateSampleData(double[] y){
        if ((this.xMatrix == null) || (y == null)) {
            throw new NullArgumentException();
        }
        if (this.xMatrix.getRowDimension() != y.length) {
            throw new DimensionMismatchException(y.length, this.xMatrix.getRowDimension());
        }
        if (this.xMatrix.getRowDimension() == 0) {  // Must be no y data either
            throw new NoDataException();
        }
        if (this.xMatrix.getColumnDimension() > this.xMatrix.getRowDimension()) {
            throw new MathIllegalArgumentException(
                    LocalizedFormats.NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS,
                    this.xMatrix.getRowDimension(), this.xMatrix.getColumnDimension());
        }
    }
    //use DiagonalMatrix
    protected void newCovarianceData(double[] omega) {
        this.Omega = new DiagonalMatrix(omega);
    }

    public RealVector calculateBeta() {
        RealMatrix XT = xMatrix.transpose();
        RealMatrix XTOTX = XT.multiply(Omega.multiply(xMatrix));
        RealMatrix inverse = (new LUDecomposition(XTOTX)).getSolver().getInverse();
        return inverse.multiply(XT).operate((Omega).operate(yVector));
    }

    protected RealMatrix calculateBetaVariance() {
        RealMatrix XTOX = xMatrix.transpose().multiply(Omega.multiply(xMatrix));
        return (new LUDecomposition(XTOX)).getSolver().getInverse();
    }
}


