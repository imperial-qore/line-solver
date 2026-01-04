/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.List;

import static jline.lib.butools.APHFrom3MomentsKt.APHFrom3Moments;

/**
 * An acyclic phase type distribution
 */
public class APH extends Markovian {

    public APH(List<Double> alpha, Matrix T) {
        this(new Matrix(alpha).transpose(), T);
    }

    public APH(Matrix alpha, Matrix T) {
        super("APH", 3);

        nPhases = alpha.getNumElements();

        this.setParam(1, "n", nPhases);
        this.setParam(2, "alpha", alpha);
        this.setParam(3, "T", T);

        MatrixCell rep = new MatrixCell();
        rep.set(0, T);
        rep.set(1, T.mult(alpha.repmat(nPhases, 1)).scale(-1));
        this.setProcess(rep);
    }

    public static APH fit(double mean, double scv, double skew) {
        if (mean <= GlobalConstants.FineTol) {
            Matrix alpha = new Matrix(1, 1, 1);
            alpha.set(0, 0, 1.0);
            return new APH(alpha, Matrix.singleton(Inf));
        } else {
            double e1 = mean;
            double e2 = (1.0 + scv) * e1 * e1;
            double e3 = -(2.0 * FastMath.pow(e1, 3.0) - 3.0 * e1 * e2 - skew * FastMath.pow(e2 - e1 * e1, 3.0 / 2.0));
            return APHFrom3Moments(new double[]{e1, e2, e3}, 100, GlobalConstants.Zero);
        }
    }

    public static APH fitCentral(double mean, double var, double skew) {
        // Fit the distribution from first three central moments (mean, variance, skewness)
        APH ex;
        if (mean <= GlobalConstants.FineTol) {
            Exp myexp = new Exp(Inf);
            ex = APH.fit(mean, myexp.getSCV(), myexp.getSkewness());
        } else {
            ex = APH.fit(mean, var / mean / mean, skew);
            ex.immediate = false;
        }
        return ex;
    }

    public static APH fitMeanAndSCV(double mean, double scv) {
        APH ex;
        if (mean <= GlobalConstants.FineTol) {
            Matrix T = new Matrix(1, 1);
            T.set(0, 0, -1.0 / mean);
            Matrix alpha = Matrix.singleton(1.0);
            ex = new APH(alpha, T);
        } else {
            double e1 = mean;
            double e2 = (1 + scv) * e1 * e1;
            double cv2 = e2 / e1 / e1 - 1.0;
            double lambda = 1.0 / e1;
            int N = FastMath.max((int) FastMath.ceil(1.0 / cv2), 2);
            double p = 1.0 / (cv2 + 1.0 + (cv2 - 1.0) / (N - 1));
            Matrix T = new Matrix(N, N);
            Matrix.eye(N).scaleEq(-lambda * p * N, T);
            for (int i = 0; i < N - 1; i++) {
                T.set(i, i + 1, -1.0 * T.get(i, i));
            }
            T.set(N - 1, N - 1, -lambda * N);
            Matrix alpha = new Matrix(1, N);
            alpha.set(0, 0, p);
            alpha.set(0, N - 1, 1.0 - p);
            ex = new APH(alpha, T);
            ex.immediate = false;
        }
        return ex;
    }

    /**
     * Fits an APH distribution from the first three raw moments.
     *
     * @param m1 First raw moment (mean)
     * @param m2 Second raw moment
     * @param m3 Third raw moment
     * @return APH distribution fitted to the given moments
     */
    public static APH fitRawMoments(double m1, double m2, double m3) {
        if (m1 <= GlobalConstants.FineTol) {
            // Return APH with very high rate (essentially immediate)
            Matrix alpha = Matrix.singleton(1.0);
            Matrix T = Matrix.singleton(-Inf);
            return new APH(alpha, T);
        }
        return APHFrom3Moments(new double[]{m1, m2, m3}, 100, GlobalConstants.Zero);
    }

    @Override
    public double evalLST(double s) {
        return super.evalLST(s);
    }

    public Matrix getInitProb() {
        return ((Matrix) this.getParam(2).getValue()).copy();
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================

    /**
     * Kotlin-style property alias for getInitProb()
     */
    public Matrix initProb() {
        return getInitProb();
    }

    // =================== CDF EVALUATION METHODS ===================

    /**
     * Evaluates the CDF at default time points and returns as a Matrix.
     * Returns a matrix with [CDF_value, time] pairs (matching MATLAB format).
     *
     * @return Matrix with numPoints rows and 2 columns [CDF, time]
     */
    public Matrix evalCDFMatrix() {
        Matrix alpha = getInitProb();
        Matrix T = (Matrix) getParam(3).getValue();
        int order = alpha.getNumCols();
        Matrix e = Matrix.ones(order, 1);

        double mean = getMean();
        double var = getVar();
        double sigma = FastMath.sqrt(var);

        int numPoints = 500;
        double maxT = mean + 10 * sigma;
        Matrix result = new Matrix(numPoints, 2);

        for (int i = 0; i < numPoints; i++) {
            double t = (numPoints > 1) ? i * maxT / (numPoints - 1) : 0;
            // F(t) = 1 - alpha * expm(T*t) * e
            Matrix Tt = T.scale(t);
            Matrix expTt = Tt.expm_higham();
            double cdfValue = 1.0 - alpha.mult(expTt).mult(e).get(0, 0);
            result.set(i, 0, cdfValue);  // CDF value
            result.set(i, 1, t);          // time
        }

        return result;
    }

    /**
     * Evaluates the CDF at specified time points and returns as a Matrix.
     *
     * @param timePoints Array of time points at which to evaluate the CDF
     * @return Matrix with timePoints.length rows and 2 columns [CDF, time]
     */
    public Matrix evalCDFMatrix(double[] timePoints) {
        Matrix alpha = getInitProb();
        Matrix T = (Matrix) getParam(3).getValue();
        int order = alpha.getNumCols();
        Matrix e = Matrix.ones(order, 1);

        Matrix result = new Matrix(timePoints.length, 2);
        for (int i = 0; i < timePoints.length; i++) {
            double t = timePoints[i];
            // F(t) = 1 - alpha * expm(T*t) * e
            Matrix Tt = T.scale(t);
            Matrix expTt = Tt.expm_higham();
            double cdfValue = 1.0 - alpha.mult(expTt).mult(e).get(0, 0);
            result.set(i, 0, cdfValue);  // CDF value
            result.set(i, 1, t);          // time
        }

        return result;
    }

}
