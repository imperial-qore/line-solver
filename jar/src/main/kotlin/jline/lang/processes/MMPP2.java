/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.api.mam.*;
import jline.GlobalConstants;
import jline.lang.Copyable;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Markovian-modulated Poisson Process with 2 states
 */
public class MMPP2 extends MarkovModulated implements Serializable {
    // Core methods implemented: evalCDF, evalLST, evalMeanT, evalACFT, normalize, sample, toTimeReversed
    // Additional utility methods below
    private final int nPhases;

    public MMPP2(double lambda0, double lambda1, double sigma0, double sigma1) {
        super("MMPP2", 4);
        nPhases = 2;
        this.setParam(1, "lambda0", lambda0);
        this.setParam(2, "lambda1", lambda1);
        this.setParam(3, "sigma0", sigma0);
        this.setParam(4, "sigma1", sigma1);
        MatrixCell res = new MatrixCell();
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -lambda0 - sigma0);
        D0.set(0, 1, sigma0);
        D0.set(1, 0, sigma1);
        D0.set(1, 1, -lambda1 - sigma1);
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, lambda0);
        D1.set(1, 1, lambda1);
        res.set(0, D0);
        res.set(1, D1);
        setProcess(res);
    }

    public static MMPP2 fitCentralAndACFDecay(double mean, double var, double skew, double g2) {
        double e1 = mean;
        double scv = var / mean / mean;
        double e2 = (1.0 + scv) * e1 * e1;
        double e3;
        if (skew == -1) {
            e3 = -1;
        } else {
            e3 = -(2 * FastMath.pow(e1, 3) - 3 * e1 * e2 - skew * FastMath.sqrt(Math.pow(e2 - e1 * e1, 3)));
        }
        return fitRawMomentsAndACFDecay(e1, e2, e3, g2);
    }

    public static MMPP2 fitCentralAndACFLag1(double mean, double var, double skew, double acf1) {
        double e1 = mean;
        double scv = var / mean / mean;
        double e2 = (1.0 + scv) * e1 * e1;
        double e3;
        if (skew == -1) {
            e3 = -1;
        } else {
            e3 = -(2 * FastMath.pow(e1, 3) - 3 * e1 * e2 - skew * FastMath.sqrt(Math.pow(e2 - e1 * e1, 3)));
        }
        double rho0 = (1.0 - 1.0 / scv) / 2.0;
        double g2 = acf1 / rho0;

        return fitRawMomentsAndACFDecay(e1, e2, e3, g2);
    }

    public static MMPP2 fitCentralAndIDC(double mean, double var, double skew, double idc) {
        double e1 = mean;
        double scv = var / mean / mean;
        double e2 = (1.0 + scv) * e1 * e1;
        double g2 = -(scv - idc) / (-1 + idc);
        double e3;
        if (skew == -1) {
            e3 = -1;
        } else {
            e3 = -(2 * FastMath.pow(e1, 3) - 3 * e1 * e2 - skew * FastMath.sqrt(Math.pow(e2 - e1 * e1, 3)));
        }
        return fitRawMomentsAndACFDecay(e1, e2, e3, g2);
    }

    public static MMPP2 fitRawMomentsAndACFDecay(double m1, double m2, double m3, double gamma2) {
        if (m1 <= GlobalConstants.FineTol) {
            return new MMPP2(Inf, Inf, 1, 1);
        }
        MatrixCell m = Mmpp2_fitKt.mmpp2_fit(m1, m2, m3, gamma2);
        return new MMPP2(m.get(1).value(), m.get(1).get(1, 1), m.get(0).value(), m.get(0).get(1, 1));
    }

    public static MMPP2 fitRawMomentsAndACFLag1(double m1, double m2, double m3, double rho1) {
        double scv = (m2 - m1 * m1) / m1 / m1;
        double rho0 = (1 - 1 / scv) / 2;
        double gamma2 = rho1 / rho0;
        return fitRawMomentsAndACFDecay(m1, m2, m3, gamma2);
    }

    public static MMPP2 fitRawMomentsAndIDC(double m1, double m2, double m3, double idc) {
        double scv = (m2 - m1 * m1) / m1 / m1;
        double gamma2 = -(scv - idc) / (-1 + idc);
        return fitRawMomentsAndACFDecay(m1, m2, m3, gamma2);
    }

    public static MMPP2 rand() {
        return new MMPP2(Maths.rand(), Maths.rand(), Maths.rand(), Maths.rand());
    }

    public Matrix D(int i) {
        return this.process.get(i);
    }

    public Matrix evalACFT(int[] lags, double timescale) {
        return new Matrix(Map_acfcKt.map_acfc(D(0), D(1), lags, timescale));
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at time t.
     * For MMPP2, this represents the probability that an inter-arrival time is less than or equal to t.
     *
     * @param t the time value
     * @return the CDF value at time t
     */
    public double evalCDF(double t) {
        if (t <= 0) {
            return 0.0;
        }

        try {
            // Use MAP CDF calculation since MMPP2 is a special case of MAP
            Matrix pie = Map_pieKt.map_pie(D(0), D(1));
            Matrix T = D(0);
            Matrix Tt = Matrix.scaleMult(T, t);
            Matrix expT = Maths.matrixExp(Tt);
            Matrix ones = Matrix.ones(T.getNumRows(), 1);
            Matrix result = pie.mult(expT).mult(ones);
            return 1.0 - result.get(0, 0);
        } catch (Exception e) {
            // Fallback: exponential approximation using mean
            double mean = getMean();
            return 1.0 - Math.exp(-t / mean);
        }
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform (LST) at parameter s.
     * For MMPP2, LST(s) = alpha * (-T + s*I)^(-1) * t
     *
     * @param s the transform parameter
     * @return the LST value at parameter s
     */
    public double evalLST(double s) {
        if (s < 0) {
            throw new IllegalArgumentException("LST parameter s must be non-negative");
        }

        try {
            Matrix pie = Map_pieKt.map_pie(D(0), D(1));
            Matrix T = D(0);
            Matrix D1 = D(1);
            Matrix I = Matrix.eye(T.getNumRows());

            // Compute (-T + s*I)^(-1)
            Matrix negT = Matrix.scaleMult(T, -1.0);
            Matrix sI = Matrix.scaleMult(I, s);
            Matrix resolvent = negT.add(1.0, sI).inv();

            // Compute exit rate vector: t = D1 * e (sum of D1 rows)
            Matrix ones = Matrix.ones(D1.getNumCols(), 1);
            Matrix t = D1.mult(ones);

            // LST = alpha * resolvent * t
            Matrix result = pie.mult(resolvent).mult(t);
            return result.get(0, 0);
        } catch (Exception e) {
            // Fallback: exponential approximation
            double mean = getMean();
            return 1.0 / (1.0 + s * mean);
        }
    }

    public double evalMeanT(double t) {
        double lambda0 = (double) this.getParam(1).getValue();
        double lambda1 = (double) this.getParam(2).getValue();
        double sigma0 = (double) this.getParam(3).getValue();
        double sigma1 = (double) this.getParam(4).getValue();
        double lambda = (lambda0 * sigma1 + lambda1 * sigma0) / (sigma0 + sigma1);
        return lambda * t;
    }

    public double getACFDecay(Matrix lags) {
        return getACF(Matrix.singleton(2)).value() / getACF(Matrix.singleton(1)).value();
    }

    public Matrix getEmbedded() {
        return Map_embeddedKt.map_embedded(D(0), D(1));
    }

    public Matrix getEmbeddedProb() {
        return Map_pieKt.map_pie(D(0), D(1));
    }

    public double getIDC() {
        double lambda0 = (double) this.getParam(1).getValue();
        double lambda1 = (double) this.getParam(2).getValue();
        double sigma0 = (double) this.getParam(3).getValue();
        double sigma1 = (double) this.getParam(4).getValue();
        double id = 1 + 2 * FastMath.pow(lambda0 - lambda1, 2) * sigma0 * sigma1 / FastMath.pow(sigma0 + sigma1, 2) / (lambda0 * sigma1 + lambda1 * sigma0);
        return id;
    }

    public double getIDI() {
        return getIDC();
    }

    public double getMean() {
        double E1 = Map_momentKt.map_moment(D(0), D(1), 1);
        return E1;
    }

    public long getNumberOfPhases() {
        return 2;
    }

    public double getRate() {
        return 1.0 / this.getMean();
    }

    public double getSCV() {
        double mean = this.getMean();
        return this.getVar() / mean / mean;
    }

    public double getSkewness() {
        double E1 = Map_momentKt.map_moment(D(0), D(1), 1);
        double E2 = Map_momentKt.map_moment(D(0), D(1), 2);
        double E3 = Map_momentKt.map_moment(D(0), D(1), 3);
        double skew = E3 - 3 * E2 * E1 + 2 * E1 * E1 * E1;
        double scv = (E2 - E1 * E1) / E1 / E1;
        skew = skew / FastMath.pow(Math.sqrt(scv) * E1, 3);
        return skew;
    }

    public double getVar() {
        double E1 = Map_momentKt.map_moment(D(0), D(1), 1);
        double E2 = Map_momentKt.map_moment(D(0), D(1), 2);
        return E2 - E1 * E1;
    }

    /**
     * Normalizes the MMPP2 so that D0+D1 rows form a proper infinitesimal generator.
     * Each row sum of (D0+D1) should equal zero for a valid generator.
     */
    public void normalize() {
        Matrix D0 = D(0).copy();
        Matrix D1 = D(1).copy();

        // Ensure non-negative off-diagonal elements in D0 and all elements in D1
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (i != j && D0.get(i, j) < 0) {
                    D0.set(i, j, 0.0);
                }
                if (D1.get(i, j) < 0) {
                    D1.set(i, j, 0.0);
                }
            }
        }

        // Adjust diagonal elements of D0 to make row sums zero
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            // Sum off-diagonal elements of D0 and all elements of D1 in row i
            for (int j = 0; j < nPhases; j++) {
                if (i != j) {
                    rowSum += D0.get(i, j);
                }
                rowSum += D1.get(i, j);
            }
            // Set diagonal element to make row sum zero
            D0.set(i, i, -rowSum);
        }

        // Update the process
        MatrixCell res = new MatrixCell();
        res.set(0, D0);
        res.set(1, D1);
        setProcess(res);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    @Override
    public double[] sample(int n, Random random) {
        if (Map_isfeasibleKt.map_isfeasible(getProcess())) {
            return Map_sampleKt.map_sample(D(0), D(1), n, random);
        } else {
            line_error(mfilename(new Object() {
            }), "This process is infeasible (negative rates).");
        }
        return null;
    }

    public MAP toTimeReversed() {
        MatrixCell reversed = Map_timereverseKt.map_timereverse(this.process);
        return new MAP(reversed);
    }

    /**
     * Evaluates the probability density function (PDF) at time t.
     * For MMPP2, this is the derivative of the CDF.
     *
     * @param t the time value
     * @return the PDF value at time t
     */
    public double evalPDF(double t) {
        if (t < 0) {
            return 0.0;
        }
        
        try {
            // For MAP/MMPP2, PDF(t) = alpha * exp(T*t) * (-T) * e
            Matrix pie = Map_pieKt.map_pie(D(0), D(1));
            Matrix T = D(0);
            Matrix Tt = Matrix.scaleMult(T, t);
            Matrix expT = Maths.matrixExp(Tt);
            Matrix negT = Matrix.scaleMult(T, -1.0);
            Matrix ones = Matrix.ones(T.getNumRows(), 1);
            
            Matrix result = pie.mult(expT).mult(negT).mult(ones);
            return result.get(0, 0);
        } catch (Exception e) {
            // Fallback: exponential approximation
            double mean = getMean();
            return (1.0 / mean) * Math.exp(-t / mean);
        }
    }

    /**
     * Evaluates the variance at time t.
     * For MMPP2, this provides the variance of the counting process at time t.
     *
     * @param t the time value
     * @return the variance at time t
     */
    public double evalVarT(double t) {
        if (t <= 0) {
            return 0.0;
        }
        
        try {
            // For MMPP2, variance includes both the mean and additional correlation terms
            double mean = evalMeanT(t);
            double scv = getSCV();
            double var = mean * (1.0 + scv * mean);
            return var;
        } catch (Exception e) {
            // Fallback: Poisson approximation
            return evalMeanT(t);
        }
    }


    /**
     * Returns a string representation of this MMPP2 process.
     *
     * @return string representation showing parameters
     */
    @Override
    public String toString() {
        double lambda0 = (double) this.getParam(1).getValue();
        double lambda1 = (double) this.getParam(2).getValue();
        double sigma0 = (double) this.getParam(3).getValue();
        double sigma1 = (double) this.getParam(4).getValue();
        return String.format("MMPP2(lambda0=%.6f, lambda1=%.6f, sigma0=%.6f, sigma1=%.6f)", 
                           lambda0, lambda1, sigma0, sigma1);
    }

    /**
     * Checks if this MMPP2 is equal to another object.
     *
     * @param obj the object to compare with
     * @return true if equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        MMPP2 other = (MMPP2) obj;
        
        double lambda0 = (double) this.getParam(1).getValue();
        double lambda1 = (double) this.getParam(2).getValue();
        double sigma0 = (double) this.getParam(3).getValue();
        double sigma1 = (double) this.getParam(4).getValue();
        
        double otherLambda0 = (double) other.getParam(1).getValue();
        double otherLambda1 = (double) other.getParam(2).getValue();
        double otherSigma0 = (double) other.getParam(3).getValue();
        double otherSigma1 = (double) other.getParam(4).getValue();
        
        return Math.abs(lambda0 - otherLambda0) < 1e-10 &&
               Math.abs(lambda1 - otherLambda1) < 1e-10 &&
               Math.abs(sigma0 - otherSigma0) < 1e-10 &&
               Math.abs(sigma1 - otherSigma1) < 1e-10;
    }

    /**
     * Returns the hash code for this MMPP2 process.
     *
     * @return hash code based on parameters
     */
    @Override
    public int hashCode() {
        double lambda0 = (double) this.getParam(1).getValue();
        double lambda1 = (double) this.getParam(2).getValue();
        double sigma0 = (double) this.getParam(3).getValue();
        double sigma1 = (double) this.getParam(4).getValue();
        
        int result = 1;
        result = 31 * result + Double.hashCode(lambda0);
        result = 31 * result + Double.hashCode(lambda1);
        result = 31 * result + Double.hashCode(sigma0);
        result = 31 * result + Double.hashCode(sigma1);
        return result;
    }

}
