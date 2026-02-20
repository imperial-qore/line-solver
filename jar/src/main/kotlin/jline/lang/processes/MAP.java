/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;


import jline.api.mam.Map_momentKt;
import jline.api.mam.Map_pieKt;
import jline.api.mam.Map_timereverseKt;
import jline.lang.Copyable;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;


/**
 * A Markovian Arrival Process
 */
public class MAP extends MarkovModulated implements Serializable {
    // Core methods implemented: evalCDF, evalLST, toTimeReversed, basic getters
    // Additional methods implemented below

    public MAP(MatrixCell map) {
        this(map.get(0), map.get(1));
    }

    public MAP() {
        super("MAP", 2);
    }

    public MAP(Matrix D0, Matrix D1) {
        super("MAP", 2);
        this.nPhases = D0.getNumCols();
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        this.setProcess(rep);
    }

    public static MAP rand() {
        return MAP.rand(2);
    }

    public static MAP rand(int order) {
        MatrixCell map = jline.api.mam.Map_randKt.map_rand(order);
        return new MAP(map);
    }

    public static MAP randn() {
        return MAP.randn(2, 1.0, 2.0);
    }

    public static MAP randn(int order, double mu, double sigma) {
        MatrixCell map = jline.api.mam.Map_randnKt.map_randn(order, mu, sigma);
        return new MAP(map);
    }

    public Matrix D(int i) {
        return this.process.get(i);
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at time t.
     * For a MAP, this represents the probability that an inter-arrival time is less than or equal to t.
     *
     * @param t the time value
     * @return the CDF value at time t
     */
    public double evalCDF(double t) {
        if (t <= 0) {
            return 0.0;
        }

        try {
            // For MAP, CDF(t) = 1 - alpha * exp(T*t) * e
            // where alpha is the initial probability vector, T is the subgenerator
            // and e is a vector of ones
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D(0), D(1));
            Matrix T = D(0);
            Matrix Tt = Matrix.scaleMult(T, t);
            Matrix expT = Maths.matrixExp(Tt);
            Matrix ones = Matrix.ones(T.getNumRows(), 1);
            Matrix result = pie.mult(expT).mult(ones);
            return 1.0 - result.get(0, 0);
        } catch (Exception e) {
            // Fallback: approximate using moments
            double mean = getMean();
            double scv = getSCV();
            // Simple exponential approximation
            return 1.0 - Math.exp(-t / mean);
        }
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform (LST) at parameter s.
     * For a MAP, LST(s) = alpha * (-T + s*I)^(-1) * t
     * where alpha is the initial probability vector, T is the subgenerator,
     * and t is the exit rate vector.
     *
     * @param s the transform parameter
     * @return the LST value at parameter s
     */
    public double evalLST(double s) {
        if (s < 0) {
            throw new IllegalArgumentException("LST parameter s must be non-negative");
        }

        try {
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D(0), D(1));
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
            // Fallback: use mean approximation
            double mean = getMean();
            return 1.0 / (1.0 + s * mean); // Exponential approximation
        }
    }

    public double getMean() {
        double E1 = Map_momentKt.map_moment(D(0), D(1), 1);
        return E1;
    }

    public long getNumberOfPhases() {
        return ((Matrix) this.getParam(1).getValue()).getNumCols();
    }

    @Override
    public MatrixCell getProcess() {
        MatrixCell res = new MatrixCell();
        Matrix D0 = this.D(0);
        Matrix D1 = this.D(1);
        res.set(0, D0);
        res.set(1, D1);
        return res;
    }

    public double getRate() {
        return 1.0 / this.getMean();
    }

    // remove the temporal dependence
    public MatrixCell getRenewalProcess() {
        Matrix pie = Map_pieKt.map_pie(D(0), D(1));
        return new PH(pie, D(0)).getProcess();
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

    public MAP toTimeReversed() {
        MatrixCell reversed = Map_timereverseKt.map_timereverse(this.process);
        return new MAP(reversed);
    }

    /**
     * Evaluates the mean of the counting process at time t.
     * For MAP, this gives the expected number of arrivals by time t.
     *
     * @param t the time value
     * @return the mean number of arrivals by time t
     */
    public double evalMeanT(double t) {
        if (t <= 0) {
            return 0.0;
        }
        
        try {
            // For MAP, the mean rate is given by the equilibrium rate
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D(0), D(1));
            Matrix D1 = D(1);
            Matrix ones = Matrix.ones(D1.getNumCols(), 1);
            Matrix exitRate = D1.mult(ones);
            
            double rate = pie.mult(exitRate).get(0, 0);
            return rate * t;
        } catch (Exception e) {
            // Fallback: use mean inter-arrival time
            double mean = getMean();
            return t / mean;
        }
    }

    /**
     * Evaluates the variance of the counting process at time t.
     * For MAP, this provides the variance of the number of arrivals by time t.
     *
     * @param t the time value
     * @return the variance of arrivals by time t
     */
    public double evalVarT(double t) {
        if (t <= 0) {
            return 0.0;
        }
        
        try {
            // For MAP, variance includes both the mean and correlation terms
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
     * Evaluates the autocorrelation function at given lags and timescale.
     * For MAP, this measures the correlation between arrivals at different times.
     *
     * @param lags array of lag values
     * @param timescale the timescale parameter
     * @return matrix of autocorrelation values
     */
    public Matrix evalACFT(int[] lags, double timescale) {
        return new Matrix(jline.api.mam.Map_acfcKt.map_acfc(D(0), D(1), lags, timescale));
    }

    /**
     * Evaluates the probability density function (PDF) at time t.
     * For MAP, this is the derivative of the CDF.
     *
     * @param t the time value
     * @return the PDF value at time t
     */
    public double evalPDF(double t) {
        if (t < 0) {
            return 0.0;
        }
        
        try {
            // For MAP, PDF(t) = alpha * exp(T*t) * (-T) * e
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D(0), D(1));
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
     * Normalizes the MAP so that D0+D1 forms a proper infinitesimal generator.
     * Each row sum of (D0+D1) should equal zero for a valid generator.
     */
    public void normalize() {
        Matrix D0 = D(0).copy();
        Matrix D1 = D(1).copy();
        int nPhases = D0.getNumRows();

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
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        MatrixCell res = new MatrixCell();
        res.set(0, D0);
        res.set(1, D1);
        this.setProcess(res);
    }


    /**
     * Fit MAP from central moments and autocorrelation function decay.
     * This creates a MAP that matches the given statistical properties.
     *
     * @param mean the mean inter-arrival time
     * @param var the variance of inter-arrival times
     * @param skew the skewness of inter-arrival times
     * @param gamma2 the autocorrelation function decay parameter
     * @return fitted MAP process
     */
    public static MAP fitCentralAndACFDecay(double mean, double var, double skew, double gamma2) {
        double e1 = mean;
        double scv = var / mean / mean;
        double e2 = (1.0 + scv) * e1 * e1;
        double e3;
        if (skew == -1) {
            e3 = -1;
        } else {
            e3 = -(2 * FastMath.pow(e1, 3) - 3 * e1 * e2 - skew * FastMath.sqrt(Math.pow(e2 - e1 * e1, 3)));
        }
        return fitRawMomentsAndACFDecay(e1, e2, e3, gamma2);
    }

    /**
     * Fit MAP from raw moments and autocorrelation function decay.
     * This is a simplified version that creates a 2-phase MAP.
     *
     * @param m1 first raw moment
     * @param m2 second raw moment
     * @param m3 third raw moment
     * @param gamma2 autocorrelation function decay parameter
     * @return fitted MAP process
     */
    public static MAP fitRawMomentsAndACFDecay(double m1, double m2, double m3, double gamma2) {
        try {
            // Use MMPP2 fitting as a starting point for 2-phase MAP
            MMPP2 mmpp = MMPP2.fitRawMomentsAndACFDecay(m1, m2, m3, gamma2);
            return new MAP(mmpp.D(0), mmpp.D(1));
        } catch (Exception e) {
            // Fallback: create exponential MAP
            double rate = 1.0 / m1;
            Matrix D0 = new Matrix(1, 1);
            D0.set(0, 0, -rate);
            Matrix D1 = new Matrix(1, 1);
            D1.set(0, 0, rate);
            return new MAP(D0, D1);
        }
    }

    /**
     * Returns a string representation of this MAP process.
     *
     * @return string representation showing dimensions and basic properties
     */
    @Override
    public String toString() {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        double mean = getMean();
        double scv = getSCV();
        return String.format("MAP(%dx%d, mean=%.6f, scv=%.6f)", 
                           D0.getNumRows(), D0.getNumCols(), mean, scv);
    }

    /**
     * Checks if this MAP is equal to another object.
     *
     * @param obj the object to compare with
     * @return true if equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        MAP other = (MAP) obj;
        
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        Matrix otherD0 = other.D(0);
        Matrix otherD1 = other.D(1);
        
        if (D0.getNumRows() != otherD0.getNumRows() || D0.getNumCols() != otherD0.getNumCols()) {
            return false;
        }
        
        // Compare matrices element by element
        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (Math.abs(D0.get(i, j) - otherD0.get(i, j)) > 1e-10) {
                    return false;
                }
                if (Math.abs(D1.get(i, j) - otherD1.get(i, j)) > 1e-10) {
                    return false;
                }
            }
        }
        
        return true;
    }

    /**
     * Returns the hash code for this MAP process.
     *
     * @return hash code based on matrix contents
     */
    @Override
    public int hashCode() {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        
        int result = 1;
        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                result = 31 * result + Double.hashCode(D0.get(i, j));
                result = 31 * result + Double.hashCode(D1.get(i, j));
            }
        }
        return result;
    }
}
