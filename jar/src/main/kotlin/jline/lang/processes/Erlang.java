/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import jline.util.RandomManager;
import java.util.Random;

import static jline.api.mam.Map_sampleKt.map_sample;

/**
 * An Erlang-n distribution with n phases.
 */

public class Erlang extends Markovian implements Serializable {

    /**
     * Creates an Erlang distribution with the specified phase rate and number of phases.
     * 
     * @param phaseRate the rate parameter for each phase
     * @param n the number of phases
     */
    public Erlang(double phaseRate, int n) {
        super("Erlang", 2);
        nPhases = n;
        this.setParam(1, "alpha", phaseRate);
        this.setParam(2, "r", nPhases);
        double alpha = (double) this.getParam(1).getValue();
        double mu = (double) this.getParam(1).getValue();
        int r = (int) this.getParam(2).getValue();
        int size = r;
        Matrix D0 = new Matrix(size, size);
        Matrix D1 = new Matrix(size, size);

        for (int i = 0; i < size - 1; i++) {
            D0.set(i, i, -mu);
            D0.set(i, i + 1, mu);
        }
        D0.set(size - 1, size - 1, -mu);
        D1.set(size - 1, 0, mu);

        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        setProcess(rep);
    }

    /**
     * Creates an Erlang distribution fitted to the specified mean and number of phases.
     * 
     * @param mean the desired mean value
     * @param numPhases the number of phases
     * @return an Erlang distribution with the specified characteristics
     */
    public static Erlang fitMeanAndOrder(double mean, long numPhases) {
        double SCV = 1.0 / numPhases;
        return fitMeanAndSCV(mean, SCV);
    }

    /**
     * Creates an Erlang distribution fitted to the specified mean and squared coefficient of variation.
     * 
     * @param mean the desired mean value
     * @param SCV the desired squared coefficient of variation
     * @return an Erlang distribution with the specified characteristics
     */
    public static Erlang fitMeanAndSCV(double mean, double SCV) {
        int r = (int) FastMath.ceil(1.0 / SCV);
        double alpha = r / mean;
        Erlang er = new Erlang(alpha, r);
        er.immediate = mean < GlobalConstants.FineTol;
        return er;
    }

    /**
     * Creates an Erlang distribution fitted to the specified mean and standard deviation.
     * 
     * @param mean the desired mean value
     * @param stdDev the desired standard deviation
     * @return an Erlang distribution with the specified characteristics
     */
    public static Erlang fitMeanAndStdDev(double mean, double stdDev) {
        return Erlang.fitMeanAndSCV(mean, (mean / Math.pow(stdDev, 2)));
    }

    /**
     * Evaluates the cumulative distribution function at the given point.
     * 
     * @param t the point at which to evaluate the CDF
     * @return the CDF value at point t
     */
    public double evalCDF(double t) {
        double alpha = (double) this.getParam(1).getValue();
        int r = (int) this.getParam(2).getValue();
        double ft = 1;

        for (int j = 0; j < r; j++) {
            int fac_j = 1;
            for (int k = 2; k <= j; k++) {
                fac_j *= k;
            }
            ft -= FastMath.exp(-alpha * t) * (alpha * t) * j / fac_j;
        }

        return ft;
    }

    /**
     * Evaluates the Laplace-Stieltjes transform at the given point.
     * 
     * @param s the transform variable
     * @return the LST value (α/(α+s))^r
     */
    public double evalLST(double s) {
        double alpha = (double) this.getParam(1).getValue();
        int r = (int) this.getParam(2).getValue();
        return FastMath.pow(alpha / (alpha + s), r);
    }

    /**
     * Gets the mean of this Erlang distribution.
     * 
     * @return the mean value (r/α)
     */
    public double getMean() {
        double alpha = (double) this.getParam(1).getValue();
        int r = (int) this.getParam(2).getValue();
        return r / alpha;
    }

    /**
     * Gets the rate of this Erlang distribution.
     * 
     * @return the rate (α/r)
     */
    public double getRate() {
        return 1.0 / getMean();
    }

    /**
     * Gets the squared coefficient of variation.
     * 
     * @return the SCV (1/r)
     */
    public double getSCV() {
        int r = (int) this.getParam(2).getValue();
        return 1.0 / r;
    }

    /**
     * Gets the skewness of this Erlang distribution.
     * 
     * @return the skewness (2/√r)
     */
    public double getSkewness() {
        int r = (int) this.getParam(2).getValue();
        return 2.0 / Math.sqrt(r);
    }

    /**
     * Gets the variance of this Erlang distribution.
     * 
     * @return the variance (r/α²)
     */
    public double getVar() {
        double alpha = (double) this.getParam(1).getValue();
        int r = (int) this.getParam(2).getValue();
        return r / Math.pow(alpha, 2);
    }

    public double[] sample(int n, Random random) {
        return map_sample(D(0), D(1), n, random);
    }

    @Override
    public double[] sample(int n) {
        return sample(n, RandomManager.getThreadRandomAsRandom());
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getMean()
     */
    public double mean() {
        return getMean();
    }
    
    /**
     * Kotlin-style property alias for getRate()
     */
    public double rate() {
        return getRate();
    }
    
    /**
     * Kotlin-style property alias for getSCV()
     */
    public double scv() {
        return getSCV();
    }
    
    /**
     * Kotlin-style property alias for getSkewness()
     */
    public double skewness() {
        return getSkewness();
    }
    
    /**
     * Kotlin-style property alias for getVar()
     */
    public double var() {
        return getVar();
    }
}
