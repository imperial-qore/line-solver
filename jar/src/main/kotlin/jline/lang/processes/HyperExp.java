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

import static jline.api.mam.Map_hyperexpKt.map_hyperexp;
import static jline.api.mam.Map_sampleKt.map_sample;

/**
 * A hyper-exponential distribution.
 */

@SuppressWarnings("unchecked")
public class HyperExp extends Markovian implements Serializable {

    private final long nPhases;

    public HyperExp(double p, double lambda1, double lambda2) {
        super("HyperExp", 1);

        this.setParam(1, "p", p);
        this.setParam(2, "lambda1", lambda1);
        this.setParam(3, "lambda2", lambda2);

        nPhases = 2;
        Matrix D0 = new Matrix(2, 2, 4);
        Matrix D1 = new Matrix(2, 2, 4);
        D0.set(0, 0, -lambda1);
        D0.set(1, 1, -lambda2);
        D1.set(0, 0, lambda1 * p);
        D1.set(0, 1, lambda1 * (1.0 - p));
        D1.set(1, 0, lambda2 * p);
        D1.set(1, 1, lambda2 * (1.0 - p));
        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        setProcess(rep);
    }

    public HyperExp(double p, double lambda) {
        this(p, lambda, lambda);
    }

    /**
     * Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
     */
    public static HyperExp fitMeanAndSCV(double mean, double scv) {
        double p, mu1, mu2;
        MatrixCell D = map_hyperexp(mean, scv, 0);
        mu1 = -D.get(0).value();
        mu2 = -D.get(0).get(1, 1);
        p = D.get(1).value() / mu1;
        HyperExp he = new HyperExp(p, mu1, mu2);
        he.immediate = mean < GlobalConstants.CoarseTol;
        return he;
    }

    /**
     * Fit distribution with given squared coefficient of variation and balanced means i.e.,
     * p/mu1 = (1-p)/mu2
     */
    public static HyperExp fitMeanAndSCVBalanced(double mean, double scv) {
        double p, mu1, mu2;
        mu1 = -(2.0 * (Math.sqrt((scv - 1) / (scv + 1)) / 2.0 - 0.5)) / mean;
        p = 0.5 - FastMath.sqrt((scv - 1) / (scv + 1)) / 2.0;
        if (mu1 < 0 || p < 0 || p > 1) {
            p = FastMath.sqrt((scv - 1) / (scv + 1)) / 2.0 + 0.5;
            mu1 = (2 * (Math.sqrt((scv - 1) / (scv + 1)) / 2.0 + 0.5)) / mean;
        }
        mu2 = (1 - p) / p * mu1;
        HyperExp he = new HyperExp(p, mu1, mu2);
        he.immediate = mean < GlobalConstants.CoarseTol;
        return he;
    }

    public double evalCDF(double t) {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return p * (1 - FastMath.exp(-mu1 * t)) + (1 - p) * (1 - FastMath.exp(-mu2 * t));
        } else {
            return super.evalCDF(t);
        }
    }

    public double evalLST(double s) {
        return super.evalLST(s);
    }

    public double getMean() {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return p / mu1 + (1 - p) / mu2;
        } else {
            return super.getMean();
        }
    }

    public long getNumberOfPhases() {
        return nPhases;
    }

    public double getRate() {
        return 1.0 / getMean();
    }

    public double getSCV() {
        if (this.nPhases == 2) {
            double p = (double) this.getParam(1).getValue();
            double mu1 = (double) this.getParam(2).getValue();
            double mu2 = (double) this.getParam(3).getValue();
            return (2 * (p / FastMath.pow(mu1, 2) + (1 - p) / FastMath.pow(mu2, 2)) - FastMath.pow(p / mu1 + (1 - p) / mu2, 2)) / FastMath.pow(p / mu1 + (1 - p) / mu2, 2);
        } else {
            return super.getSCV();
        }
    }

    public double getSkewness() {
        return super.getSkewness();
    }

    public double getVar() {
        return this.getSCV() * FastMath.pow(this.getMean(), 2);
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    @Override
    public double[] sample(int n, Random random) {
        return map_sample(D(0), D(1), n, random);
    }

    public String toString() {
        return String.format("jline.HyperExp(%f)", this.getRate());
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
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numberOfPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numPhases() {
        return getNumberOfPhases();
    }
}
