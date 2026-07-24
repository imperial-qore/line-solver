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

import static jline.api.mam.Map_hyperexp.map_hyperexp;
import static jline.api.mam.Map_sample.map_sample;

/**
 * A hyper-exponential distribution.
 */

@SuppressWarnings("unchecked")
public class HyperExp extends Markovian implements Serializable {

    // see _kb/01-model-classes.md (Java process-construction notes) for rationale

    // Branch probabilities and phase rates, in phase order. Both have length
    // nPhases and are set by every constructor, so that a HyperExp of any phase
    // count can be serialized without going through the 2-phase
    // (p, lambda1, lambda2) named parameters.
    private final double[] pVector;
    private final double[] lambdaVector;

    public HyperExp(double p, double lambda1, double lambda2) {
        super("HyperExp", 1);

        this.setParam(1, "p", p);
        this.setParam(2, "lambda1", lambda1);
        this.setParam(3, "lambda2", lambda2);

        this.pVector = new double[]{p, 1.0 - p};
        this.lambdaVector = new double[]{lambda1, lambda2};

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
     * Creates an n-phase hyper-exponential distribution: with probability p[i]
     * the sample is exponential with rate lambda[i].
     * <p>
     * This mirrors the MATLAB HyperExp representation exactly:
     * D0 = -diag(lambda), D1 = -D0*p(:)*ones(1,n), i.e. D1(i,j) = lambda(i)*p(j).
     * For n = 2 this coincides with the (p, lambda1, lambda2) constructor.
     *
     * @param p      branch probabilities, one per phase (must sum to 1)
     * @param lambda phase rates, one per phase
     */
    public HyperExp(double[] p, double[] lambda) {
        super("HyperExp", 1);

        if (p == null || lambda == null) {
            throw new IllegalArgumentException("HyperExp: p and lambda must not be null");
        }
        if (p.length != lambda.length) {
            throw new IllegalArgumentException(
                    "HyperExp: p has " + p.length + " entries but lambda has " + lambda.length);
        }
        int n = p.length;
        if (n < 1) {
            throw new IllegalArgumentException("HyperExp: at least one phase is required");
        }
        double psum = 0.0;
        for (int i = 0; i < n; i++) {
            if (lambda[i] <= 0.0) {
                throw new IllegalArgumentException(
                        "HyperExp: rate of phase " + (i + 1) + " must be positive, got " + lambda[i]);
            }
            if (p[i] < 0.0) {
                throw new IllegalArgumentException(
                        "HyperExp: probability of phase " + (i + 1) + " must be non-negative, got " + p[i]);
            }
            psum += p[i];
        }
        if (FastMath.abs(psum - 1.0) > GlobalConstants.CoarseTol) {
            throw new IllegalArgumentException("HyperExp: branch probabilities sum to " + psum + ", not 1");
        }

        this.pVector = new double[n];
        this.lambdaVector = new double[n];
        System.arraycopy(p, 0, this.pVector, 0, n);
        System.arraycopy(lambda, 0, this.lambdaVector, 0, n);

        // see _kb/01-model-classes.md (Java process-construction notes) for rationale
        if (n == 2) {
            this.setParam(1, "p", p[0]);
            this.setParam(2, "lambda1", lambda[0]);
            this.setParam(3, "lambda2", lambda[1]);
        } else {
            this.setParam(1, "p", this.pVector);
            this.setParam(2, "lambda", this.lambdaVector);
        }

        nPhases = n;
        Matrix D0 = new Matrix(n, n, n);
        Matrix D1 = new Matrix(n, n, n * n);
        for (int i = 0; i < n; i++) {
            D0.set(i, i, -lambda[i]);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D1.set(i, j, lambda[i] * p[j]);
            }
        }
        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        setProcess(rep);
    }

    /**
     * Creates an n-phase hyper-exponential distribution from Matrix vectors.
     *
     * @param p      branch probabilities, one per phase (must sum to 1)
     * @param lambda phase rates, one per phase
     */
    public HyperExp(Matrix p, Matrix lambda) {
        this(toArray(p, "p"), toArray(lambda, "lambda"));
    }

    private static double[] toArray(Matrix m, String what) {
        if (m == null) {
            throw new IllegalArgumentException("HyperExp: " + what + " must not be null");
        }
        if (m.getNumRows() != 1 && m.getNumCols() != 1) {
            throw new IllegalArgumentException(
                    "HyperExp: " + what + " must be a vector, got " + m.getNumRows() + "x" + m.getNumCols());
        }
        return m.toArray1D();
    }

    /**
     * Gets the branch probabilities, one per phase.
     *
     * @return a copy of the branch probability vector
     */
    public double[] getP() {
        return pVector.clone();
    }

    /**
     * Gets the phase rates, one per phase.
     *
     * @return a copy of the rate vector
     */
    public double[] getLambda() {
        return lambdaVector.clone();
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
            // F(t) = sum_i p_i * (1 - exp(-lambda_i * t))
            double cdf = 0.0;
            for (int i = 0; i < pVector.length; i++) {
                cdf += pVector[i] * (1 - FastMath.exp(-lambdaVector[i] * t));
            }
            return cdf;
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
            // E[X] = sum_i p_i / lambda_i
            return firstMoment();
        }
    }

    /**
     * Exact first moment of the n-phase hyper-exponential: sum_i p_i / lambda_i.
     */
    private double firstMoment() {
        double m1 = 0.0;
        for (int i = 0; i < pVector.length; i++) {
            m1 += pVector[i] / lambdaVector[i];
        }
        return m1;
    }

    /**
     * Exact second moment of the n-phase hyper-exponential: sum_i 2 * p_i / lambda_i^2.
     */
    private double secondMoment() {
        double s = 0.0;
        for (int i = 0; i < pVector.length; i++) {
            s += pVector[i] / FastMath.pow(lambdaVector[i], 2);
        }
        return 2 * s;
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
            // SCV = E[X^2]/E[X]^2 - 1
            double m1 = firstMoment();
            return (secondMoment() - FastMath.pow(m1, 2)) / FastMath.pow(m1, 2);
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
