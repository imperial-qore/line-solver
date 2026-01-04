/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Random;

import static jline.api.mam.Map_normalizeKt.map_normalize;
import static jline.api.mam.Map_sampleKt.map_sample;

/**
 * A general phase-type (PH) distribution
 */
@SuppressWarnings("unchecked")
public class PH extends Markovian {

    public PH(Matrix alpha, Matrix T) {
        super("PH", 3);

        // Ensure alpha is a row vector (1 x n)
        // If passed as column vector (n x 1), transpose it
        Matrix alphaRow = alpha;
        if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
            alphaRow = alpha.transpose();
        }

        nPhases = alphaRow.getNumElements();

        this.setParam(1, "n", nPhases);
        this.setParam(2, "alpha", alphaRow);
        this.setParam(3, "T", T);
        MatrixCell rep = new MatrixCell();
        rep.set(0, T);
        rep.set(1, T.mult(alphaRow.repmat(nPhases, 1)).scale(-1));
        this.setProcess(rep);
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at time t.
     * For a PH distribution, CDF(t) = 1 - alpha * exp(T*t) * e
     * where alpha is the initial probability vector, T is the subgenerator,
     * and e is a vector of ones.
     *
     * @param t the time value
     * @return the CDF value at time t
     */
    @Override
    public double evalCDF(double t) {
        if (t <= 0) {
            return 0.0;
        }

        try {
            Matrix alpha = getInitProb();
            Matrix T = getSubgenerator();

            // Compute exp(T*t)
            Matrix Tt = Matrix.scaleMult(T, t);
            Matrix expTt = Maths.matrixExp(Tt);

            // Create vector of ones
            Matrix ones = Matrix.ones(T.getNumRows(), 1);

            // CDF(t) = 1 - alpha * exp(T*t) * e
            Matrix result = alpha.mult(expTt).mult(ones);
            return 1.0 - result.get(0, 0);
        } catch (Exception e) {
            // Fallback: approximate using exponential distribution with same mean
            double mean = getMean();
            if (mean > 0) {
                return 1.0 - Math.exp(-t / mean);
            } else {
                return t > 0 ? 1.0 : 0.0;
            }
        }
    }

    @Override
    public double evalLST(double s) {
        return super.evalLST(s);
    }

    public Matrix getInitProb() {
        return (Matrix) this.getParam(2).getValue();
    }

    @Override
    public double getMean() {
        return super.getMean();
    }

    @Override
    public long getNumberOfPhases() {
        return (long) this.getParam(1).getValue();
    }

    @Override
    public MatrixCell getProcess() {
        MatrixCell res;
        Matrix T = getSubgenerator();

        Matrix ones = new Matrix(T.getNumCols(), 1, T.getNumCols());
        Matrix Te = new Matrix(0, 0, 0);
        Matrix Tepie = new Matrix(0, 0, 0);
        ones.fill(1.0);
        T.mult(ones, Te);
        Te.mult(this.getInitProb(), Tepie);
        //Tepie.removeZeros(0);
        Tepie.changeSign();

        res = map_normalize(T, Tepie);
        return res;
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    @Override
    public double getSCV() {
        return super.getSCV();
    }

    @Override
    public double getSkewness() {
        return super.getSkewness();
    }

    public Matrix getSubgenerator() {
        return (Matrix) this.getParam(3).getValue();
    }

    @Override
    public double getVar() {
        return this.getSCV() * Math.pow(this.getMean(), 2);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    @Override
    public double[] sample(int n, Random random) {
        return map_sample(D(0), D(1), n, random);
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
     * Kotlin-style property alias for getInitProb()
     */
    public Matrix initProb() {
        return getInitProb();
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
    
    /**
     * Kotlin-style property alias for getProcess()
     */
    public MatrixCell process() {
        return getProcess();
    }
    
    /**
     * Kotlin-style property alias for getSubgenerator()
     */
    public Matrix subgenerator() {
        return getSubgenerator();
    }
}
