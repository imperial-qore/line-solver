/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

/**
 * A Gamma distribution.
 */

public class Gamma extends ContinuousDistribution implements Serializable {

    public Gamma(double shape, double scale) {
        super("Gamma", 2, new Pair<Double, Double>(0.0, Inf));
        this.setParam(1, "alpha", shape);
        this.setParam(2, "beta", scale);
    }

    public static Gamma fitMeanAndSCV(double mean, double scv) {
        double shape = 1.0 / scv;
        double scale = mean / shape;
        return new Gamma(shape, scale);
    }

    @Override
    public double evalCDF(double t) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        return gammaDistribution.cumulativeProbability(t);
    }

    /**
     * Evaluates the probability density function (PDF) at the given point.
     *
     * @param t the point at which to evaluate the PDF
     * @return the PDF value at point t
     */
    public double evalPDF(double t) {
        if (t <= 0) {
            return 0.0;
        }
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        return gammaDistribution.density(t);
    }

    @Override
    public double evalLST(double s) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return FastMath.pow((1.0 / scale) / (s + 1.0 / scale), shape);
    }

    @Override
    public double getMean() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return shape * scale;
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    @Override
    public double getSCV() {
        double shape = (double) this.getParam(1).getValue();
        return 1.0 / shape;
    }

    @Override
    public double getSkewness() {
        double shape = (double) this.getParam(1).getValue();
        return 2.0 / FastMath.sqrt(shape);
    }

    @Override
    public double getVar() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return shape * scale * scale;
    }

    public double[] sample(int n, Random random) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        double[] sampleList = new double[(int) n];
        for (int i = 0; i < n; i++) {
            sampleList[i] = gammaDistribution.inverseCumulativeProbability(random.nextDouble());
        }
        return sampleList;
    }

    @Override
    public double[] sample(int n) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        return gammaDistribution.sample(n);
    }

    public MatrixCell getProcess() {
        // Returns {shape, scale} for gamma distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // shape (alpha)
        representation.set(1, Matrix.singleton((Double) this.getParam(2).getValue()));  // scale (beta)
        return representation;
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getProcess()
     */
    public MatrixCell process() {
        return getProcess();
    }
    
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

