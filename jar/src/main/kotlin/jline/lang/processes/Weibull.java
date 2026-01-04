/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Maths;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Weibull distribution
 */
public class Weibull extends ContinuousDistribution implements Serializable {

    public Weibull(double shape, double scale) {
        super("Weibull", 2, new Pair<Double, Double>(0.0, Inf));
        if (shape < 0) {
            line_error(mfilename(new Object() {
            }), "shape parameter must be >= 0.0");
        }
        this.setParam(1, "alpha", scale);
        this.setParam(2, "r", shape);
    }

    public static Weibull fitMeanAndSCV(double mean, double scv) {
        // Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
        double c = FastMath.sqrt(scv);
        double r = FastMath.pow(c, -1.086); //Justus approximation (1976)
        double alpha = mean / Maths.gammaFunction(1 + 1.0 / r);
        return new Weibull(r, alpha);
    }

    @Override
    public double evalCDF(double t) {
        // Evaluate the cumulative distribution function at t
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        if (t <= 0.0) {
            return 0.0;
        } else {
            return 1 - FastMath.exp(-Math.pow(t / alpha, r));
        }
    }

    /**
     * Evaluates the probability density function (PDF) at the given point.
     * Weibull PDF: (r/α)(x/α)^(r-1) * exp(-(x/α)^r)
     *
     * @param t the point at which to evaluate the PDF
     * @return the PDF value at point t
     */
    public double evalPDF(double t) {
        if (t <= 0) {
            return 0.0;
        }
        double alpha = (double) this.getParam(1).getValue(); // scale parameter
        double r = (double) this.getParam(2).getValue(); // shape parameter
        return (r / alpha) * FastMath.pow(t / alpha, r - 1) * FastMath.exp(-FastMath.pow(t / alpha, r));
    }

    @Override
    public double evalLST(double s) {
        // The LST of Weibull distribution doesn't have a simple closed form.
        // We compute it numerically using the definition: E[e^(-sX)] = ∫₀^∞ e^(-sx) f(x) dx
        // where f(x) is the Weibull PDF
        
        double alpha = (double) this.getParam(1).getValue(); // scale parameter
        double r = (double) this.getParam(2).getValue(); // shape parameter
        
        // Use numerical integration
        // For practical purposes, we integrate up to a reasonable upper bound
        double upperBound = alpha * FastMath.pow(-FastMath.log(1e-10), 1.0/r); // covers most of the distribution
        int n = 1000;
        double dx = upperBound / n;
        double sum = 0.0;
        
        for (int i = 1; i <= n; i++) {
            double x = i * dx;
            if (x > 0) {
                // Weibull PDF: (r/α)(x/α)^(r-1) * exp(-(x/α)^r)
                double pdf = (r / alpha) * FastMath.pow(x / alpha, r - 1) * 
                           FastMath.exp(-FastMath.pow(x / alpha, r));
                sum += FastMath.exp(-s * x) * pdf;
            }
        }
        
        return sum * dx;
    }

    @Override
    public double getMean() {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        // Mean of Weibull: scale * Gamma(1 + 1/shape)
        return alpha * Maths.gammaFunction(1 + 1.0 / r);
    }

    @Override
    public double getRate() {
        return 1.0 / this.getMean();
    }

    @Override
    public double getSCV() {
        // Get distribution squared coefficient of variation (SCV = variance / mean^2)
        return getVar() / FastMath.pow(getMean(), 2);
    }

    @Override
    public double getSkewness() {
        double r = (double) this.getParam(2).getValue();
        double gamma1 = Maths.gammaFunction(1 + 1.0 / r);
        double gamma2 = Maths.gammaFunction(1 + 2.0 / r);
        double gamma3 = Maths.gammaFunction(1 + 3.0 / r);

        double numerator = gamma3 - 3 * gamma2 * gamma1 + 2 * FastMath.pow(gamma1, 3);
        double denominator = FastMath.pow(gamma2 - FastMath.pow(gamma1, 2), 1.5);

        return numerator / denominator;
    }

    @Override
    public double getVar() {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();
        return FastMath.pow(alpha, 2) * (Maths.gammaFunction(1 + 2 / r) - FastMath.pow(Maths.gammaFunction(1 + 1.0 / r), 2));
    }

    @Override
    public double[] sample(int n, Random random) {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        double[] samples = new double[(int) n];
        // Apache Commons Math WeibullDistribution constructor takes (shape, scale) order
        WeibullDistribution weibullDistribution = new WeibullDistribution(r, alpha);
        for (int i = 0; i < n; i++) {
            samples[i] = weibullDistribution.inverseCumulativeProbability(random.nextDouble());
        }
        return samples;
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

    public MatrixCell getProcess() {
        // Returns {shape (r), scale (alpha)} for Weibull distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(2).getValue()));  // shape (r)
        representation.set(1, Matrix.singleton((Double) this.getParam(1).getValue()));  // scale (alpha)
        return representation;
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
     * Kotlin-style property alias for getProcess()
     */
    public MatrixCell process() {
        return getProcess();
    }
}

