/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Pareto distribution
 */
public class Pareto extends ContinuousDistribution implements Serializable {

    public Pareto(double shape, double scale) {
        super("Pareto", 2, new Pair<Double, Double>(0.0, Inf));
        if (shape < 2) {
            line_error(mfilename(new Object() {
            }), "shape parameter must be >= 2.0");
        }
        this.setParam(1, "alpha", shape);
        this.setParam(2, "k", scale);
    }

    public static Pareto fitMeanAndSCV(double mean, double scv) {
        // Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
        // For Pareto distribution with shape alpha and scale k:
        //   Mean = alpha*k / (alpha-1)
        //   SCV = 1 / (alpha*(alpha-2))
        //
        // Solving for alpha from SCV:
        //   alpha*(alpha-2) = 1/SCV
        //   alpha^2 - 2*alpha - 1/SCV = 0
        //   alpha = 1 + sqrt(1 + 1/SCV)  (taking positive root, need alpha > 2)
        //
        // Then scale k from mean:
        //   k = mean * (alpha-1) / alpha
        if (mean <= 0 || scv <= 0) {
            line_error(mfilename(new Object() {
            }), "Mean and SCV of the Pareto distribution must be positive.");
        }
        double shape = 1 + FastMath.sqrt(1 + 1.0 / scv);
        double scale = mean * (shape - 1) / shape;
        return new Pareto(shape, scale);
    }

    public static double gpcdf(double x, double k, double sigma, double theta) {
        // Check for invalid sigma
        if (sigma <= 0) {
            return Double.NaN;
        }

        // Calculate the (x - theta) / sigma term
        double z = (x - theta) / sigma;

        // Return 0 for out-of-range values
        if (z < 0) {
            return 0.0;
        }

        // Handle the k == 0 case
        if (Math.abs(k) < FastMath.ulp(1.0)) {
            return -Math.expm1(-z);
        }

        // Compute the CDF value for k != 0
        double t = z * k;

        // When k < 0, the support is 0 <= x/sigma <= -1/k.
        if (t <= -1 && k < -Math.ulp(1.0)) {
            return 1.0;
        }

        return -Math.expm1((-1.0 / k) * FastMath.log1p(t));
    }

    @Override
    public double evalCDF(double t) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        double k = 1 / shape;
        double sigma = scale * k;
        return gpcdf(t, k, sigma, sigma / k);
    }

    /**
     * Evaluates the probability density function (PDF) at the given point.
     * Pareto PDF: α*k^α / x^(α+1) for x >= k
     *
     * @param t the point at which to evaluate the PDF
     * @return the PDF value at point t
     */
    public double evalPDF(double t) {
        double alpha = (double) this.getParam(1).getValue(); // shape parameter
        double k = (double) this.getParam(2).getValue(); // scale parameter (minimum value)
        if (t < k) {
            return 0.0;
        }
        return alpha * FastMath.pow(k, alpha) / FastMath.pow(t, alpha + 1);
    }

    @Override
    public double evalLST(double s) {
        // The LST of Pareto distribution doesn't have a simple closed form.
        // We compute it numerically using the definition: E[e^(-sX)] = ∫_{k}^∞ e^(-sx) f(x) dx
        // where f(x) is the Pareto PDF: α*k^α / x^(α+1) for x >= k
        
        double alpha = (double) this.getParam(1).getValue(); // shape parameter
        double k = (double) this.getParam(2).getValue(); // scale parameter
        
        // Use numerical integration
        // For practical purposes, we integrate up to a reasonable upper bound
        double upperBound = k * FastMath.pow(1000, 1.0/alpha); // covers most of the distribution
        int n = 1000;
        double dx = (upperBound - k) / n;
        double sum = 0.0;
        
        for (int i = 1; i <= n; i++) {
            double x = k + i * dx;
            if (x >= k) {
                // Pareto PDF: α*k^α / x^(α+1)
                double pdf = alpha * FastMath.pow(k, alpha) / FastMath.pow(x, alpha + 1);
                sum += FastMath.exp(-s * x) * pdf;
            }
        }
        
        return sum * dx;
        
        /* OLD IMPLEMENTATION (using RombergIntegrator):
        double alpha = (double) this.getParam(1).getValue();
        double k = (double) this.getParam(2).getValue();
        RombergIntegrator integrator = new RombergIntegrator();
        UnivariateFunction function = new UnivariateFunction() {
            public double value(double t) {
                return FastMath.pow(t, -alpha - 1) * FastMath.exp(-t);
            }
        };

        double IL = integrator.integrate(10000, function, s * k, Inf);
        return alpha * FastMath.pow(k, alpha) * IL * FastMath.pow(s, 1 + alpha) / s;
        */
    }

    @Override
    public double getMean() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();

        if (shape <= 1) {
            return Inf;
        } else {
            return shape * scale / (shape - 1);
        }
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    @Override
    public double getSCV() {
        double shape = (double) this.getParam(1).getValue();
        if (shape <= 2) {
            return Inf;
        } else {
            double scale = (double) this.getParam(2).getValue();
            double var = FastMath.pow(scale, 2) * shape / FastMath.pow(shape - 1, 2) / (shape - 2);
            double ex = shape * scale / (shape - 1);
            return var / FastMath.pow(ex, 2);
        }
    }

    @Override
    public double getSkewness() {
        double shape = (double) this.getParam(1).getValue();
        if (shape <= 3) {
            return Inf;
        } else {
            return 2 * (1 + shape) / (shape - 3) * FastMath.sqrt((shape - 2) / shape);
        }
    }

    @Override
    public double getVar() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();

        if (shape <= 2) {
            return Inf;
        } else {
            return (scale * scale * shape) / ((shape - 1) * (shape - 1) * (shape - 2));
        }
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[(int) n];
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();

        for (int i = 0; i < n; i++) {
            double u = random.nextDouble();
            double sample = scale / FastMath.pow(u, 1.0 / shape);
            samples[i] = sample;
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
        // Returns {shape (alpha), scale (k)} for Pareto distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // shape (alpha)
        representation.set(1, Matrix.singleton((Double) this.getParam(2).getValue()));  // scale (k)
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

