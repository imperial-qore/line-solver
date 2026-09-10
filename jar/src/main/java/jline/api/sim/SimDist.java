/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

/**
 * Normal and Student t quantiles used by the output-analysis routines.
 *
 * <p>Thin wrappers over commons-math3, kept in one place so the MATLAB twins
 * (which reach the same values through {@code erfc}, {@code erfinv} and
 * {@code betaincinv} to avoid a toolbox dependency) and the Python twins (which
 * use SciPy) have a single point of comparison. Agreement across the three is to
 * within 1e-10.
 *
 * @since LINE 3.1.0
 */
public final class SimDist {
    private static final NormalDistribution STD_NORMAL = new NormalDistribution(0.0, 1.0);

    private SimDist() {}

    /**
     * Standard normal cumulative distribution function.
     *
     * @param z the argument
     * @return Phi(z)
     */
    public static double normcdf(double z) {
        return STD_NORMAL.cumulativeProbability(z);
    }

    /**
     * Standard normal quantile function.
     *
     * @param p probability in [0,1]
     * @return Phi inverse of p, infinite at the endpoints
     */
    public static double norminv(double p) {
        if (p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("Probability p=" + p + " must lie in [0,1]");
        }
        if (p == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return STD_NORMAL.inverseCumulativeProbability(p);
    }

    /**
     * Quantile function of Student's t distribution.
     *
     * @param p  probability in [0,1]
     * @param nu degrees of freedom, positive
     * @return the p-quantile of t with nu degrees of freedom
     */
    public static double tinv(double p, double nu) {
        if (!(nu > 0.0)) {
            throw new IllegalArgumentException("Degrees of freedom nu=" + nu + " must be positive");
        }
        if (p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("Probability p=" + p + " must lie in [0,1]");
        }
        if (p == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return new TDistribution(nu).inverseCumulativeProbability(p);
    }
}
