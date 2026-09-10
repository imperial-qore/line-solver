/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

/**
 * Rectangle probability of a multivariate normal, the cell integral behind
 * {@code SolverFluid.getProbAggr} under the moment-closure methods. Java twin of
 * the MATLAB {@code fluid_mvn_rectangle} and of the Python {@code mvn_rectangle}.
 *
 * <p>The integral has no closed form beyond one dimension, so it is evaluated by
 * the separation-of-variables transformation of Genz (1992): the Cholesky factor
 * of C turns the rectangle into an iterated integral over the unit cube whose
 * integrand is a product of normal-CDF differences, and the first coordinate is
 * integrated exactly. The remaining cube is integrated with a DETERMINISTIC
 * Richtmyer lattice rule, frac(k*sqrt(p_j)) over the first primes, averaged with
 * its antithetic reflection. Determinism is required here, not merely
 * convenient: the MATLAB, Java and Python twins must return the same number, and
 * a randomized rule would make them agree only in distribution.</p>
 *
 * <p>C may be SINGULAR, which is the common case: a closed population fixes the
 * sum of the station coordinates, so the covariance of a station holding a whole
 * class is rank deficient. A coordinate whose CONDITIONAL variance vanishes is
 * not integrated; it is a hard constraint, contributing 1 when the conditional
 * mean falls inside its interval and 0 otherwise.</p>
 *
 * @see FluidLyapunov
 */
public final class MvnRectangle {

    /** Lattice points per antithetic pair. */
    public static final int DEFAULT_POINTS = 4096;

    /**
     * First 100 primes, listed rather than sieved so that the MATLAB and Python
     * twins generate the identical lattice.
     */
    private static final int[] PRIMES = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
            179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
            283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541};

    private MvnRectangle() {
    }

    /**
     * Probability that a normal vector falls in the rectangle [a,b].
     *
     * @param m mean vector, length d
     * @param C covariance matrix, d-by-d, symmetric positive semi-definite
     * @param a lower corner, -infinity allowed
     * @param b upper corner, +infinity allowed
     * @return the probability in [0,1]
     */
    public static double probability(double[] m, double[][] C, double[] a, double[] b) {
        return probability(m, C, a, b, DEFAULT_POINTS);
    }

    /**
     * Probability that a normal vector falls in the rectangle [a,b].
     *
     * @param m       mean vector, length d
     * @param C       covariance matrix, d-by-d, symmetric positive semi-definite
     * @param a       lower corner, -infinity allowed
     * @param b       upper corner, +infinity allowed
     * @param npoints lattice points per antithetic pair
     * @return the probability in [0,1]
     */
    public static double probability(double[] m, double[][] C, double[] a, double[] b, int npoints) {
        int d = m.length;
        if (d == 0) {
            return 1.0;
        }
        double[] al = new double[d];
        double[] bu = new double[d];
        for (int i = 0; i < d; i++) {
            al[i] = a[i] - m[i];
            bu[i] = b[i] - m[i];
            if (bu[i] <= al[i]) {
                return 0.0;
            }
        }

        // scale-relative tolerances: dtol decides which coordinate carries noise,
        // ctol whether a deterministic coordinate satisfies its constraint
        double scale = 1.0;
        for (int i = 0; i < d; i++) {
            scale = FastMath.max(scale, FastMath.abs(C[i][i]));
        }
        double dtol = 1e-12 * scale;
        double ctol = 1e-6 * FastMath.sqrt(scale);

        double[][] L = cholPsd(C, dtol);
        int nInt = 0;
        int lastInt = -1;
        for (int i = 0; i < d; i++) {
            if (L[i][i] > 0) {
                nInt++;
                lastInt = i;
            }
        }
        int nw = FastMath.max(0, nInt - 1);
        if (nw > PRIMES.length) {
            throw new RuntimeException(String.format(
                    "The lattice rule carries generators for at most %d integration dimensions, but this "
                            + "rectangle has %d. Aggregate classes before evaluating the cell.",
                    PRIMES.length, nw));
        }
        double[] alpha = new double[nw];
        for (int j = 0; j < nw; j++) {
            alpha[j] = FastMath.sqrt((double) PRIMES[j]);
        }

        int npairs = (nw == 0) ? 1 : npoints;
        int nEval = (nw == 0) ? 1 : 2 * npoints;
        double acc = 0.0;
        double[] w = new double[nw];
        double[] y = new double[d];
        for (int k = 1; k <= npairs; k++) {
            for (int j = 0; j < nw; j++) {
                double v = k * alpha[j];
                w[j] = v - FastMath.floor(v);
            }
            acc += evaluate(L, al, bu, w, y, lastInt, ctol, false);
            if (nw > 0) {
                acc += evaluate(L, al, bu, w, y, lastInt, ctol, true); // antithetic reflection
            }
        }

        double p = acc / nEval;
        return FastMath.min(FastMath.max(p, 0.0), 1.0);
    }

    /**
     * Natural logarithm of {@link #probability(double[], double[][], double[], double[])},
     * negative infinity when the probability is zero.
     *
     * @param m mean vector
     * @param C covariance matrix
     * @param a lower corner
     * @param b upper corner
     * @return log of the rectangle probability
     */
    public static double logProbability(double[] m, double[][] C, double[] a, double[] b) {
        double p = probability(m, C, a, b);
        return p > 0 ? FastMath.log(p) : Double.NEGATIVE_INFINITY;
    }

    /** One point of the transformed integrand, Genz's recursion over the coordinates. */
    private static double evaluate(double[][] L, double[] al, double[] bu, double[] w, double[] y,
                                   int lastInt, double ctol, boolean antithetic) {
        int d = al.length;
        double f = 1.0;
        int kw = 0;
        for (int i = 0; i < d; i++) {
            double s = 0.0;
            for (int j = 0; j < i; j++) {
                s += L[i][j] * y[j];
            }
            if (L[i][i] > 0) {
                double dd = phi((al[i] - s) / L[i][i]);
                double ee = phi((bu[i] - s) / L[i][i]);
                f *= FastMath.max(0.0, ee - dd);
                if (f == 0.0) {
                    return 0.0;
                }
                if (i != lastInt) {
                    double wk = antithetic ? 1.0 - w[kw] : w[kw];
                    kw++;
                    double u = dd + wk * (ee - dd);
                    // the inverse CDF is evaluated strictly inside the unit interval
                    u = FastMath.min(FastMath.max(u, 1e-15), 1.0 - 1e-15);
                    y[i] = phiInv(u);
                }
            } else {
                // zero conditional variance: the coordinate is pinned at s, so the
                // cell is either met or not
                if (s < al[i] - ctol || s > bu[i] + ctol) {
                    return 0.0;
                }
                y[i] = 0.0;
            }
        }
        return f;
    }

    /**
     * Cholesky factor of a symmetric positive SEMI-definite matrix. A vanishing
     * pivot leaves a zero row/column, which the caller reads as a deterministic
     * coordinate rather than as a failure.
     */
    private static double[][] cholPsd(double[][] C, double dtol) {
        int d = C.length;
        double[][] L = new double[d][d];
        for (int i = 0; i < d; i++) {
            double v = C[i][i];
            for (int j = 0; j < i; j++) {
                v -= L[i][j] * L[i][j];
            }
            if (v > dtol) {
                L[i][i] = FastMath.sqrt(v);
                for (int r = i + 1; r < d; r++) {
                    double s = C[r][i];
                    for (int j = 0; j < i; j++) {
                        s -= L[r][j] * L[i][j];
                    }
                    L[r][i] = s / L[i][i];
                }
            } else {
                L[i][i] = 0.0;
                for (int r = i + 1; r < d; r++) {
                    L[r][i] = 0.0;
                }
            }
        }
        return L;
    }

    /** Standard normal CDF. */
    private static double phi(double x) {
        if (Double.isInfinite(x)) {
            return x > 0 ? 1.0 : 0.0;
        }
        return 0.5 * Erf.erfc(-x / FastMath.sqrt(2.0));
    }

    /** Standard normal quantile. */
    private static double phiInv(double u) {
        return FastMath.sqrt(2.0) * Erf.erfInv(2.0 * u - 1.0);
    }
}
