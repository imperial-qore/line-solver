/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.petri;

import jline.util.matrix.Matrix;

/**
 * The min-normal moment closures a Petri-net transition mode needs. Java twin of
 * the MATLAB {@code fluid_min_closure} and {@code fluid_minmulti_closure}.
 *
 * <p>A queueing station needs the TWO-argument {@code min(n,c)}. A transition
 * mode needs the many-argument one: its enabling degree is
 *
 * <pre>    e(m) = min_a ( m_a / w_a )</pre>
 *
 * over every input arc a, and the rate law then caps that at the mode's server
 * count. There is no closed form for the expectation of a min of more than two
 * correlated normals, so this uses the recursion of Clark (1961): the running
 * min is replaced at each step by the normal with its exact first two moments,
 * and the next argument is folded in with the exact bivariate formulas. The
 * deterministic cap is folded in last, by the two-argument closure itself, so a
 * single-arc mode reduces EXACTLY to the closure the queueing methods already
 * use and no second code path exists for it.
 *
 * <p>THE RECURSION IS ORDER DEPENDENT, as Clark's approximation always is: only
 * the first two moments of the running min are kept, so folding the arcs in a
 * different order gives a slightly different answer. The order here is the
 * caller's, i.e. increasing state coordinate, which the layout fixes and which
 * therefore reproduces across the four codebases.
 *
 * <p>Reference: C. E. Clark, "The greatest of a finite set of random
 * variables", Operations Research 9(2):145-162, 1961.
 */
public final class PetriClosures {

    /** The band inside which two arguments are treated as equal at zero variance. */
    public static final double FINE_TOL = 1e-8;

    private PetriClosures() {
    }

    /** E[min(X,Y)] and dE/dE[X] for jointly normal X, Y. */
    public static double[] minClosure(double n, double c, double s2, double vc, double cov) {
        if (Double.isInfinite(c)) {
            return new double[]{n, 1.0};
        }
        double th2 = s2 - 2.0 * cov + vc;
        if (th2 <= 0.0) {
            // A degenerate pair: the min is the smaller of the two exactly. The
            // band matches the other three codebases, so two of them stopping
            // either side of the kink read the same indicator.
            if (c - n > FINE_TOL * Math.max(1.0, Math.abs(n))) {
                return new double[]{n, 1.0};
            }
            return new double[]{c, 0.0};
        }
        double th = Math.sqrt(th2);
        double al = (n - c) / th;
        double phi = normPdf(al);
        double p = 1.0 - normCdf(al);
        double h = n * p + c * (1.0 - p) - th * phi;
        return new double[]{h, p};
    }

    /** The result of the many-argument closure. */
    public static final class MinMulti {
        /** E[min(X_1,...,X_A,c)]. */
        public final double h;
        /** dH/dMU(a): the probability that arc a is the binding one. */
        public final double[] g;
        /** Var[min] before the cap, carried for the caller's report. */
        public final double v;

        MinMulti(double h, double[] g, double v) {
            this.h = h;
            this.g = g;
            this.v = v;
        }
    }

    /**
     * Min-normal closure of E[min(X_1,...,X_A,c)] by Clark's recursion.
     *
     * @param mu means of the arguments, already scaled by the arc weights
     * @param S  covariance of the arguments, symmetric positive semi-definite
     * @param c  deterministic cap (the mode's server count); infinite for none
     */
    public static MinMulti minMultiClosure(double[] mu, Matrix S, double c) {
        int A = mu == null ? 0 : mu.length;
        // A mode with no input arc is enabled at degree one, the convention the
        // exact engines use (Solver_ssa_nrm's spnEnDegree, State.afterGlobalEvent).
        if (A == 0) {
            return new MinMulti(Math.min(1.0, c), new double[0], 0.0);
        }
        double[][] s = new double[A][A];
        if (S != null && S.getNumRows() >= A && S.getNumCols() >= A) {
            for (int i = 0; i < A; i++) {
                for (int j = 0; j < A; j++) {
                    s[i][j] = S.get(i, j);
                }
            }
        }

        double mz = mu[0];
        double vz = s[0][0];
        double[] covz = new double[A];
        for (int i = 0; i < A; i++) {
            covz[i] = s[0][i];
        }
        double[] p = new double[A];
        p[0] = 1.0;

        for (int k = 1; k < A; k++) {
            double th2 = vz + s[k][k] - 2.0 * covz[k];
            if (th2 < 0.0) {
                th2 = 0.0; // a covariance beyond the Cauchy-Schwarz bound is not admissible
            }
            double th = Math.sqrt(th2);
            double pk;
            double mw;
            double vw;
            double[] cw = new double[A];
            if (th > 0.0) {
                double al = (mz - mu[k]) / th;
                double Phi = normCdf(al);
                double phi = normPdf(al);
                pk = 1.0 - Phi;
                mw = mz * pk + mu[k] * Phi - th * phi;
                double e2 = (mz * mz + vz) * pk + (mu[k] * mu[k] + s[k][k]) * Phi
                        - (mz + mu[k]) * th * phi;
                vw = e2 - mw * mw;
                // Clark's moment match can leave a negative variance where the
                // two arguments are nearly identical; the min of two equal
                // normals has the variance of either, which the clamp restores.
                if (vw < 0.0) {
                    vw = 0.0;
                }
                for (int i = 0; i < A; i++) {
                    cw[i] = covz[i] * pk + s[i][k] * Phi;
                }
            } else {
                pk = (mu[k] - mz > FINE_TOL * Math.max(1.0, Math.abs(mz))) ? 1.0 : 0.0;
                mw = Math.min(mz, mu[k]);
                vw = pk * vz + (1.0 - pk) * s[k][k];
                for (int i = 0; i < A; i++) {
                    cw[i] = covz[i] * pk + s[i][k] * (1.0 - pk);
                }
            }
            p[k] = pk;
            mz = mw;
            vz = vw;
            covz = cw;
        }

        double v = vz;
        // The cap, by the two-argument closure itself: c is deterministic, so
        // its variance and its covariance with the running min are both zero.
        double[] capped = minClosure(mz, c, vz, 0.0, 0.0);

        double[] g = new double[A];
        double tail = capped[1];
        for (int a = A - 1; a >= 1; a--) {
            g[a] = (1.0 - p[a]) * tail;
            tail = tail * p[a];
        }
        g[0] = tail;
        return new MinMulti(capped[0], g, v);
    }

    /** The standard normal CDF, without a statistics dependency. */
    public static double normCdf(double z) {
        return 0.5 * org.apache.commons.math3.special.Erf.erfc(-z / Math.sqrt(2.0));
    }

    /** The standard normal density. */
    public static double normPdf(double z) {
        return Math.exp(-0.5 * z * z) / Math.sqrt(2.0 * Math.PI);
    }
}
