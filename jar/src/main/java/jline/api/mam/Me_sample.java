package jline.api.mam;

import java.util.List;
import java.util.Random;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Matrix Exponential (ME) sampling by numerical inversion of the exact CDF.
 *
 * <p>An ME distribution with representation (alpha, A) has
 * {@code F(t) = 1 - alpha*expm(A*t)*e} and density
 * {@code f(t) = -alpha*expm(A*t)*A*e}. Unlike a phase-type distribution, alpha
 * may contain negative entries and A may have negative off-diagonal entries, so
 * the CTMC walk used by {@link Map_sample} is not applicable: the walk assumes a
 * probabilistic interpretation of the (alpha, A) pair and silently produces
 * garbage when that interpretation does not hold. Inversion of F is valid for
 * every ME representation.</p>
 *
 * <p>The inversion table is built once per {@link MeSampler} instance: a single
 * matrix exponential {@code expm(A*h)} generates the whole grid by repeated
 * vector-matrix products, the grid is extended until the survival function is
 * below {@value #TAIL_MASS_TOL}, and quantiles above the last tabulated CDF
 * value are obtained by exponential extrapolation with the dominant eigenvalue
 * of A rather than being clamped to the grid endpoint. Each quantile is refined
 * by Newton steps on the exact CDF/density, so the result is near-exact rather
 * than piecewise linear.</p>
 *
 * @since LINE 3.0
 */
public final class Me_sample {
    private Me_sample() {}

    /** Number of tabulated grid points (including t = 0). */
    static final int GRID_POINTS = 1001;
    /** Maximum number of doublings applied when extending the grid horizon. */
    static final int MAX_DOUBLINGS = 40;
    /** Survival mass below which the grid horizon is considered sufficient. */
    static final double TAIL_MASS_TOL = 1e-12;
    /** Number of Newton refinement steps applied to each inverted quantile. */
    static final int NEWTON_STEPS = 3;
    /** Maximum value of {@code ||A||*d} handled by a single Taylor substep. */
    static final double TAYLOR_THETA = 0.5;
    /** Maximum number of Taylor substeps before falling back to expm. */
    static final int TAYLOR_MAX_STEPS = 64;
    /** Maximum number of Taylor terms per substep. */
    static final int TAYLOR_MAX_TERMS = 30;

    /**
     * Stateful ME sampler that builds the inverse-CDF table once and draws
     * independent variates from it.
     *
     * <p>The inter-event times of an ME renewal process are i.i.d., so unlike
     * {@link Map_sample.MapSampler} this sampler carries no phase across calls;
     * the state it holds is the inversion table, which is what makes repeated
     * single-sample calls (the LDES usage pattern) cheap.</p>
     */
    public static final class MeSampler {
        private final int n;
        private final Matrix aMatrix;
        private final double[][] a;
        private final double[] rowSumA;
        private final double normA;
        private final double[] tGrid;
        private final double[] cdfGrid;
        private final double[][] wGrid;
        private final double h;
        private final int gridSize;
        private final double tEnd;
        private final double sEnd;
        private final double eta;
        private final boolean exponential;
        private final double expRate;

        /**
         * Builds the inversion table for the ME representation (alpha, A).
         *
         * @param alpha the initial vector, as a row or column vector
         * @param A     the ME matrix parameter (square, eigenvalues with negative real part)
         */
        public MeSampler(Matrix alpha, Matrix A) {
            Matrix alphaRow = alpha;
            if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
                alphaRow = alpha.transpose();
            }
            this.n = A.getNumRows();
            this.aMatrix = A;
            this.a = toArray2D(A);
            this.rowSumA = new double[n];
            double maxRowAbs = 0.0;
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                double abs = 0.0;
                for (int j = 0; j < n; j++) {
                    s += a[i][j];
                    abs += Math.abs(a[i][j]);
                }
                rowSumA[i] = s;
                if (abs > maxRowAbs) maxRowAbs = abs;
            }
            this.normA = maxRowAbs;

            double[] alphaVec = new double[n];
            for (int j = 0; j < n; j++) {
                alphaVec[j] = alphaRow.get(0, j);
            }

            boolean isExp = (n == 1) && Math.abs(alphaVec[0] - 1.0) < 1e-12 && a[0][0] < 0.0;
            this.exponential = isExp;
            this.expRate = isExp ? -a[0][0] : 0.0;

            if (isExp) {
                this.gridSize = 0;
                this.h = 0.0;
                this.tGrid = null;
                this.cdfGrid = null;
                this.wGrid = null;
                this.tEnd = 0.0;
                this.sEnd = 0.0;
                this.eta = -expRate;
                return;
            }

            double mean = Me_mean.me_mean(alphaRow, A);
            double var = Me_var.me_var(alphaRow, A);
            double sigma = (var > 0.0 && !Double.isNaN(var)) ? Math.sqrt(var) : 0.0;
            double horizon = mean + 10.0 * sigma;
            if (!(horizon > 0.0) || Double.isInfinite(horizon) || Double.isNaN(horizon)) {
                horizon = 1.0;
            }
            for (int k = 0; k < MAX_DOUBLINGS; k++) {
                if (survivalByExpm(alphaVec, horizon) < TAIL_MASS_TOL) {
                    break;
                }
                horizon *= 2.0;
            }

            this.gridSize = GRID_POINTS;
            this.h = horizon / (gridSize - 1);
            this.tGrid = new double[gridSize];
            this.cdfGrid = new double[gridSize];
            this.wGrid = new double[gridSize][];

            double[][] eh = toArray2D(A.scale(h).expm_higham());
            wGrid[0] = alphaVec;
            tGrid[0] = 0.0;
            cdfGrid[0] = clampUnit(1.0 - sum(alphaVec));
            for (int i = 1; i < gridSize; i++) {
                wGrid[i] = vecTimesMat(wGrid[i - 1], eh);
                tGrid[i] = i * h;
                double c = clampUnit(1.0 - sum(wGrid[i]));
                if (c < cdfGrid[i - 1]) {
                    c = cdfGrid[i - 1];
                }
                cdfGrid[i] = c;
            }
            this.tEnd = tGrid[gridSize - 1];
            double survival = 1.0 - cdfGrid[gridSize - 1];
            this.sEnd = (survival > 0.0) ? survival : 0.0;
            this.eta = dominantRate(A, mean);
        }

        /** Survival probability at t computed directly from a matrix exponential. */
        private double survivalByExpm(double[] alphaVec, double t) {
            double[][] e = toArray2D(aMatrix.scale(t).expm_higham());
            double[] w = vecTimesMat(alphaVec, e);
            double s = sum(w);
            return (s > 0.0) ? s : 0.0;
        }

        /**
         * Draws one variate by inverting the ME cumulative distribution function.
         *
         * @param random the uniform random source
         * @return a sample from the ME distribution
         */
        public double next(Random random) {
            double u = random.nextDouble();
            if (exponential) {
                return -Math.log(1.0 - u) / expRate;
            }
            if (u <= cdfGrid[0]) {
                return 0.0;
            }
            if (u >= cdfGrid[gridSize - 1]) {
                double tailProb = 1.0 - u;
                if (sEnd <= 0.0 || tailProb <= 0.0 || !(eta < 0.0)) {
                    return tEnd;
                }
                double x = tEnd + Math.log(sEnd / tailProb) / (-eta);
                return (x > tEnd) ? x : tEnd;
            }

            int lo = 0;
            int hi = gridSize - 1;
            while (hi - lo > 1) {
                int mid = (lo + hi) >>> 1;
                if (cdfGrid[mid] <= u) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            double denom = cdfGrid[lo + 1] - cdfGrid[lo];
            double x = (denom > 0.0)
                    ? tGrid[lo] + (u - cdfGrid[lo]) / denom * h
                    : tGrid[lo];
            double left = tGrid[lo];
            double right = tGrid[lo] + h;

            for (int k = 0; k < NEWTON_STEPS; k++) {
                double[] w = expmPropagate(wGrid[lo], a, aMatrix, normA, x - left);
                double surv = sum(w);
                double f = -dot(w, rowSumA);
                if (!(f > 0.0)) {
                    break;
                }
                double err = (1.0 - surv) - u;
                if (Math.abs(err) < 1e-14) {
                    break;
                }
                double xn = x - err / f;
                if (!(xn > left) || !(xn < right)) {
                    break;
                }
                boolean converged = Math.abs(xn - x) <= 1e-15 * Math.max(1.0, Math.abs(x));
                x = xn;
                if (converged) {
                    break;
                }
            }
            return x;
        }
    }

    /**
     * Generates random samples from a Matrix Exponential (ME) distribution by
     * numerical inversion of its exact cumulative distribution function.
     *
     * @param alpha  the initial vector of the ME representation
     * @param A      the matrix parameter of the ME representation
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of n i.i.d. samples
     */
    public static double[] me_sample(Matrix alpha, Matrix A, long n, Random random) {
        MeSampler sampler = new MeSampler(alpha, A);
        double[] samples = new double[(int) n];
        for (int i = 0; i < (int) n; i++) {
            samples[i] = sampler.next(random);
        }
        return samples;
    }

    /**
     * Generates random samples from an ME distribution stored in a MatrixCell
     * {D0 = A, D1 = -A*e*alpha}.
     *
     * <p>D1 of an ME process is the rank-one matrix {@code (-A*e)*alpha}, so
     * alpha is recovered as the row of D1 with the largest {@code |(-A*e)|}
     * entry, divided by that entry. If the recovered alpha does not reproduce
     * D1 the pair is a general Rational Arrival Process rather than an ME
     * renewal process, and sampling is delegated to {@link Rap_sample}.</p>
     *
     * @param ME     the process as a MatrixCell {D0, D1}
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of n samples
     */
    public static double[] me_sample(MatrixCell ME, long n, Random random) {
        Matrix alpha = me_alpha(ME);
        if (alpha != null) {
            return me_sample(alpha, ME.get(0), n, random);
        }
        return Rap_sample.rap_sample(ME, n, random);
    }

    /**
     * Recovers the initial vector alpha of an ME process stored as a MatrixCell
     * {D0 = A, D1 = -A*e*alpha}.
     *
     * <p>D1 is rank one for an ME renewal process, so alpha is read off the row
     * with the largest {@code |(-A*e)|} entry and validated by reconstructing
     * the whole of D1. A null return means the pair is not an ME renewal
     * process (the caller must then treat it as a general RAP).</p>
     *
     * @param ME the process as a MatrixCell {D0, D1}
     * @return the recovered alpha as a 1 x n row vector, or null if D1 is not
     *         the rank-one matrix implied by an ME representation
     */
    public static Matrix me_alpha(MatrixCell ME) {
        Matrix A = ME.get(0);
        Matrix D1 = ME.get(1);
        int nPhases = A.getNumRows();
        if (D1 == null || D1.getNumRows() != nPhases || D1.getNumCols() != nPhases) {
            return null;
        }

        double[][] am = toArray2D(A);
        double[][] d1 = toArray2D(D1);
        double[] exitRate = new double[nPhases];
        int best = -1;
        double bestAbs = 0.0;
        double d1Scale = 1.0;
        for (int i = 0; i < nPhases; i++) {
            double s = 0.0;
            for (int j = 0; j < nPhases; j++) {
                s += am[i][j];
                d1Scale = Math.max(d1Scale, Math.abs(d1[i][j]));
            }
            exitRate[i] = -s;
            if (Math.abs(exitRate[i]) > bestAbs) {
                bestAbs = Math.abs(exitRate[i]);
                best = i;
            }
        }
        if (best < 0 || bestAbs <= 1e-14) {
            return null;
        }

        double[] alphaVec = new double[nPhases];
        Matrix alpha = new Matrix(1, nPhases);
        for (int j = 0; j < nPhases; j++) {
            alphaVec[j] = d1[best][j] / exitRate[best];
            alpha.set(0, j, alphaVec[j]);
        }
        double tol = 1e-9 * d1Scale;
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (Math.abs(d1[i][j] - exitRate[i] * alphaVec[j]) > tol) {
                    return null;
                }
            }
        }
        return alpha;
    }

    // ------------------------------------------------------------------
    // Shared numerical helpers (package-private: reused by Rap_sample).
    // ------------------------------------------------------------------

    /** Copies a Matrix into a dense row-major array for fast element access. */
    static double[][] toArray2D(Matrix m) {
        int rows = m.getNumRows();
        int cols = m.getNumCols();
        double[][] out = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                out[i][j] = m.get(i, j);
            }
        }
        return out;
    }

    /** Row-vector times matrix product. */
    static double[] vecTimesMat(double[] v, double[][] m) {
        int rows = v.length;
        int cols = m[0].length;
        double[] out = new double[cols];
        for (int i = 0; i < rows; i++) {
            double vi = v[i];
            if (vi == 0.0) continue;
            double[] mi = m[i];
            for (int j = 0; j < cols; j++) {
                out[j] += vi * mi[j];
            }
        }
        return out;
    }

    /** Sum of the entries of a vector. */
    static double sum(double[] v) {
        double s = 0.0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
        }
        return s;
    }

    /** Inner product of two vectors of equal length. */
    static double dot(double[] x, double[] y) {
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += x[i] * y[i];
        }
        return s;
    }

    /** Clamps a probability to [0,1], absorbing roundoff of order 1e-16. */
    static double clampUnit(double p) {
        if (p < 0.0) return 0.0;
        if (p > 1.0) return 1.0;
        return p;
    }

    /**
     * Computes {@code w*expm(A*d)} without forming a new matrix exponential
     * whenever {@code ||A||*d} is small, by a scaled Taylor series applied to
     * the row vector. Falls back to {@code expm_higham} when the propagation
     * distance is too large for the substep budget.
     *
     * @param w        the row vector to propagate
     * @param a        dense copy of A
     * @param aMatrix  A as a Matrix (used only by the fallback path)
     * @param normA    the infinity norm of A
     * @param d        the propagation distance (must be non-negative)
     * @return the propagated row vector
     */
    static double[] expmPropagate(double[] w, double[][] a, Matrix aMatrix,
                                  double normA, double d) {
        int n = w.length;
        if (!(d > 0.0)) {
            return w.clone();
        }
        double theta = normA * d;
        int steps = 1;
        if (theta > TAYLOR_THETA) {
            steps = (int) Math.ceil(theta / TAYLOR_THETA);
            if (steps > TAYLOR_MAX_STEPS) {
                return vecTimesMat(w, toArray2D(aMatrix.scale(d).expm_higham()));
            }
        }
        double ds = d / steps;
        double[] acc = w.clone();
        for (int s = 0; s < steps; s++) {
            double[] term = acc.clone();
            double[] next = acc.clone();
            for (int k = 1; k <= TAYLOR_MAX_TERMS; k++) {
                term = vecTimesMat(term, a);
                double c = ds / k;
                double maxTerm = 0.0;
                double maxAcc = 0.0;
                for (int j = 0; j < n; j++) {
                    term[j] *= c;
                    next[j] += term[j];
                    double at = Math.abs(term[j]);
                    if (at > maxTerm) maxTerm = at;
                    double an = Math.abs(next[j]);
                    if (an > maxAcc) maxAcc = an;
                }
                if (maxTerm <= 1e-18 * Math.max(maxAcc, 1e-300)) {
                    break;
                }
            }
            acc = next;
        }
        return acc;
    }

    /**
     * Returns the dominant (largest real part) eigenvalue of A, which governs
     * the exponential decay of the tail. Falls back to {@code -1/mean} when the
     * eigenvalue decomposition is unavailable or does not yield a negative rate.
     *
     * @param A    the matrix whose spectrum is required
     * @param mean the distribution mean, used for the fallback rate
     * @return a negative decay rate
     */
    static double dominantRate(Matrix A, double mean) {
        double eta = Double.NEGATIVE_INFINITY;
        try {
            List<org.apache.commons.math3.complex.Complex> ev = A.eig();
            for (int i = 0; i < ev.size(); i++) {
                double re = ev.get(i).getReal();
                if (re > eta) {
                    eta = re;
                }
            }
        } catch (RuntimeException e) {
            eta = Double.NEGATIVE_INFINITY;
        }
        if (!(eta < 0.0) || Double.isInfinite(eta) || Double.isNaN(eta)) {
            eta = (mean > 0.0) ? -1.0 / mean : -1.0;
        }
        return eta;
    }
}
