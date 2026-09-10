/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.complex.Complex;

import static jline.io.InputOutput.line_warning;
import static jline.lib.butools.ph.CheckMERepresentation.checkMERepresentation;

/**
 * A Matrix Exponential (ME) distribution.
 *
 * ME distributions are characterized by an initial vector alpha and a matrix parameter A.
 * They generalize Phase-Type (PH) distributions by allowing alpha to have entries outside [0,1]
 * and A to have arbitrary structure (not necessarily a valid sub-generator).
 *
 * Representation:
 * - alpha: initial vector (may have negative entries or sum != 1)
 * - A: matrix parameter (must have all eigenvalues with negative real parts)
 * - Dominant eigenvalue of A must be negative and real
 *
 * The distribution function is: F(t) = 1 - alpha * exp(A*t) * e (under certain conditions)
 * Moments: m_k = k! * alpha * (-A)^(-k) * e
 */
public class ME extends Markovian {

    /**
     * Constants of the direct density scan, see {@link #scanNegativeDensity}.
     * They match MATLAB ME.m and native Python so that the three codebases reach
     * the same verdict.
     */
    private static final double ME_SCAN_TAIL = 1e-12;        // residual mass beyond the horizon
    private static final double ME_SCAN_HORIZON_CAP = 1e4;   // cap for near-degenerate A
    private static final int ME_SCAN_PER_PERIOD = 20;        // samples per period of oscillation
    private static final int ME_SCAN_MIN_PTS = 2001;
    private static final int ME_SCAN_MAX_PTS = 200001;
    private static final double ME_SCAN_RELTOL = 1e-10;      // negative only below -reltol*max|f|

    /**
     * Outcome of the density scan: whether a negative value was found, and where.
     */
    public static class NegativeDensityScan {
        private final boolean negative;
        private final double fmin;
        private final double tmin;

        NegativeDensityScan(boolean negative, double fmin, double tmin) {
            this.negative = negative;
            this.fmin = fmin;
            this.tmin = tmin;
        }

        /** @return true if the scan found a point where the density is negative */
        public boolean isNegative() { return negative; }

        /** @return the smallest density value seen */
        public double getFmin() { return fmin; }

        /** @return the time at which the smallest density value was seen */
        public double getTmin() { return tmin; }
    }

    /**
     * Formats a double the way the C and MATLAB %g conversion does, so that the
     * warning text is identical in the three codebases. Java's %g keeps trailing
     * zeros, which %g in MATLAB and Python does not.
     *
     * @param x the value to format
     * @return the shortest %g-style rendering of x
     */
    private static String fmtG(double x) {
        String s = String.format("%.6g", x);
        if (s.indexOf('.') < 0) {
            return s;
        }
        int e = s.indexOf('e');
        if (e < 0) {
            e = s.indexOf('E');
        }
        String mant = (e < 0) ? s : s.substring(0, e);
        String suf = (e < 0) ? "" : s.substring(e);
        while (mant.endsWith("0")) {
            mant = mant.substring(0, mant.length() - 1);
        }
        if (mant.endsWith(".")) {
            mant = mant.substring(0, mant.length() - 1);
        }
        return mant + suf;
    }

    /**
     * Searches the density f(t) = -alpha*expm(A*t)*A*e for a negative value.
     *
     * <p>A negative value found here is a witness: it proves that the
     * representation is not a distribution. Finding none proves nothing, so the
     * caller must not report the converse.</p>
     *
     * <p>This is used in place of
     * {@link jline.lib.butools.ph.CheckMEPositiveDensity} as the trigger for the
     * construction-time warning. That routine searches for a Markovian
     * monocyclic equivalent, which is a sufficient condition only, and its
     * verdict depends on the representation rather than on the distribution: for
     * {@code alpha = [1,0,0]}, {@code A = [[-0.5,0,0],[0,-1,w],[0,-w,-1]]} the
     * distribution is Exp(0.5) for every w, yet the search fails once
     * w &gt;= 2*pi. It also costs of the order of a second per call at search
     * order 1000, which is far too slow for a constructor.</p>
     *
     * <p>The horizon covers all but {@code ME_SCAN_TAIL} of the mass, using the
     * dominant (least negative) eigenvalue of A; the sampling rate resolves the
     * fastest oscillation present, taken from the largest imaginary part.</p>
     *
     * @param alphaRow the initial vector, as a 1 x n row
     * @param A the matrix parameter
     * @return the scan outcome
     */
    public static NegativeDensityScan scanNegativeDensity(Matrix alphaRow, Matrix A) {
        int n = A.getNumRows();

        double decay = Double.NEGATIVE_INFINITY;
        double wmax = 0.0;
        for (Complex ev : A.eig()) {
            if (ev.getReal() > decay) {
                decay = ev.getReal();
            }
            double im = Math.abs(ev.getImaginary());
            if (im > wmax) {
                wmax = im;
            }
        }
        if (!(decay < 0.0) || Double.isNaN(decay) || Double.isInfinite(decay)) {
            // Not a valid ME; checkMERepresentation has already rejected it.
            return new NegativeDensityScan(false, 0.0, 0.0);
        }
        double horizon = Math.min(-Math.log(ME_SCAN_TAIL) / Math.abs(decay), ME_SCAN_HORIZON_CAP);

        int npts = ME_SCAN_MIN_PTS;
        if (wmax > 0.0) {
            long wanted = (long) Math.ceil(horizon * wmax * ME_SCAN_PER_PERIOD / (2.0 * Math.PI)) + 1L;
            if (wanted > npts) {
                npts = (int) Math.min(wanted, (long) ME_SCAN_MAX_PTS);
            }
        }
        npts = Math.min(npts, ME_SCAN_MAX_PTS);

        double step = horizon / (npts - 1);
        // One matrix exponential, then propagate: v_k = alpha*expm(A*k*step).
        Matrix E = A.copy().scale(step).expm();
        Matrix ve = A.mult(Matrix.ones(n, 1)).scale(-1);

        Matrix v = alphaRow.copy();
        double fmin = Double.POSITIVE_INFINITY;
        double tmin = 0.0;
        double fabsmax = 0.0;
        for (int k = 0; k < npts; k++) {
            double f = v.mult(ve).get(0, 0);
            if (f < fmin) {
                fmin = f;
                tmin = k * step;
            }
            if (Math.abs(f) > fabsmax) {
                fabsmax = Math.abs(f);
            }
            v = v.mult(E);
        }

        return new NegativeDensityScan(fmin < -ME_SCAN_RELTOL * fabsmax, fmin, tmin);
    }

    /**
     * Creates a Matrix Exponential distribution with specified initial vector and matrix parameter.
     *
     * @param alpha the initial vector (row vector as Matrix)
     * @param A the matrix parameter (must be square with negative real eigenvalues)
     * @throws IllegalArgumentException if the representation is invalid
     */
    public ME(Matrix alpha, Matrix A) {
        this(alpha, A, true);
    }

    /**
     * Creates a Matrix Exponential distribution, optionally skipping the density scan.
     *
     * @param alpha the initial vector (row vector as Matrix)
     * @param A the matrix parameter (must be square with negative real eigenvalues)
     * @param checkDensity scan the density for a negative value. Set to false only by
     *        subclasses whose representation is a density by construction, such as
     *        {@link CME}, where the scan cannot fire and costs O(1e5) propagations of a
     *        large matrix.
     * @throws IllegalArgumentException if the representation is invalid
     */
    protected ME(Matrix alpha, Matrix A, boolean checkDensity) {
        super("ME", 2);

        // Ensure alpha is a row vector (1 x n)
        // If passed as column vector (n x 1), transpose it
        Matrix alphaRow = alpha;
        if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
            alphaRow = alpha.transpose();
        }

        // Strict validation using BuTools
        if (!checkMERepresentation(alphaRow, A, 1e-14)) {
            throw new IllegalArgumentException("Invalid ME representation: " +
                    "Check that A is square, alpha and A have compatible dimensions, " +
                    "all eigenvalues of A have negative real parts, and " +
                    "the dominant eigenvalue is real.");
        }

        // see _kb/01-model-classes.md (ME/CME sections) for rationale
        NegativeDensityScan scan = checkDensity
                ? scanNegativeDensity(alphaRow, A)
                : new NegativeDensityScan(false, 0.0, 0.0);
        if (scan.isNegative()) {
            // see _kb/01-model-classes.md (Java process-construction notes) for rationale
            line_warning("ME",
                    "The ME representation has a negative density: f(t) = -alpha*expm(A*t)*A*e "
                            + "reaches %s at t = %s. Moments and transforms remain well defined, "
                            + "but evalPDF returns negative values and sample() will not reproduce "
                            + "a proper distribution.",
                    fmtG(scan.getFmin()), fmtG(scan.getTmin()));
        }

        nPhases = alphaRow.getNumElements();

        // Store parameters (always store as row vector)
        this.setParam(1, "alpha", alphaRow);
        this.setParam(2, "A", A);

        // Build MatrixCell process representation for compatibility with map_* functions
        // Process format: {D0=A, D1=-A*e*alpha'}
        MatrixCell rep = new MatrixCell();
        rep.set(0, A);  // D0 = A

        // D1 = -A * e * alpha' = outer product of (-A*e) and alpha
        // where e is column vector of ones
        Matrix ones = Matrix.ones(nPhases, 1);
        Matrix Ae = A.mult(ones).scale(-1);  // -A * e (column vector, n x 1)
        Matrix D1 = Ae.mult(alphaRow);  // outer product: (n x 1) * (1 x n) = (n x n)
        rep.set(1, D1);

        this.setProcess(rep);
    }

    /**
     * Gets the initial vector alpha.
     *
     * @return the initial vector as a Matrix
     */
    public Matrix getAlpha() {
        return (Matrix) this.getParam(1).getValue();
    }

    /**
     * Gets the matrix parameter A.
     *
     * @return the matrix parameter A
     */
    public Matrix getA() {
        return (Matrix) this.getParam(2).getValue();
    }

    @Override
    public long getNumberOfPhases() {
        return nPhases;
    }

    @Override
    public MatrixCell getProcess() {
        return this.process;
    }

    /**
     * Generates samples from this ME distribution.
     *
     * <p>Sampling inverts the exact CDF F(t) = 1 - alpha*exp(A*t)*e. The CTMC
     * walk inherited from {@link Markovian} (map_sample) is not used because it
     * presumes a phase-type interpretation of (alpha, A), which fails whenever
     * alpha has negative entries or A has negative off-diagonal entries.</p>
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of n i.i.d. samples
     */
    @Override
    public double[] sample(int n, java.util.Random random) {
        java.util.Random rng = (random != null)
                ? random : jline.util.RandomManager.getThreadRandomAsRandom();
        return jline.api.mam.Me_sample.me_sample(getAlpha(), getA(), n, rng);
    }

    // see _kb/01-model-classes.md (ME/RAP: getMu/getPhi inheritance rationale) for rationale

    /**
     * Gets the mean, m1 = -alpha*A^(-1)*e.
     *
     * <p>Overrides the Markovian formula 1/map_lambda, which obtains the rate from the
     * stationary vector of D0+D1. That vector is a probabilistic object of a Markovian
     * process, and computing it for an ME means solving a linear system whose
     * conditioning degrades with the oscillation of A: a CME of order 101 came out with a
     * relative error of 4e-8, where the definition below is exact to 1e-13. Native Python
     * ME.getMean already uses the definition.</p>
     *
     * @return the mean
     */
    @Override
    public double getMean() {
        Matrix alpha = getAlpha();
        Matrix A = getA();
        Matrix e = Matrix.ones(A.getNumRows(), 1);
        return -alpha.mult(A.inv()).mult(e).value();
    }

    /**
     * Gets the variance, m2 - m1^2 with m2 = 2*alpha*A^(-2)*e.
     *
     * @return the variance
     */
    @Override
    public double getVar() {
        Matrix alpha = getAlpha();
        Matrix A = getA();
        Matrix e = Matrix.ones(A.getNumRows(), 1);
        Matrix Ainv = A.inv();
        double m1 = -alpha.mult(Ainv).mult(e).value();
        double m2 = 2.0 * alpha.mult(Ainv).mult(Ainv).mult(e).value();
        return m2 - m1 * m1;
    }

    /**
     * Gets the squared coefficient of variation, var/mean^2.
     *
     * @return the squared coefficient of variation
     */
    @Override
    public double getSCV() {
        double m1 = getMean();
        return getVar() / (m1 * m1);
    }

    /**
     * Creates an ME distribution by fitting the given moments.
     * Uses BuTools MEFromMoments algorithm.
     *
     * @param moments array of moments (requires 2*M-1 moments for order M ME distribution)
     * @return an ME distribution matching the given moments
     * @throws IllegalArgumentException if moments are invalid or fitting fails
     */
    public static ME fitMoments(double[] moments) {
        // Use the MEFromMoments algorithm to fit the distribution
        jline.lib.butools.ph.MERepresentation rep = jline.lib.butools.ph.MEFromMoments.meFromMoments(moments);
        return new ME(rep.getAlpha(), rep.getA());
    }

    /**
     * Creates an ME distribution from an exponential distribution.
     * This is a convenience method showing that Exp is a special case of ME.
     *
     * @param rate the rate parameter (lambda)
     * @return an ME distribution equivalent to Exp(rate)
     */
    public static ME fromExp(double rate) {
        Matrix alpha = new Matrix(new double[]{1.0});
        Matrix A = new Matrix(new double[][]{{-rate}});
        return new ME(alpha, A);
    }

    /**
     * Creates an ME distribution from an Erlang distribution.
     * This is a convenience method showing that Erlang is a special case of ME.
     *
     * @param k number of phases
     * @param rate rate parameter for each phase
     * @return an ME distribution equivalent to Erlang(k, rate)
     */
    public static ME fromErlang(int k, double rate) {
        Matrix alpha = new Matrix(1, k);
        alpha.set(0, 0, 1.0);  // alpha = [1, 0, 0, ..., 0]

        Matrix A = new Matrix(k, k);
        for (int i = 0; i < k; i++) {
            A.set(i, i, -rate);  // diagonal
            if (i < k - 1) {
                A.set(i, i + 1, rate);  // super-diagonal
            }
        }

        return new ME(alpha, A);
    }

    /**
     * Creates an ME distribution from a HyperExponential distribution.
     * This is a convenience method showing that HyperExp is a special case of ME.
     *
     * @param p array of probabilities for each branch
     * @param rates array of rates for each branch
     * @return an ME distribution equivalent to HyperExp(p, rates)
     */
    public static ME fromHyperExp(double[] p, double[] rates) {
        if (p.length != rates.length) {
            throw new IllegalArgumentException("p and rates must have the same length");
        }

        int k = p.length;
        // Create alpha as row vector (1 x k), not column vector
        Matrix alpha = new Matrix(1, k);
        for (int i = 0; i < k; i++) {
            alpha.set(0, i, p[i]);
        }

        Matrix A = new Matrix(k, k);
        for (int i = 0; i < k; i++) {
            A.set(i, i, -rates[i]);  // diagonal matrix of rates
        }

        return new ME(alpha, A);
    }
}
