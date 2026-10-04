/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Per-coordinate service share of the closing fluid ODE and its analytic
 * Jacobian. Java twin of the MATLAB {@code ode_rates_closing_factors} and of
 * the rate-factor part of {@code fluid_drift_jacobian}.
 *
 * <p>The two live in one class on purpose. The drift the ODE integrates and the
 * Jacobian the covariance (Lyapunov) equation reads must agree branch by branch,
 * or the fixed point of the mean solve is not the point at which the covariance
 * is linearised; keeping them in separate files is what makes that divergence
 * possible. {@link jline.solvers.fluid.handlers.PassageTimeODE} and
 * {@link FluidMomentTerms} both delegate here, so the closing method and the
 * moment-closure methods evaluate literally the same shares.</p>
 *
 * <p>{@code sigma2} is the per-station closure variance (0: first-order
 * closure), {@code lldscaling} is {@code sn.lldscaling} (null: no load
 * dependence) and {@code covblk} carries the per-station coordinate covariance
 * blocks that close the capacity-share ratio at second order. All three leave
 * the legacy code path bit-identical when absent, so the untouched methods are
 * unaffected.</p>
 *
 * @see FluidClosures
 * @see FluidMomentTerms
 */
public final class FluidRateFactors {

    /** Same stability margin the Lyapunov gate uses, so both read one scale. */
    private static final double SQRT_EPS = FastMath.sqrt(2.220446049250313e-16);

    private final int M;
    private final int K;
    private final boolean[][] enabled;
    private final Matrix qIndices;
    private final Matrix Kic;
    private final Matrix nservers;
    private final Matrix w;
    private final SchedStrategy[] sched;
    private final Matrix lldscaling;

    /**
     * @param M          number of stations
     * @param K          number of classes
     * @param enabled    per-(station,class) service flag
     * @param qIndices   per-(station,class) first state coordinate (0-based)
     * @param Kic        per-(station,class) phase count
     * @param nservers   per-station server count
     * @param w          per-(station,class) DPS/GPS weight, 1 elsewhere
     * @param sched      per-station scheduling strategy
     * @param lldscaling {@code sn.lldscaling}, null when the model has none
     */
    public FluidRateFactors(int M, int K, boolean[][] enabled, Matrix qIndices, Matrix Kic,
                            Matrix nservers, Matrix w, SchedStrategy[] sched, Matrix lldscaling) {
        this.M = M;
        this.K = K;
        this.enabled = enabled;
        this.qIndices = qIndices;
        this.Kic = Kic;
        this.nservers = nservers;
        this.w = w;
        this.sched = sched;
        this.lldscaling = (lldscaling == null || lldscaling.isEmpty()) ? null : lldscaling;
    }

    /**
     * Scheduling disciplines with a branch below. Anything else would keep
     * {@code rates = x}, i.e. be integrated as an INFINITE SERVER, and the
     * answer would be wrong without any warning: on Delay(Z=1) -&gt; Queue(c=1),
     * N=4, exact Q2 = 3.0154, the fall-through returns 2.0000.
     */
    public static boolean isSupported(SchedStrategy s) {
        return s == SchedStrategy.INF || s == SchedStrategy.EXT || s == SchedStrategy.PS
                || s == SchedStrategy.FCFS || s == SchedStrategy.DPS || s == SchedStrategy.GPS;
    }

    /**
     * Refuses any station whose discipline has no drift branch, naming it.
     *
     * <p>Called before the drift is built rather than left to a fall-through:
     * the metric readers accept SIRO as FCFS, so an unguarded SIRO station was
     * integrated as INF while its metrics were read as if it shared the server,
     * and every closing-family method returned a wrong answer silently.</p>
     */
    public void checkSupported() {
        for (int i = 0; i < M; i++) {
            if (!isSupported(sched[i])) {
                line_error(mfilename(new Object() {
                }), String.format("Station %d uses %s, which has no fluid drift branch in the closing family. "
                        + "Use options.method=\"matrix\", which builds a PS drift for every queueing station.",
                        i + 1, sched[i]));
            }
        }
    }

    /** Per-station load-dependence row, null when the station has none. */
    private double[] lldRow(int i) {
        if (lldscaling == null || i >= lldscaling.getNumRows()) {
            return null;
        }
        int L = lldscaling.getNumCols();
        double[] row = new double[L];
        boolean allOne = true;
        for (int j = 0; j < L; j++) {
            row[j] = lldscaling.get(i, j);
            if (row[j] != 1.0) {
                allOne = false;
            }
        }
        return allOne ? null : row; // a station without load dependence keeps the plain branch
    }

    private static Matrix covOf(Matrix[] covblk, int i) {
        if (covblk == null || i >= covblk.length) {
            return null;
        }
        return covblk[i];
    }

    private int blockLo(int i) {
        return (int) qIndices.get(i, 0);
    }

    private int blockHi(int i) {
        return (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
    }

    /**
     * First station whose population sits ON the saturation kink n_i = c_i of
     * the first-order rate factor, or -1 when none does.
     *
     * <p>With sigma2 = 0 the occupancy factor is min(n_i, c_i), whose derivative
     * is the indicator of the unsaturated region: slope 1 below c_i, slope 0
     * above, and NO derivative at c_i itself. {@link #jacobian} resolves the tie
     * with a strict {@code n_i > c_i}, as do the MATLAB and Python twins, so it
     * silently returns the left slope there. That one-sided value is not the a.e.
     * derivative, and which side a fixed point lands on is decided by the
     * integrator's rounding residue rather than by the model: the same network at
     * N=4 converges to n_i = c_i exactly here and to c_i + 3e-6 in Python, which
     * flipped a Jacobian eigenvalue between -0.5 and 0 and hence the hyperbolicity
     * verdict. Callers that need a differentiable drift must consult this instead
     * of trusting the tie-break.</p>
     *
     * <p>Only the branches that take the indicator derivative can sit on a kink:
     * a positive sigma2 or a load-dependent row makes the closure smooth, and an
     * infinite server never saturates.</p>
     *
     * @param x      phase-resolved state
     * @param sigma2 per-station population variance, null or all-zero for the
     *               first-order closure
     * @return the station index, or -1 when every station is off its kink
     */
    public int driftKinkStation(double[] x, double[] sigma2) {
        int[] all = driftKinkStations(x, sigma2);
        return all.length == 0 ? -1 : all[0];
    }

    /**
     * Every station sitting on the saturation kink, empty when none does. Same
     * test as {@link #driftKinkStation(double[], double[])}, which is its first
     * element.
     *
     * @param x      phase-resolved state
     * @param sigma2 per-station population variance, null or all-zero for the
     *               first-order closure
     * @return the station indices, in increasing order
     */
    public int[] driftKinkStations(double[] x, double[] sigma2) {
        if (anyPositive(sigma2)) {
            return new int[0]; // the Gaussian closure is smooth, it has no kink
        }
        java.util.List<Integer> hits = new java.util.ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (sched[i] == SchedStrategy.INF || sched[i] == SchedStrategy.EXT) {
                continue;
            }
            if (lldRow(i) != null) {
                continue; // psi is piecewise quadratic and closed smoothly
            }
            double c = nservers.get(i, 0);
            if (Double.isInfinite(c) || c <= 0) {
                continue;
            }
            double ni = sum(x, blockLo(i), blockHi(i));
            if (ni <= 0) {
                continue; // g = x on an empty station, no saturation term
            }
            if (FastMath.abs(ni - c) <= SQRT_EPS * FastMath.max(1.0, c)) {
                hits.add(Integer.valueOf(i));
            }
        }
        int[] out = new int[hits.size()];
        for (int t = 0; t < out.length; t++) {
            out[t] = hits.get(t).intValue();
        }
        return out;
    }

    /**
     * A copy of {@code x} with every station sitting on its saturation kink moved
     * to {@code c_i*(1+rel)}, i.e. strictly onto one side of it.
     *
     * <p>The point of the copy is that the kink is a measure-zero event whose
     * TIE-BREAK is arbitrary while its two one-sided Jacobians are both perfectly
     * well defined. Evaluating the drift at the two nudged points is how a caller
     * asks whether the side matters: if both sides agree on the verdict it wants,
     * the tie-break is inconsequential and the reference's strict
     * {@code n_i > c_i} may be used; if they disagree, the answer would be decided
     * by the integrator's rounding residue instead of by the model.</p>
     *
     * <p>The station's coordinates are scaled together, so the phase mix and every
     * other station are untouched.</p>
     *
     * @param x   phase-resolved state
     * @param rel signed relative offset from the server count, e.g. 1e-6
     * @return the nudged state, or {@code x} itself when no station is on a kink
     */
    public double[] nudgedOffKink(double[] x, double[] sigma2, double rel) {
        int[] kinks = driftKinkStations(x, sigma2);
        if (kinks.length == 0) {
            return x;
        }
        double[] y = x.clone();
        for (int t = 0; t < kinks.length; t++) {
            int i = kinks[t];
            int lo = blockLo(i);
            int hi = blockHi(i);
            double ni = sum(y, lo, hi);
            if (!(ni > 0)) {
                continue;
            }
            double scale = nservers.get(i, 0) * (1.0 + rel) / ni;
            for (int j = lo; j < hi; j++) {
                y[j] *= scale;
            }
        }
        return y;
    }

    /**
     * Per-coordinate service share, before event indexing and before the
     * constant rate factors are applied.
     *
     * @param x       phase-resolved state
     * @param sigma2  per-station population variance, null or zero for the first-order closure
     * @param covblk  per-station coordinate covariance blocks, null for the plug-in share
     * @return one share per state coordinate
     */
    public Matrix factors(double[] x, double[] sigma2, Matrix[] covblk) {
        boolean gaussian = anyPositive(sigma2);

        // DENSE: every entry is written below, and a CSC set() is O(nnz). This
        // vector is rebuilt on every fluid ODE right-hand-side evaluation.
        Matrix rates = Matrix.dense(x.length, 1);
        for (int j = 0; j < x.length; j++) {
            rates.set(j, 0, x[j]); // basic vector valid for INF and for PS with min(ni,c) = ni
        }

        for (int i = 0; i < M; i++) {
            double[] lldrow = lldRow(i);
            double s2i = (sigma2 == null) ? 0.0 : sigma2[i];
            switch (sched[i]) {
                case INF: {
                    // each job is served at its own rate and the share is the identity;
                    // alpha(n_i) scales the whole station
                    if (lldrow != null) {
                        int lo = blockLo(i);
                        int hi = blockHi(i);
                        double ni = sum(x, lo, hi);
                        if (ni > 0) {
                            double h = FluidClosures.capacityClosure(ni, nservers.get(i, 0), s2i, lldrow, true).value;
                            for (int idx = lo; idx < hi; idx++) {
                                rates.set(idx, 0, x[idx] / ni * h);
                            }
                        }
                    }
                    break;
                }
                case EXT: {
                    // a delay except that the local population must conserve unit mass
                    for (int k = 0; k < K; k++) {
                        if (!enabled[i][k]) {
                            continue;
                        }
                        int lo = (int) qIndices.get(i, k);
                        int hi = lo + (int) Kic.get(i, k);
                        rates.set(lo, 0, 1 - sum(x, lo + 1, hi));
                    }
                    break;
                }
                case PS:
                case FCFS: {
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    int n = hi - lo;
                    double ni = sum(x, lo, hi);
                    if ((gaussian || lldrow != null) && ni > 0) {
                        double h = FluidClosures.capacityClosure(ni, nservers.get(i, 0), s2i, lldrow, false).value;
                        Matrix ci = covOf(covblk, i);
                        if (ci == null) {
                            for (int idx = lo; idx < hi; idx++) {
                                rates.set(idx, 0, x[idx] / ni * h);
                            }
                        } else {
                            // THE SHARE AND THE CAPACITY ARE CLOSED JOINTLY. What the
                            // station clears is S_j*psi(N), and both factors move with N,
                            // so the product needs Cov(S_j,N)*psi'(n) on top of the two
                            // separate closures; see FluidClosures#shareClosure. With unit
                            // weights this is the DPS branch below.
                            Matrix xb = block(x, lo, hi);
                            Matrix wv = ones(n);
                            FluidClosures.ValueDeriv cap =
                                    FluidClosures.capacityClosure(ni, nservers.get(i, 0), s2i, lldrow, false);
                            FluidClosures.ShareResult sr =
                                    FluidClosures.shareClosure(xb, wv, ci, false, true);
                            double[] rb = new double[n];
                            double[] xbv = new double[n];
                            for (int idx = 0; idx < n; idx++) {
                                rb[idx] = sr.s.get(idx, 0) * h + cap.deriv * sr.cn.get(idx, 0);
                                xbv[idx] = xb.get(idx, 0);
                            }
                            FluidClosures.projectRate(rb, xbv, lldrow == null || lldrow.length == 0, h);
                            for (int idx = lo; idx < hi; idx++) {
                                rates.set(idx, 0, rb[idx - lo]);
                            }
                        }
                    } else if (ni > nservers.get(i, 0)) { // case min = ni handled by rates = x
                        for (int idx = lo; idx < hi; idx++) {
                            rates.set(idx, 0, x[idx] / ni * nservers.get(i, 0));
                        }
                    }
                    break;
                }
                case DPS: {
                    // DPS is PS with a weighted share: the class-k coordinates get
                    // w_k*x/ni of the station capacity psi(xi) instead of x/xi of it
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    int n = hi - lo;
                    double[] wi = normalisedWeights(i);
                    Matrix wv = Matrix.dense(n, 1);
                    for (int k = 0; k < K; k++) {
                        if (!enabled[i][k]) {
                            continue;
                        }
                        int klo = (int) qIndices.get(i, k);
                        int khi = klo + (int) Kic.get(i, k);
                        for (int idx = klo; idx < khi; idx++) {
                            wv.set(idx - lo, 0, wi[k]);
                        }
                    }
                    double xi = sum(x, lo, hi);
                    double wx = 0;
                    for (int idx = lo; idx < hi; idx++) {
                        wx += wv.get(idx - lo, 0) * x[idx];
                    }
                    if (xi > 0 && wx > 0) {
                        FluidClosures.ValueDeriv cap =
                                FluidClosures.capacityClosure(xi, nservers.get(i, 0), s2i, lldrow, false);
                        Matrix xb = block(x, lo, hi);
                        FluidClosures.ShareResult sr =
                                FluidClosures.shareClosure(xb, wv, covOf(covblk, i), false, true);
                        int nb = hi - lo;
                        double[] rb = new double[nb];
                        double[] xbv = new double[nb];
                        for (int idx = 0; idx < nb; idx++) {
                            rb[idx] = sr.s.get(idx, 0) * cap.value + cap.deriv * sr.cn.get(idx, 0);
                            xbv[idx] = xb.get(idx, 0);
                        }
                        FluidClosures.projectRate(rb, xbv, lldrow == null || lldrow.length == 0, cap.value);
                        for (int idx = lo; idx < hi; idx++) {
                            rates.set(idx, 0, rb[idx - lo]);
                        }
                    }
                    break;
                }
                case GPS: {
                    // GPS splits the server by WEIGHT among the BACKLOGGED classes, then
                    // equally among that class's own jobs; no capacity term multiplies it
                    if (nservers.get(i, 0) > 1) {
                        line_error(mfilename(new Object() {
                        }), "Multi-server GPS stations are not supported yet.");
                    }
                    ClassMoments cm = classMoments(x, i, covOf(covblk, i));
                    Matrix sk = FluidClosures.gpsShare(cm.xk, weightColumn(i), cm.vk, false).s;
                    double a = 1;
                    if (lldrow != null) {
                        a = FluidClosures.lldScaling(lldrow, sum(x, blockLo(i), blockHi(i))).value;
                    }
                    for (int k = 0; k < K; k++) {
                        if (!enabled[i][k] || cm.xk.get(k, 0) <= 0) {
                            continue;
                        }
                        int klo = (int) qIndices.get(i, k);
                        int khi = klo + (int) Kic.get(i, k);
                        for (int idx = klo; idx < khi; idx++) {
                            rates.set(idx, 0, x[idx] / cm.xk.get(k, 0) * sk.get(k, 0) * a);
                        }
                    }
                    break;
                }
                default:
                    line_error(mfilename(new Object() {
                    }), String.format("Station %d uses %s, which has no fluid drift branch in the closing family.",
                            i + 1, sched[i]));
            }
        }
        return rates;
    }

    /**
     * Analytic Jacobian dg/dx of the rate factors, mirroring {@link #factors}
     * branch by branch.
     *
     * <p>With sigma2 = 0 the derivative of the occupancy factor is the indicator
     * of the unsaturated region, i.e. the a.e. derivative of the first-order
     * closure. With sigma2 &gt; 0 it is the smooth derivative of the Gaussian
     * closure.</p>
     *
     * @param x      phase-resolved state
     * @param sigma2 per-station population variance
     * @param covblk per-station coordinate covariance blocks
     * @return the (n x n) Jacobian of the rate factors
     */
    public Matrix jacobian(double[] x, double[] sigma2, Matrix[] covblk) {
        boolean gaussian = anyPositive(sigma2);
        int n = x.length;
        // DENSE: the diagonal alone fills it, and the per-station blocks below are
        // full, so element-wise insertion into a CSC would be O(n^4).
        Matrix G = Matrix.dense(n, n);
        for (int j = 0; j < n; j++) {
            G.set(j, j, 1.0); // INF, EXT phases 2..end, and every policy without a case below
        }

        for (int i = 0; i < M; i++) {
            double[] lldrow = lldRow(i);
            double s2i = (sigma2 == null) ? 0.0 : sigma2[i];
            switch (sched[i]) {
                case INF: {
                    if (lldrow == null) {
                        break; // g = x, identity already in place
                    }
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    double ni = sum(x, lo, hi);
                    if (ni <= 0) {
                        break;
                    }
                    FluidClosures.ValueDeriv vd =
                            FluidClosures.capacityClosure(ni, nservers.get(i, 0), s2i, lldrow, true);
                    setScaledShareJacobian(G, x, lo, hi, vd.value, vd.deriv, ni);
                    break;
                }
                case EXT: {
                    for (int k = 0; k < K; k++) {
                        if (!enabled[i][k]) {
                            continue;
                        }
                        int lo = (int) qIndices.get(i, k);
                        int hi = lo + (int) Kic.get(i, k);
                        zeroRow(G, lo);
                        for (int m = lo + 1; m < hi; m++) {
                            G.set(lo, m, -1.0);
                        }
                    }
                    break;
                }
                case PS:
                case FCFS: {
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    double ni = sum(x, lo, hi);
                    if (ni <= 0) {
                        break; // g = x on an empty station
                    }
                    double h;
                    double dh;
                    if (gaussian || lldrow != null) {
                        FluidClosures.ValueDeriv vd =
                                FluidClosures.capacityClosure(ni, nservers.get(i, 0), s2i, lldrow, false);
                        Matrix ci = covOf(covblk, i);
                        if (ci != null) {
                            // the share closure is reached only from here, exactly as in factors()
                            int nb = hi - lo;
                            FluidClosures.ShareResult sr =
                                    FluidClosures.shareClosure(block(x, lo, hi), ones(nb), ci, true, true);
                            setShareJacobian(G, lo, hi, sr, vd.value, vd.deriv, vd.deriv2);
                            break;
                        }
                        h = vd.value;
                        dh = vd.deriv;
                    } else if (ni > nservers.get(i, 0)
                            - GlobalConstants.FineTol * FastMath.max(1.0, ni)) {
                        // THE SATURATION TEST CARRIES A BAND, and it is a cross-codebase
                        // requirement: a saturated fixed point sits exactly at ni = c, and
                        // each engine's ODE stops on its own residual (MATLAB 1.0004, the
                        // C++ port 1 - 1.8e-13 on the same model). A strict ni > c reads
                        // saturated in one and unsaturated in the other, which flips this
                        // whole station block between a zero row and the identity, and with
                        // it the hyperbolicity verdict. The two-sided driftKinkStation guard
                        // nudges by 1e-6, two orders wider, so it still resolves both sides.
                        h = nservers.get(i, 0);
                        dh = 0;
                    } else {
                        break; // g = x, identity already in place
                    }
                    setScaledShareJacobian(G, x, lo, hi, h, dh, ni);
                    break;
                }
                case DPS: {
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    int nb = hi - lo;
                    double[] wi = normalisedWeights(i);
                    Matrix wv = Matrix.dense(nb, 1);
                    for (int k = 0; k < K; k++) {
                        if (!enabled[i][k]) {
                            continue;
                        }
                        int klo = (int) qIndices.get(i, k);
                        int khi = klo + (int) Kic.get(i, k);
                        for (int idx = klo; idx < khi; idx++) {
                            wv.set(idx - lo, 0, wi[k]);
                        }
                    }
                    double xi = sum(x, lo, hi);
                    double wx = 0;
                    for (int idx = lo; idx < hi; idx++) {
                        wx += wv.get(idx - lo, 0) * x[idx];
                    }
                    if (xi <= 0 || wx <= 0) {
                        break; // g = x on an empty station
                    }
                    FluidClosures.ValueDeriv vd =
                            FluidClosures.capacityClosure(xi, nservers.get(i, 0), s2i, lldrow, false);
                    FluidClosures.ShareResult sr =
                            FluidClosures.shareClosure(block(x, lo, hi), wv, covOf(covblk, i), true, true);
                    setShareJacobian(G, lo, hi, sr, vd.value, vd.deriv, vd.deriv2);
                    break;
                }
                case GPS: {
                    if (nservers.get(i, 0) > 1) {
                        line_error(mfilename(new Object() {
                        }), "Multi-server GPS stations are not supported yet.");
                    }
                    int lo = blockLo(i);
                    int hi = blockHi(i);
                    ClassMoments cm = classMoments(x, i, covOf(covblk, i));
                    FluidClosures.ShareResult sr =
                            FluidClosures.gpsShare(cm.xk, weightColumn(i), cm.vk, true);
                    double a = 1;
                    double da = 0;
                    if (lldrow != null) {
                        FluidClosures.ValueDeriv vd = FluidClosures.lldScaling(lldrow, sum(x, lo, hi));
                        a = vd.value;
                        da = vd.deriv;
                    }
                    for (int r = lo; r < hi; r++) {
                        zeroRow(G, r);
                    }
                    for (int k = 0; k < K; k++) {
                        double xkk = cm.xk.get(k, 0);
                        if (!enabled[i][k] || xkk <= 0) {
                            continue;
                        }
                        int klo = (int) qIndices.get(i, k);
                        int khi = klo + (int) Kic.get(i, k);
                        double sa = sr.s.get(k, 0) * a;
                        for (int r = klo; r < khi; r++) {
                            for (int c = klo; c < khi; c++) {
                                double add = (r == c ? sa / xkk : 0.0) - (sa / (xkk * xkk)) * x[r];
                                G.set(r, c, G.get(r, c) + add);
                            }
                            for (int m = 0; m < K; m++) {
                                if (!enabled[i][m]) {
                                    continue;
                                }
                                int mlo = (int) qIndices.get(i, m);
                                int mhi = mlo + (int) Kic.get(i, m);
                                double add = (a * sr.ds.get(k, m) / xkk) * x[r];
                                for (int c = mlo; c < mhi; c++) {
                                    G.set(r, c, G.get(r, c) + add);
                                }
                            }
                            if (da != 0) {
                                double add = (sr.s.get(k, 0) * da / xkk) * x[r];
                                for (int c = lo; c < hi; c++) {
                                    G.set(r, c, G.get(r, c) + add);
                                }
                            }
                        }
                    }
                    break;
                }
                default:
                    line_error(mfilename(new Object() {
                    }), String.format("Station %d uses %s, which has no fluid drift branch in the closing family.",
                            i + 1, sched[i]));
            }
        }
        return G;
    }

    /** g_j = x_j*h(ni)/ni, so dg_j/dx_m = delta_jm*f + x_j*f' with f = h/ni. */
    private static void setScaledShareJacobian(Matrix G, double[] x, int lo, int hi, double h, double dh, double ni) {
        double f = h / ni;
        double fp = (ni * dh - h) / (ni * ni);
        for (int r = lo; r < hi; r++) {
            zeroRow(G, r);
            for (int c = lo; c < hi; c++) {
                G.set(r, c, (r == c ? f : 0.0) + x[r] * fp);
            }
        }
    }

    /** g = s(x_blk)*psi(xi), so dg_j/dx_m = ds_j/dx_m*psi + s_j*dpsi. */
    /**
     * g = s(x_blk)*psi(ni) + psi'(ni)*cn(x_blk), the joint closure of the share and
     * the capacity; differentiating it with the covariance held fixed adds psi'*dcn
     * and psi''*cn to the product rule.
     */
    private static void setShareJacobian(Matrix G, int lo, int hi, FluidClosures.ShareResult sr,
                                         double psi, double dpsi, double d2psi) {
        for (int r = lo; r < hi; r++) {
            zeroRow(G, r);
            for (int c = lo; c < hi; c++) {
                double val = sr.ds.get(r - lo, c - lo) * psi + sr.s.get(r - lo, 0) * dpsi;
                if (sr.cn != null) {
                    val += dpsi * sr.dcn.get(r - lo, c - lo) + sr.cn.get(r - lo, 0) * d2psi;
                }
                G.set(r, c, val);
            }
        }
    }

    private static void zeroRow(Matrix G, int r) {
        for (int c = 0; c < G.getNumCols(); c++) {
            G.set(r, c, 0.0);
        }
    }

    /** Per-class population and variance at a station, plus its coordinate blocks. */
    private static final class ClassMoments {
        final Matrix xk;
        final Matrix vk;

        ClassMoments(Matrix xk, Matrix vk) {
            this.xk = xk;
            this.vk = vk;
        }
    }

    private ClassMoments classMoments(double[] x, int i, Matrix ci) {
        Matrix xk = Matrix.dense(K, 1);
        Matrix vk = Matrix.dense(K, 1);
        int lo = blockLo(i);
        for (int k = 0; k < K; k++) {
            if (!enabled[i][k]) {
                continue;
            }
            int klo = (int) qIndices.get(i, k);
            int khi = klo + (int) Kic.get(i, k);
            xk.set(k, 0, sum(x, klo, khi));
            if (ci != null) {
                double acc = 0;
                for (int r = klo; r < khi; r++) {
                    for (int c = klo; c < khi; c++) {
                        acc += ci.get(r - lo, c - lo);
                    }
                }
                vk.set(k, 0, FastMath.max(0, acc));
            }
        }
        return new ClassMoments(xk, vk);
    }

    private double[] normalisedWeights(int i) {
        double[] wi = new double[K];
        double sw = 0;
        for (int k = 0; k < K; k++) {
            wi[k] = w.get(i, k);
            sw += wi[k];
        }
        if (sw > 0) {
            for (int k = 0; k < K; k++) {
                wi[k] /= sw;
            }
        }
        return wi;
    }

    private Matrix weightColumn(int i) {
        Matrix wk = Matrix.dense(K, 1);
        for (int k = 0; k < K; k++) {
            wk.set(k, 0, w.get(i, k));
        }
        return wk;
    }

    private static double sum(double[] x, int lo, int hi) {
        double acc = 0;
        for (int j = lo; j < hi; j++) {
            acc += x[j];
        }
        return acc;
    }

    private static Matrix block(double[] x, int lo, int hi) {
        Matrix b = Matrix.dense(hi - lo, 1);
        for (int j = lo; j < hi; j++) {
            b.set(j - lo, 0, x[j]);
        }
        return b;
    }

    private static Matrix ones(int n) {
        Matrix o = Matrix.dense(n, 1);
        for (int j = 0; j < n; j++) {
            o.set(j, 0, 1.0);
        }
        return o;
    }

    private static boolean anyPositive(double[] v) {
        if (v == null) {
            return false;
        }
        for (int j = 0; j < v.length; j++) {
            if (v[j] > 0) {
                return true;
            }
        }
        return false;
    }
}
