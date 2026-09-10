/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Moment closures of the fluid drift, backing {@code options.method='minnormal'}
 * of SolverFluid. Java twin of the MATLAB {@code fluid_lld_scaling},
 * {@code fluid_min_closure}, {@code fluid_capacity_closure},
 * {@code fluid_share_closure} and {@code fluid_gps_share}.
 *
 * <p>The mean-field ODEs of SolverFluid close the moment hierarchy at first
 * order, replacing {@code E[min(X,c)]} by {@code min(E[X],c)} and the capacity
 * share {@code w_j X_j / sum_m w_m X_m} by that ratio at the mean. Both are
 * exact only where the map is locally linear, so their error peaks at the kink,
 * i.e. exactly where the two means meet (rho ~ 1 at a queueing station). The
 * functions here return instead the expectation under a normal marginal whose
 * variance is produced by the covariance (Lyapunov) equation, so mean and
 * covariance are solved self-consistently by SolverFluidMoments.</p>
 *
 * <p>Every closure collapses to the first-order expression when its variance
 * argument is zero, so callers share a single code path and the untouched
 * methods are bit-identical.</p>
 *
 * @see FluidMomentTerms
 * @see jline.solvers.fluid.analyzers.MinNormalAnalyzer
 */
public final class FluidClosures {

    private static final double SQRT2 = FastMath.sqrt(2.0);
    private static final double SQRT2PI = FastMath.sqrt(2.0 * FastMath.PI);

    private FluidClosures() {
    }

    /** A closure value together with its first two derivatives with respect to the mean. */
    public static final class ValueDeriv {
        public final double value;
        public final double deriv;
        public final double deriv2;

        public ValueDeriv(double value, double deriv) {
            this(value, deriv, 0.0);
        }

        public ValueDeriv(double value, double deriv, double deriv2) {
            this.value = value;
            this.deriv = deriv;
            this.deriv2 = deriv2;
        }
    }

    /**
     * A vector of capacity shares together with its Jacobian, and the joint-closure
     * covariance Cov(S_j,N) and its Jacobian. Each is null when not requested.
     */
    public static final class ShareResult {
        public final Matrix s;
        public final Matrix ds;
        public final Matrix cn;
        public final Matrix dcn;

        public ShareResult(Matrix s, Matrix ds) {
            this(s, ds, null, null);
        }

        public ShareResult(Matrix s, Matrix ds, Matrix cn, Matrix dcn) {
            this.s = s;
            this.ds = ds;
            this.cn = cn;
            this.dcn = dcn;
        }
    }

    /** Standard normal cumulative distribution function. */
    public static double normcdf(double z) {
        return 0.5 * Erf.erfc(-z / SQRT2);
    }

    /** Standard normal probability density function. */
    public static double normpdf(double z) {
        return FastMath.exp(-0.5 * z * z) / SQRT2PI;
    }

    /**
     * Limited load-dependent rate scaling at a CONTINUOUS population, mirroring
     * MATLAB {@code fluid_lld_scaling}.
     *
     * <p>{@code sn.lldscaling(i,:)} tabulates the station rate multiplier at
     * integer populations 1..lldlimit, and the discrete solvers read it as
     * {@code lldscaling(i, min(n, lldlimit))}. The fluid state is continuous, so
     * the table is read here by linear interpolation between consecutive
     * entries, clamped to the first entry below n=1 and to the last entry above
     * the table end, matching the clamping the CTMC already applies.</p>
     *
     * @param lldrow rate multipliers at populations 1..lldlimit; null or empty
     *               means no load dependence and returns a = 1, da = 0
     * @param n      population, may be non-integer or negative
     * @return alpha(n) and d alpha / d n, the latter zero on the clamped tails
     */
    public static ValueDeriv lldScaling(double[] lldrow, double n) {
        if (lldrow == null || lldrow.length == 0) {
            return new ValueDeriv(1.0, 0.0);
        }
        int L = lldrow.length;
        if (n <= 1) {
            return new ValueDeriv(lldrow[0], 0.0);
        }
        if (n >= L) {
            return new ValueDeriv(lldrow[L - 1], 0.0);
        }
        int k = (int) FastMath.floor(n);
        double frac = n - k;
        double a = lldrow[k - 1] * (1 - frac) + lldrow[k] * frac;
        return new ValueDeriv(a, lldrow[k] - lldrow[k - 1]);
    }

    /**
     * Min-normal moment closure of {@code E[min(X,Y)]} for jointly normal X, Y,
     * mirroring MATLAB {@code fluid_min_closure}.
     *
     * <p>This is the closure of Guenther, Stefanek and Bradley (EPEW/UKPEW 2012,
     * LNCS 7587:32-47, eq. 4), in their general two-population form</p>
     *
     * <pre>
     *   E[min(X,Y)] = E[X]*Phi((E[Y]-E[X])/th) + E[Y]*Phi((E[X]-E[Y])/th)
     *                 - th*phi((E[Y]-E[X])/th)
     *   th = (Var[X] - 2*Cov[X,Y] + Var[Y])^(1/2)
     * </pre>
     *
     * <p>and the derivative with respect to E[X] is {@code Phi((E[Y]-E[X])/th)},
     * the probability that X is the smaller of the two. SolverFluid only ever
     * needs the specialisation Y = c, the deterministic server count: vc = 0 and
     * covNc = 0 give th = sqrt(s2) and the expression collapses to the
     * truncated-normal form. The general arguments are kept so this IS the
     * published closure rather than one instance of it, and so a future
     * state-dependent capacity needs no new derivation.</p>
     *
     * <p>With th = 0 the expressions collapse to {@code min(n,c)} and to the
     * indicator {@code 1{n < c}}, recovering the first-order closure exactly.
     * The derivative at the kink is taken as 0, the right derivative of min(),
     * which is the convention already implied by the strict inequality test in
     * the closing rate factors.</p>
     *
     * @param n     mean of the first argument
     * @param c     mean of the second argument (server count when deterministic)
     * @param s2    variance of the first argument
     * @param vc    variance of the second argument (0 when deterministic)
     * @param covNc covariance of the two arguments
     * @return E[min(X,Y)] and P(X &lt; Y) under the normal marginal
     */
    public static ValueDeriv minClosure(double n, double c, double s2, double vc, double covNc) {
        double th2 = s2 - 2 * covNc + vc;
        if (th2 < 0) {
            th2 = 0; // a covariance beyond the Cauchy-Schwarz bound is not admissible
        }
        if (!(th2 > 0) || Double.isInfinite(c)) {
            // THE INDICATOR CARRIES A BAND, and it is a cross-codebase requirement rather
            // than a modelling choice. A saturated fluid fixed point sits exactly AT n = c,
            // and each engine's ODE stops on its own residual: MATLAB lands at 1.0004 and
            // the C++ port at 1 - 1.8e-13 on the same model, so a strict n < c reads
            // saturated in one and unsaturated in the other. That flips a whole Jacobian
            // row between zero and unit, and with it the hyperbolicity verdict that decides
            // whether SolverFLD answers with 'minnormal' or falls back to the first-order
            // method. A population within FineTol of the server count IS at the kink.
            boolean below = (c - n) > GlobalConstants.FineTol * FastMath.max(1.0, n);
            return new ValueDeriv(FastMath.min(n, c), below ? 1.0 : 0.0);
        }
        double th = FastMath.sqrt(th2);
        double d = (n - c) / th;
        double phid = normcdf(d);
        double dphid = normpdf(d);
        double h = n * (1 - phid) + c * phid - th * dphid;
        double dh = 1 - phid;
        // min() is piecewise linear, so its second derivative is carried entirely by
        // the atom at X = Y; smoothing over the normal marginal turns that atom into
        // the density. It stays zero on the degenerate branch, which smooths nothing.
        double d2h = -dphid / th;
        // A POPULATION IS NONNEGATIVE AND THE NORMAL MARGINAL IS NOT. For X >= 0
        // pathwise min(X,c) >= 0, and min() being concave Jensen puts
        // E[min(X,c)] <= min(E[X],c), so the value belongs to [0, min(n,c)]. The
        // normal marginal has no such support, and the mass it places below zero
        // drags the expectation out of that range once the mean falls to about one
        // standard deviation: n = 0, c = 1, th = 0.664 returns -0.019, which is a
        // station that CREATES work. Project onto the admissible range and take the
        // derivatives of the bound that binds, so the Jacobian still matches h.
        double hi = FastMath.min(n, c);
        boolean atLo = h < 0;
        boolean atHi = h > hi;
        h = FastMath.min(FastMath.max(h, 0.0), hi);
        if (atLo) {
            dh = 0.0;
        } else if (atHi) {
            dh = n < c ? 1.0 : 0.0;
        }
        if (atLo || atHi) {
            d2h = 0.0;
        }
        return new ValueDeriv(h, dh, d2h);
    }

    /** Specialisation of {@link #minClosure(double, double, double, double, double)} to a deterministic second argument. */
    public static ValueDeriv minClosure(double n, double c, double s2) {
        return minClosure(n, c, s2, 0.0, 0.0);
    }

    /**
     * Moment closure of the station capacity term psi(X) and of its derivative,
     * mirroring MATLAB {@code fluid_capacity_closure}.
     *
     * <p>Every scheduling branch of the closing rate factors scales the
     * coordinates of a station by psi(n_i)/n_i, where psi is how much work the
     * station clears at population n_i:</p>
     *
     * <pre>
     *   psi(n) = min(n,c) * alpha(n)   at a queueing station
     *   psi(n) = n        * alpha(n)   at an infinite server
     * </pre>
     *
     * <p>and alpha is the limited load-dependent scaling (1 when the station has
     * none). This returns E[psi(X)] for X ~ Normal(n, s2) together with
     * d/dn E[psi(X)]. s2 = 0 gives psi(n) itself, so the first-order closure is
     * the same code path.</p>
     *
     * <p>Without load dependence the expectation is the closed form of
     * {@link #minClosure}. With a tabulated alpha, {@link #lldScaling} makes
     * alpha piecewise linear on the integer lattice, so psi is piecewise
     * QUADRATIC with breakpoints at the integers and at c, and the expectation
     * is integrated segment by segment against the normal density using the
     * truncated moments M0, M1, M2. The derivative is E[psi'(X)] by
     * differentiation under the integral sign, valid because psi is Lipschitz.</p>
     *
     * <p>The integration must be exact, NOT a fixed quadrature rule.
     * Gauss-Hermite with fixed nodes applied to a piecewise-linear integrand
     * does not smooth its kinks, it relocates them: the resulting estimate of
     * E[psi] is itself piecewise linear in n, so its second derivative is zero
     * almost everywhere and the 1/N refinement silently returns a null
     * correction. That is what the segment-wise closed form below avoids.</p>
     *
     * <p>psi is extended by zero below n = 0, the only physically admissible
     * continuation, since a station holding no jobs clears no work.</p>
     *
     * @param n      mean population at the station
     * @param c      number of servers (ignored when isInf)
     * @param s2     population variance, 0 for the first-order closure
     * @param lldrow load-dependent scaling, null or empty when absent
     * @param isInf  true at an infinite-server station
     * @return E[psi(X)] and d/dn E[psi(X)]
     */
    public static ValueDeriv capacityClosure(double n, double c, double s2, double[] lldrow, boolean isInf) {
        if (lldrow == null || lldrow.length == 0) {
            if (isInf) {
                return new ValueDeriv(n, 1.0);
            }
            ValueDeriv m = minClosure(n, c, s2);
            // the normal marginal puts mass below zero, where min(X,c) = X < 0, and a
            // negative capacity destroys mass through the non-negativity clamp
            if (m.value < 0) {
                return new ValueDeriv(0.0, 0.0, 0.0);
            }
            return m;
        }

        int L = lldrow.length;
        if (!(s2 > 0)) {
            return psiAt(n, c, lldrow, isInf);
        }

        double s = FastMath.sqrt(s2);
        // breakpoints of psi: the lattice of the table, the origin, and the saturation point c
        int nbp = L + 1;
        boolean addC = !isInf && !Double.isInfinite(c) && c > L;
        double[] bps = new double[addC ? nbp + 1 : nbp];
        for (int k = 0; k <= L; k++) {
            bps[k] = k;
        }
        if (addC) {
            bps[nbp] = c;
        }

        double h = 0;
        double dh = 0;
        double d2h = 0;
        // psi is extended by zero below the first breakpoint, so psi' jumps there too
        // and that atom belongs in psi'' exactly like the interior ones
        double bPrev = 0;
        double cPrev = 0;
        for (int k = 0; k < bps.length; k++) {
            double p = bps[k];
            double q = (k < bps.length - 1) ? bps[k + 1] : Double.POSITIVE_INFINITY;
            double[] abc = segmentCoefficients(p, q, c, lldrow, L, isInf);
            double[] mom = truncatedMoments(p, q, n, s);
            h += abc[0] * mom[0] + abc[1] * mom[1] + abc[2] * mom[2];
            dh += abc[1] * mom[0] + 2 * abc[2] * mom[1];
            d2h += 2 * abc[2] * mom[0];
            // psi' jumps across this breakpoint, so psi'' carries an atom there; the
            // segment sum above sees only the quadratic part and would miss it
            double jump = (abc[1] + 2 * abc[2] * p) - (bPrev + 2 * cPrev * p);
            d2h += jump * normpdf((p - n) / s) / s;
            bPrev = abc[1];
            cPrev = abc[2];
        }
        return new ValueDeriv(h, dh, d2h);
    }

    /** psi(u) = A + B*u + C*u^2 on [p,q], from base(u)*alpha(u). */
    private static double[] segmentCoefficients(double p, double q, double c, double[] lldrow, int L, boolean isInf) {
        double a0;
        double a1;
        // alpha is constant on [0,1] and on [L,inf), linear on each unit interval
        if (p >= L) {
            a0 = lldrow[L - 1];
            a1 = 0;
        } else if (p < 1) {
            a0 = lldrow[0];
            a1 = 0;
        } else {
            int k = (int) FastMath.floor(p);
            a1 = lldrow[k] - lldrow[k - 1];
            a0 = lldrow[k - 1] - a1 * k;
        }
        // the breakpoints guarantee the segment lies entirely on one side of the saturation point
        if (isInf || Double.isInfinite(c) || q <= c) {
            return new double[]{0, a0, a1};
        }
        return new double[]{c * a0, c * a1, 0};
    }

    /** Truncated moments E[X^j * 1{p &lt; X &lt; q}] for X ~ Normal(n, s^2). */
    private static double[] truncatedMoments(double p, double q, double n, double s) {
        double zp = (p - n) / s;
        double pp = normpdf(zp);
        double pP = normcdf(zp);
        double qP;
        double pq;
        if (Double.isInfinite(q)) {
            qP = 1.0;
            pq = 0.0;
        } else {
            qP = normcdf((q - n) / s);
            pq = normpdf((q - n) / s);
        }
        double m0 = qP - pP;
        double m1 = n * m0 + s * (pp - pq);
        double m2;
        if (Double.isInfinite(q)) {
            m2 = (n * n + s * s) * m0 + s * ((p + n) * pp);
        } else {
            m2 = (n * n + s * s) * m0 + s * ((p + n) * pp - (q + n) * pq);
        }
        return new double[]{m0, m1, m2};
    }

    /** psi and its derivative at a continuous population, extended by zero below 0. */
    private static ValueDeriv psiAt(double u, double c, double[] lldrow, boolean isInf) {
        if (u <= 0) {
            return new ValueDeriv(0.0, 0.0);
        }
        ValueDeriv a = lldScaling(lldrow, u);
        double base = isInf ? u : FastMath.min(u, c);
        double dbase = isInf ? 1.0 : (u < c ? 1.0 : 0.0);
        // base and alpha are both piecewise linear, so psi is piecewise quadratic and
        // psi'' = 2*base'*alpha' away from the breakpoints; the atoms AT them are not
        // representable without a marginal, and this first-order branch smooths none.
        return new ValueDeriv(base * a.value, dbase * a.value + base * a.deriv,
                2.0 * dbase * a.deriv);
    }

    /**
     * How much of the second-order correction the series admits at this point, and
     * d tau / d ratio, mirroring MATLAB {@code local_expansion_weight}.
     *
     * <p>Every second-order term of the share closure is a term of the series for
     * E[1/v], whose successive terms are in the ratio {@code Var(v)/v^2}, so the
     * truncation is meaningful below 1 and the terms GROW above it. Nothing in the
     * algebra notices: at a near-empty station the corrections come back larger than
     * the quantity they correct, and the drift that follows is not integrable.</p>
     *
     * <p>One on [0,1], zero from 4 up, and the C^1 smoothstep between. Both ends
     * matter. The lower one has to be EXACTLY one on the whole convergent region, so
     * every model already inside it is bit-identical; the upper one has to be reached
     * with a vanishing derivative, because the drift is integrated and a kink in it is
     * what collapses the step size. The thresholds are the series, not a tuning: at
     * ratio 1 successive terms stop shrinking, and at ratio 4 the standard deviation of
     * v is twice its mean, where a non-negative v has essentially no mass near the
     * point being expanded about.</p>
     *
     * @param ratio Var(v)/v^2 at the point being expanded about
     * @return {tau, d tau / d ratio}
     */
    public static double[] expansionWeight(double ratio) {
        double lo = 1.0;
        double hi = 4.0;
        if (ratio <= lo) {
            return new double[] {1.0, 0.0};
        }
        if (ratio >= hi) {
            return new double[] {0.0, 0.0};
        }
        double t = (ratio - lo) / (hi - lo);
        return new double[] {1.0 - t * t * (3.0 - 2.0 * t), -6.0 * t * (1.0 - t) / (hi - lo)};
    }

    /**
     * Project a jointly closed per-coordinate service share onto the set it has to
     * live in, mirroring MATLAB {@code local_project_rate}: {@code r >= 0},
     * {@code r <= xb} where that bound applies, and {@code sum(r) = tot}.
     *
     * <p>THE JOINT CLOSURE IS AN EXPANSION AND CAN LEAVE THAT SET. {@code r = s*psi +
     * psi'*Cov(S,N)} adds a term that sums to ZERO over the coordinates, so it moves
     * mass between them and its entries can push one past either bound; the
     * first-order share {@code x_j/n_i*psi} cannot, being x_j scaled by
     * {@code psi/n_i <= 1}. Either breach ends the same way, because the integrator
     * holds every coordinate non-negative: {@code r_j > x_j} drains coordinate j
     * faster than it holds, the state goes negative and the clamp INJECTS mass.</p>
     *
     * <p>THE UPPER BOUND HOLDS ONLY WITHOUT LOAD DEPENDENCE, which is what CAPPED
     * selects. r is an expected NUMBER in service so {@code r_j <= x_j}, but
     * {@code psi(n) = min(n,c)*alpha(n)} folds the load-dependent scaling into the
     * same variable, and with alpha &gt; 1 the first-order share itself exceeds
     * x_j.</p>
     *
     * <p>Clip, then move the residual onto the coordinates that still have slack in
     * proportion to it, so {@code sum(r) = tot} survives and the station still clears
     * what its capacity closure says it clears. It is a NO-OP whenever the expansion
     * stayed inside the set, which is why models already inside it are
     * bit-identical.</p>
     *
     * @param r      the jointly closed rates, projected in place
     * @param xb     the coordinate means, the upper bound when capped
     * @param capped true without load dependence, where r_j &lt;= x_j holds
     * @param tot    what the capacity closure says the station clears
     */
    public static void projectRate(double[] r, double[] xb, boolean capped, double tot) {
        double zt = GlobalConstants.Zero;
        boolean inside = true;
        for (int j = 0; j < r.length && inside; j++) {
            if (!(r[j] >= -zt) || (capped && !(r[j] <= xb[j] + zt))) {
                inside = false;
            }
        }
        if (inside) {
            return;
        }
        for (int j = 0; j < r.length; j++) {
            r[j] = FastMath.max(r[j], 0.0);
            if (capped) {
                r[j] = FastMath.min(r[j], xb[j]);
            }
        }
        for (int it = 0; it <= r.length; it++) {
            double sum = 0;
            for (int j = 0; j < r.length; j++) {
                sum += r[j];
            }
            double d = tot - sum;
            if (FastMath.abs(d) <= zt) {
                break;
            }
            double[] slack = new double[r.length];
            double tsl = 0;
            for (int j = 0; j < r.length; j++) {
                slack[j] = d > 0 ? (capped ? xb[j] - r[j] : 1.0) : r[j];
                tsl += slack[j];
            }
            if (tsl <= zt) {
                break;
            }
            for (int j = 0; j < r.length; j++) {
                r[j] = FastMath.max(r[j] + d * slack[j] / tsl, 0.0);
                if (capped) {
                    r[j] = FastMath.min(r[j], xb[j]);
                }
            }
        }
    }

    /** The share closure without the joint-closure covariance. */
    public static ShareResult shareClosure(Matrix x, Matrix wv, Matrix cov, boolean wantJac) {
        return shareClosure(x, wv, cov, wantJac, false);
    }

    /**
     * Second-order closure of the capacity share of a sharing discipline,
     * mirroring MATLAB {@code fluid_share_closure}.
     *
     * <p>A DPS station gives coordinate j the fraction
     * {@code S_j = w_j X_j / sum_m w_m X_m} of its capacity, and PS is the same
     * expression with unit weights. The first-order fluid closure evaluates that
     * ratio at the mean, which is not E[S_j]: the map is a ratio, so Jensen's
     * inequality biases it towards the coordinates carrying the LARGER weight.
     * Writing u_j = w_j X_j and v = sum_m u_m, the delta method gives</p>
     *
     * <pre>
     *   E[u_j/v] = mu_j/v - Cov(u_j,v)/v^2 + mu_j*Var(v)/v^3 + O(sigma^3)
     * </pre>
     *
     * <p>evaluated at the means. The correction is exactly capacity-conserving:
     * {@code sum_j Cov(u_j,v) = Cov(v,v) = Var(v)}, so the two correction terms
     * cancel in the sum and {@code sum_j S_j = 1} identically, as it must for a
     * work-conserving discipline. That identity is the invariant to check on any
     * change here.</p>
     *
     * <p>The expansion is local and fails when v is small against its own
     * standard deviation, where the exact expectation is a Cauchy-like integral
     * with no finite mean. There a raw share can come out negative; it is
     * clipped at zero and the survivors renormalised, which preserves the
     * conservation identity.</p>
     *
     * @param x       coordinate means of one station block (n x 1)
     * @param wv      per-coordinate weight, constant within a class (n x 1)
     * @param cov     covariance of the same coordinates (n x n); null or zero
     *                selects the first-order (plug-in) share
     * @param wantJac whether to build the (n x n) Jacobian ds_j/dx_m
     * @return the expected shares, summing to one, and their Jacobian
     */
    public static ShareResult shareClosure(Matrix x, Matrix wv, Matrix cov, boolean wantJac,
            boolean wantCov) {
        int n = x.getNumRows();
        // DENSE, because both blocks are filled ENTRY BY ENTRY and both are full.
        // `new Matrix(r,c)` is sparse (CSC) and its set() is O(nnz), so filling an
        // n-by-n Jacobian this way is O(n^4) and its get() is a binary search. This
        // runs inside the fluid ODE's right-hand side, once per step and once per
        // Jacobian column, which is where LN+Fluid on test_LQN_5 spent its hours.
        Matrix s = Matrix.dense(n, 1);
        Matrix ds = wantJac ? Matrix.dense(n, n) : null;
        Matrix cn = wantCov ? Matrix.dense(n, 1) : null;
        Matrix dcn = (wantCov && wantJac) ? Matrix.dense(n, n) : null;

        double[] u = new double[n];
        double v = 0;
        for (int j = 0; j < n; j++) {
            u[j] = wv.get(j, 0) * x.get(j, 0);
            v += u[j];
        }
        if (v <= 0) {
            return new ShareResult(s, ds, cn, dcn);
        }

        for (int j = 0; j < n; j++) {
            s.set(j, 0, u[j] / v);
        }
        if (wantJac) {
            for (int j = 0; j < n; j++) {
                for (int m = 0; m < n; m++) {
                    double val = (j == m ? wv.get(j, 0) / v : 0.0) - (u[j] / (v * v)) * wv.get(m, 0);
                    ds.set(j, m, val);
                }
            }
        }

        if (cov == null || cov.getNumRows() != n || !hasNonZero(cov)) {
            return new ShareResult(s, ds, cn, dcn);
        }

        double[] cuv = new double[n]; // Cov(u_j, v)
        double cvv = 0;               // Var(v)
        for (int j = 0; j < n; j++) {
            double acc = 0;
            for (int m = 0; m < n; m++) {
                acc += cov.get(j, m) * wv.get(m, 0);
            }
            cuv[j] = wv.get(j, 0) * acc;
            cvv += wv.get(j, 0) * acc;
        }

        double v2 = v * v;
        double v3 = v2 * v;
        double v4 = v3 * v;

        // How far into the series this point sits, and how much of the second-order
        // correction that leaves admissible. tau depends on x through v alone, since
        // cov is held fixed here exactly as the Jacobians are.
        double[] tw = expansionWeight(cvv / v2);
        double tau = tw[0];
        if (tau <= 0 && !wantJac) {
            return new ShareResult(s, ds, cn, dcn); // first-order share, already exact
        }
        double dtauScale = tw[1] * (-2.0 * cvv / v3); // d tau / d x_m = this * wv_m

        if (wantCov) {
            // Cov(S_j, N) at the means, from grad(S_j)'*cov*1: S_j = w_j X_j / V, so
            // dS_j/dX_a = w_j*delta_aj/V - u_j*w_a/V^2 and the two pieces contract
            // against cov*1 = Cov(X,N) and w'*cov*1 = Cov(V,N).
            double[] an = new double[n];
            double cvn = 0;
            for (int j = 0; j < n; j++) {
                double acc = 0;
                for (int m = 0; m < n; m++) {
                    acc += cov.get(j, m);
                }
                an[j] = acc;
                cvn += wv.get(j, 0) * acc;
            }
            double[] cn0 = new double[n];
            for (int j = 0; j < n; j++) {
                cn0[j] = wv.get(j, 0) * an[j] / v - u[j] * (cvn / v2);
                cn.set(j, 0, tau * cn0[j]);
            }
            if (dcn != null) {
                for (int j = 0; j < n; j++) {
                    for (int m = 0; m < n; m++) {
                        double val = tau * (-(wv.get(j, 0) * an[j]) * (wv.get(m, 0) / v2)
                                - (j == m ? (cvn / v2) * wv.get(j, 0) : 0.0)
                                + (2.0 * cvn / v3) * u[j] * wv.get(m, 0))
                                + cn0[j] * dtauScale * wv.get(m, 0);
                        dcn.set(j, m, val);
                    }
                }
            }
        }

        double[] scorr = new double[n];
        for (int j = 0; j < n; j++) {
            scorr[j] = -cuv[j] / v2 + (u[j] * cvv) / v3;
            s.set(j, 0, s.get(j, 0) + tau * scorr[j]);
        }
        if (wantJac) {
            for (int j = 0; j < n; j++) {
                for (int m = 0; m < n; m++) {
                    double val = ds.get(j, m)
                            + tau * ((2.0 / v3) * cuv[j] * wv.get(m, 0)
                                    + (j == m ? (cvv / v3) * wv.get(j, 0) : 0.0)
                                    - (3.0 * cvv / v4) * u[j] * wv.get(m, 0))
                            + scorr[j] * dtauScale * wv.get(m, 0);
                    ds.set(j, m, val);
                }
            }
        }

        // A COORDINATE CARRYING NO MASS MUST NOT DECIDE THE CLIP. With u_j = 0 the
        // plug-in share is zero and the correction leaves s_j = -Cov(u_j,v)/v^2, a
        // quantity of the order of rounding whose SIGN is not meaningful. Letting it
        // select the branch below zeroes that coordinate's whole Jacobian row, and a
        // zero row is an exact zero eigenvalue: the Lyapunov step then reads the fixed
        // point as non-hyperbolic and SolverFLD silently drops from 'minnormal' to the
        // first-order method. Clip only a share that is negative BEYOND the numerical
        // zero.
        boolean allNonNegative = true;
        for (int j = 0; j < n; j++) {
            if (s.get(j, 0) < -GlobalConstants.Zero) {
                allNonNegative = false;
                break;
            }
        }
        if (allNonNegative) {
            for (int j = 0; j < n; j++) {
                if (s.get(j, 0) < 0) {
                    s.set(j, 0, 0);
                }
            }
            return new ShareResult(s, ds, cn, dcn);
        }

        // the expansion has left its region of validity for at least one coordinate;
        // clip and renormalise the survivors so the shares still sum to one
        boolean[] act = new boolean[n];
        double total = 0;
        boolean anyActive = false;
        for (int j = 0; j < n; j++) {
            act[j] = s.get(j, 0) > 0;
            if (act[j]) {
                total += s.get(j, 0);
                anyActive = true;
            }
        }
        if (!anyActive) {
            for (int j = 0; j < n; j++) {
                s.set(j, 0, u[j] / v);
            }
            if (wantJac) {
                for (int j = 0; j < n; j++) {
                    for (int m = 0; m < n; m++) {
                        ds.set(j, m, (j == m ? wv.get(j, 0) / v : 0.0) - (u[j] / (v * v)) * wv.get(m, 0));
                    }
                }
            }
            return new ShareResult(s, ds, cn, dcn);
        }

        Matrix sNew = Matrix.dense(n, 1);
        Matrix dsNew = wantJac ? Matrix.dense(n, n) : null;
        double[] dT = null;
        if (wantJac) {
            dT = new double[n];
            for (int m = 0; m < n; m++) {
                double acc = 0;
                for (int j = 0; j < n; j++) {
                    if (act[j]) {
                        acc += ds.get(j, m);
                    }
                }
                dT[m] = acc;
            }
        }
        for (int j = 0; j < n; j++) {
            if (!act[j]) {
                continue;
            }
            sNew.set(j, 0, s.get(j, 0) / total);
            if (wantJac) {
                for (int m = 0; m < n; m++) {
                    dsNew.set(j, m, ds.get(j, m) / total - (s.get(j, 0) / (total * total)) * dT[m]);
                }
            }
        }
        return new ShareResult(sNew, dsNew, cn, dcn);
    }

    /**
     * Expected capacity share of a GPS station under a normal marginal,
     * mirroring MATLAB {@code fluid_gps_share}.
     *
     * <p>GPS divides the server by WEIGHT among the classes that are BACKLOGGED,
     * and then equally among that class's own jobs (State.afterEventStation, GPS
     * branch: {@code cir = min(nir,1)}, share {@code = w_r/(w*cir')}). The share
     * therefore depends on the backlog INDICATOR vector, not on the populations,
     * and that is what makes GPS unreachable for a first-order closure: with
     * continuous x_k &gt; 0 every class is always backlogged, the indicator is
     * identically one, and the share collapses to the constant
     * {@code w_r/sum_j w_j} regardless of load. That constant is the
     * heavy-traffic limit and is wrong everywhere else, so for GPS the second
     * moment is not a correction, it is the entire mechanism.</p>
     *
     * <p>The closure is an EXACT enumeration rather than an expansion. The share
     * is piecewise CONSTANT over the 2^K backlog patterns, so</p>
     *
     * <pre>
     *   E[S_r] = sum_{A contains r} P(backlog set = A) * w_r / sum_{j in A} w_j
     * </pre>
     *
     * <p>carries no truncation error once the pattern probabilities are given.
     * Those come from the marginals, {@code P(N_k >= 1) = Phi((x_k - 1/2)/sigma_k)}
     * with the continuity correction for an integer population, multiplied as if
     * the backlogs were independent. That independence is the one approximation
     * here and it is not innocuous: in a closed network the station coordinates
     * are NEGATIVELY correlated through population conservation, so the exact
     * treatment would need multivariate-normal orthant probabilities.</p>
     *
     * <p>The empty pattern contributes zero share, so
     * {@code sum_r E[S_r] = 1 - P(all classes empty)} rather than 1. That is
     * deliberate and is how the idle server is represented: GPS is single-server,
     * so the backlog indicator plays the role that min(n,c) plays at a PS
     * station, and no separate capacity term is applied by the caller.</p>
     *
     * <p>Unlike the DPS ratio closure the expansion is NOT perturbative in
     * sigma: as sigma_k -&gt; 0 the probability tends to a step at x_k = 1 and its
     * derivative phi(.)/sigma_k diverges, so the Jacobian stiffens at low
     * variance. vk = 0 falls back to the hard indicator with zero derivative,
     * which is the correct mean-field starting point for the outer iteration.</p>
     *
     * @param xk      per-class populations at the station (K x 1)
     * @param wk      per-class GPS weights, normalised internally (K x 1)
     * @param vk      per-class population variances; 0 selects the hard indicator
     * @param wantJac whether to build the (K x K) Jacobian ds_r/dx_m
     * @return the expected shares, summing to 1 - P(station empty), and their Jacobian
     */
    public static ShareResult gpsShare(Matrix xk, Matrix wk, Matrix vk, boolean wantJac) {
        int K = xk.getNumRows();
        Matrix s = Matrix.dense(K, 1);
        Matrix ds = wantJac ? Matrix.dense(K, K) : null;

        // the enumeration is 2^K, so refuse rather than crawl; K here is the number
        // of classes AT ONE STATION, which is small in every practical model
        if (K > 12) {
            line_error(mfilename(new Object() {
            }), String.format("GPS closes its capacity share by enumerating the 2^K backlog patterns of a station, "
                    + "and this station carries %d classes. Above 12 the enumeration is no longer tractable; use "
                    + "options.method=\"closing\" with a DPS station instead.", K));
        }

        double sw = 0;
        for (int k = 0; k < K; k++) {
            sw += wk.get(k, 0);
        }
        if (sw <= 0) {
            return new ShareResult(s, ds);
        }
        double[] w = new double[K];
        for (int k = 0; k < K; k++) {
            w[k] = wk.get(k, 0) / sw;
        }

        // P(class k backlogged) = P(N_k >= 1) for an INTEGER population, so the normal
        // approximation needs the continuity correction P(N_k > 1/2)
        double[] p = new double[K];
        double[] dp = new double[K];
        for (int k = 0; k < K; k++) {
            if (vk.get(k, 0) > 0) {
                double sd = FastMath.sqrt(vk.get(k, 0));
                double z = (xk.get(k, 0) - 0.5) / sd;
                p[k] = normcdf(z);
                dp[k] = normpdf(z) / sd;
            } else {
                p[k] = xk.get(k, 0) > 0 ? 1.0 : 0.0;
                dp[k] = 0.0;
            }
        }

        double[][] dsdp = wantJac ? new double[K][K] : null;
        int npat = 1 << K;
        for (int mask = 1; mask < npat; mask++) {
            double W = 0;
            for (int k = 0; k < K; k++) {
                if ((mask & (1 << k)) != 0) {
                    W += w[k];
                }
            }
            if (W <= 0) {
                continue; // every backlogged class in this pattern carries zero weight
            }
            double[] q = new double[K];
            double prod = 1.0;
            for (int k = 0; k < K; k++) {
                q[k] = ((mask & (1 << k)) != 0) ? p[k] : 1 - p[k];
                prod *= q[k];
            }
            for (int k = 0; k < K; k++) {
                if ((mask & (1 << k)) == 0) {
                    continue;
                }
                s.set(k, 0, s.get(k, 0) + prod * w[k] / W);
            }
            if (!wantJac) {
                continue;
            }
            for (int m = 0; m < K; m++) {
                double prodm = 1.0; // product over j != m
                for (int j = 0; j < K; j++) {
                    if (j != m) {
                        prodm *= q[j];
                    }
                }
                double sign = ((mask & (1 << m)) != 0) ? 1.0 : -1.0;
                for (int k = 0; k < K; k++) {
                    if ((mask & (1 << k)) != 0) {
                        dsdp[k][m] += sign * prodm * w[k] / W;
                    }
                }
            }
        }

        if (wantJac) {
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < K; m++) {
                    ds.set(k, m, dsdp[k][m] * dp[m]); // chain rule through p_m
                }
            }
        }
        return new ShareResult(s, ds);
    }

    private static boolean hasNonZero(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (m.get(i, j) != 0.0) {
                    return true;
                }
            }
        }
        return false;
    }
}
