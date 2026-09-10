/**
 * @file Explicit closed-form normalizing constant of a multiclass limited load-dependent network
 *
 * Port at parity of matlab/src/api/pfqn/pfqn_explicit_ld.m, i.e. Theorem 1 of
 * G. Casale, P. G. Harrison, W. H. Ong, "Facilitating Load-Dependent Queueing Analysis
 * Through Factorization", Perform. Eval. 2021, carried over the divided-difference form
 * of G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic
 * Methods", ACM SIGMETRICS 2017, Corollary 3.2.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.nc.Pfqn_explicit;
import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Pfqn_explicit_ld {
    private Pfqn_explicit_ld() {}

    /** Result of the closed form, mirroring [lG, G, method, lossDigits]. */
    public static final class Result {
        /** Logarithm of the normalizing constant. */
        public final double lG;
        /** The normalizing constant. */
        public final double G;
        /** Expression used for g_sigma, "distinct" (Eq. 15) or "repeated" (Eq. 16). */
        public final String method;
        /** Decimal digits lost to cancellation. */
        public final double lossDigits;

        public Result(double lG, double G, String method, double lossDigits) {
            this.lG = lG;
            this.G = G;
            this.method = method;
            this.lossDigits = lossDigits;
        }
    }

    /**
     * Explicit closed-form normalizing constant of a multiclass limited load-dependent network.
     *
     * <p>Load-dependent counterpart of Pfqn_explicit. It evaluates the same divided-difference
     * form of Casale (SIGMETRICS 2017), Corollary 3.2,
     *
     * <pre>
     *   G(N) = sum_{0&lt;=t&lt;=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r C(N_r,t_r) h_t(|N|)
     * </pre>
     *
     * but substitutes for the single-class constant h_t(|N|) the LIMITED LOAD-DEPENDENT closed
     * form of Casale, Harrison and Ong (Perform. Eval. 2021), Theorem 1, Eq. (8),
     *
     * <pre>
     *   h_theta(N) = sum_{0&lt;=v&lt;s} g_sigma(N-|v|) prod_k phi_k(v_k)
     *   phi_k(v_k) = theta_k^v_k / prod_{t=1..v_k} alpha_k(t) * (1 - alpha_k(v_k)/alpha_k(s_k))
     * </pre>
     *
     * at the induced demands theta_k(t) = sum_r t_r L(k,r). Here alpha_k(.) = mu(k,.) is the
     * load-dependent scaling of station k, s_k the population past which it stays constant,
     * sigma_k = theta_k/alpha_k(s_k) the SCALED demands, and g_sigma the FIXED-RATE single-class
     * constant at those scaled demands, which is exactly what Pfqn_explicit evaluates in closed
     * form (Eqs. 15 and 16). The result is therefore explicit throughout, with no recursion over
     * population.</p>
     *
     * <p>Two conventions of Theorem 1 are not those of the equilibrium distribution and are easy
     * to get wrong. alpha_k(0) is taken as ZERO inside the bracket of phi_k, so that
     * phi_k(0) = 1, even though the state probabilities use alpha_k(0) = 1; and
     * g_sigma(n) = 0 for n &lt; 0, which caps the outer sum at |v| &lt;= |N|. With
     * alpha_k(n) = min(n,s_k) the expression collapses to Gordon's multi-server formula,
     * Oper. Res. 38(5), 1990, Eq. (29), but unlike that one it needs neither a multi-server shape
     * nor distinct scaled demands.</p>
     *
     * <p>LIMITED LOAD DEPENDENCE. Theorem 1 holds for any s_k with
     * alpha_k(n) = alpha_k(s_k) for all n &gt;= s_k, and a LARGER s_k is always admissible, so
     * s_k is detected here as the smallest index whose value the tail of mu(k,:) repeats to
     * within tol. A station whose rates never settle (an infinite server, mu(k,n) = n) gets
     * s_k = |N|, which is still exact: populations above |N| do not occur, so redefining alpha_k
     * there changes nothing. It is merely expensive, since the inner sum costs prod_k s_k terms,
     * capped by |v| &lt;= |N|. Think time is not admissible: a delay would have to enter g_sigma,
     * whose closed form covers queues only.</p>
     *
     * <p>NUMERICS. Both sums alternate in sign with terms far larger than the result, so they are
     * evaluated as signed log-sum-exps. phi_k is sign-definite when alpha_k increases, as a
     * multi-server station does, and changes sign where alpha_k decreases, so a decreasing rate
     * function costs digits in the inner sum too.</p>
     *
     * <p>SINGLE CLASS. At R=1 the divided difference is the identity, since h_theta(N) is
     * homogeneous of degree N in theta exactly as in the fixed-rate case, so the outer sum is
     * skipped and Theorem 1 is evaluated once at theta = L.</p>
     *
     * @param L       service demand matrix (M x R)
     * @param N       population vector (1 x R)
     * @param mu      load-dependent rate matrix (M x sum(N)), alpha_i(j) = mu(i,j-1); null for
     *                all ones
     * @param tol     relative tolerance declaring two scaled demands redundant, and the rate tail
     *                constant
     * @param method  "auto", "distinct" (force Eq. 15) or "repeated" (force Eq. 16)
     * @param maxloss cancellation budget in decimal digits; a finite value turns the warnings
     *                into a silent REFUSAL (lG = NaN) once the budget is exceeded, for callers
     *                that hold a fallback; Double.POSITIVE_INFINITY keeps the warnings
     * @return the constant, its logarithm, the expression used and the digits lost
     */
    public static Result pfqn_explicit_ld(Matrix L, Matrix N, Matrix mu, double tol,
                                          String method, double maxloss) {
        int R = N.length();
        if (Double.isNaN(tol)) {
            tol = Math.ulp(1.0); // machine precision, the tolerance is relative to max(sigma)
        }
        if (method == null || method.isEmpty()) {
            method = "auto";
        }
        if (!"auto".equals(method) && !"distinct".equals(method) && !"repeated".equals(method)) {
            throw new IllegalArgumentException(
                    "pfqn_explicit_ld: unrecognized method, use 'auto', 'distinct' (Eq. 15) or "
                            + "'repeated' (Eq. 16).");
        }
        double Ntd = 0.0;
        for (int r = 0; r < R; r++) {
            Ntd += N.get(r);
        }
        if (Ntd < 0) {
            return new Result(Double.NEGATIVE_INFINITY, 0.0, "distinct", 0.0);
        }
        if (Ntd == 0) {
            return new Result(0.0, 1.0, "distinct", 0.0);
        }
        if (L == null || L.isEmpty()) {
            return new Result(Double.NEGATIVE_INFINITY, 0.0, "distinct", 0.0);
        }
        if (L.getNumCols() != R) {
            throw new IllegalArgumentException(
                    "pfqn_explicit_ld: the demand matrix must have one column per class of N.");
        }
        int M = L.getNumRows();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) < 0) {
                    throw new IllegalArgumentException(
                            "pfqn_explicit_ld: the demand matrix must be nonnegative.");
                }
            }
        }
        int Nt = (int) FastMath.rint(Ntd);
        double[][] alpha = new double[M][Nt];
        if (mu == null || mu.isEmpty()) {
            for (int i = 0; i < M; i++) {
                for (int n = 0; n < Nt; n++) {
                    alpha[i][n] = 1.0;
                }
            }
        } else {
            if (mu.getNumRows() != M) {
                throw new IllegalArgumentException("pfqn_explicit_ld: the load-dependent rate "
                        + "matrix must have one row per station of L.");
            }
            if (mu.getNumCols() < Nt) {
                throw new IllegalArgumentException("pfqn_explicit_ld: the load-dependent rate "
                        + "matrix must have at least sum(N) columns.");
            }
            for (int i = 0; i < M; i++) {
                for (int n = 0; n < Nt; n++) {
                    alpha[i][n] = mu.get(i, n);
                    if (alpha[i][n] <= 0) {
                        throw new IllegalArgumentException("pfqn_explicit_ld: the load-dependent "
                                + "rates must be strictly positive.");
                    }
                }
            }
        }

        // ---- s_k: the smallest index whose value the tail of the rate row repeats ----
        // Any larger s_k also satisfies alpha_k(n)=alpha_k(s_k) for n>=s_k, so a missed tie only
        // adds terms; a false tie would be a wrong answer, hence the strict tol.
        int[] s = new int[M];
        double[] alphaS = new double[M];
        for (int i = 0; i < M; i++) {
            s[i] = Nt;
            double tail = alpha[i][Nt - 1];
            for (int n = Nt; n > 1; n--) {
                if (Math.abs(alpha[i][n - 2] - tail) <= tol * Math.max(Math.abs(tail), 1.0)) {
                    s[i] = n - 1;
                } else {
                    break;
                }
            }
            alphaS[i] = alpha[i][s[i] - 1];
        }

        // ---- per-station phi tables, in the log domain, indexed by v_k = 0..s_k-1 ----
        double[][] lcum = new double[M][]; // sum_{t=1..v} log alpha_k(t)
        double[][] lbr = new double[M][];  // log|1 - alpha_k(v)/alpha_k(s_k)|, alpha_k(0) := 0
        double[][] sbr = new double[M][];
        int[] vcap = new int[M];
        for (int i = 0; i < M; i++) {
            lcum[i] = new double[s[i]];
            lbr[i] = new double[s[i]];
            sbr[i] = new double[s[i]];
            for (int v = 0; v < s[i]; v++) {
                lcum[i][v] = (v == 0) ? 0.0 : lcum[i][v - 1] + FastMath.log(alpha[i][v - 1]);
                double br = 1.0 - ((v == 0) ? 0.0 : alpha[i][v - 1]) / alphaS[i];
                lbr[i][v] = (br == 0) ? Double.NEGATIVE_INFINITY : FastMath.log(Math.abs(br));
                sbr[i][v] = Math.signum(br);
            }
            // g_sigma vanishes below zero population, Eq. (8) caps |v| <= |N|
            vcap[i] = Math.min(s[i] - 1, Nt);
        }

        // ---- redundancy scan: are the SCALED induced demands pairwise distinct? ----
        // The scan MUST form sigma exactly as hlld does, (L*t)/alphaS and not
        // (L/alphaS)*t: the two orderings differ in the last ulp, so an exact tie can
        // clear an eps-relative gap under one and not the other, and Eq. (15) would then
        // divide by that ulp.
        int[] Nv = new int[R];
        for (int r = 0; r < R; r++) {
            Nv[r] = (int) FastMath.rint(N.get(r));
        }
        boolean isRedundant;
        if (R == 1) {
            // the scaled demands at t are t*sigma, so both the tie structure and the
            // relative tolerance are those of sigma itself, at every t at once
            double[] th = new double[M];
            for (int i = 0; i < M; i++) {
                th[i] = L.get(i, 0) / alphaS[i];
            }
            isRedundant = Pfqn_explicit.redundantAt(th, tol);
        } else {
            isRedundant = false;
            int[] t = new int[R];
            while (true) {
                int ts = 0;
                for (int r = 0; r < R; r++) {
                    ts += t[r];
                }
                if (ts > 0 && Pfqn_explicit.redundantAt(scaled(L, t, M, R, alphaS), tol)) {
                    isRedundant = true;
                    break;
                }
                if (!Pfqn_explicit.nextLattice(t, Nv)) {
                    break;
                }
            }
        }
        String expr = method;
        if ("auto".equals(expr)) {
            expr = isRedundant ? "repeated" : "distinct";
        } else if ("distinct".equals(expr) && isRedundant) {
            throw new IllegalArgumentException("pfqn_explicit_ld: Eq. (15) requires pairwise "
                    + "distinct scaled demands, but two of them agree to within tol. Use 'auto' "
                    + "or 'repeated'.");
        }

        Pfqn_explicit.SignedLse total;
        if (R == 1) {
            // ---- single class: the divided difference is the identity ----
            double[] th = new double[M];
            for (int i = 0; i < M; i++) {
                th[i] = L.get(i, 0);
            }
            total = hlld(th, M, Nt, alphaS, vcap, lcum, lbr, sbr, expr, tol);
        } else {
            // ---- outer divided-difference sum over 0 <= t <= N ----
            List<Double> lterm = new ArrayList<Double>();
            List<Double> sterm = new ArrayList<Double>();
            double innerLoss = 0.0;
            int[] t = new int[R];
            while (true) {
                int ts = 0;
                for (int r = 0; r < R; r++) {
                    ts += t[r];
                }
                if (ts > 0) {
                    double[] th = Pfqn_explicit.induced(L, t, M, R);
                    double thmax = 0.0;
                    for (int i = 0; i < M; i++) {
                        thmax = Math.max(thmax, th[i]);
                    }
                    if (thmax > 0) {
                        Pfqn_explicit.SignedLse h =
                                hlld(th, M, Nt, alphaS, vcap, lcum, lbr, sbr, expr, tol);
                        innerLoss = Math.max(innerLoss, h.lossDigits);
                        if (h.sgn != 0) {
                            double l = h.lS;
                            for (int r = 0; r < R; r++) {
                                l -= Pfqn_explicit.factln(t[r]);
                                l -= Pfqn_explicit.factln(Nv[r] - t[r]);
                            }
                            lterm.add(l);
                            sterm.add(h.sgn * (((Nt - ts) % 2 == 0) ? 1.0 : -1.0));
                        }
                    }
                }
                if (!Pfqn_explicit.nextLattice(t, Nv)) {
                    break;
                }
            }
            total = Pfqn_explicit.signedLogSumExp(Pfqn_explicit.toArray(lterm),
                    Pfqn_explicit.toArray(sterm));
            total = new Pfqn_explicit.SignedLse(total.lS, total.sgn,
                    Math.max(total.lossDigits, innerLoss));
        }

        // A caller that named a cancellation budget has a fallback and wants a verdict, not a
        // warning: refuse quietly. lossDigits is infinite when the sum vanished identically,
        // which is a total loss rather than a legitimate G = 0.
        if (!Double.isInfinite(maxloss) && (total.sgn < 0 || total.lossDigits > maxloss)) {
            return new Result(Double.NaN, Double.NaN, expr, total.lossDigits);
        }
        if (total.sgn == 0) {
            return new Result(Double.NEGATIVE_INFINITY, 0.0, expr, total.lossDigits);
        }
        if (total.sgn < 0) {
            InputOutput.line_warning("pfqn_explicit_ld",
                    "The explicit expression returned a negative value, double precision is "
                            + "exhausted by cancellation (%.1f digits lost). Multiprecision "
                            + "arithmetic is required.\n", total.lossDigits);
            return new Result(Double.NaN, Double.NaN, expr, total.lossDigits);
        }
        if (total.lossDigits > 15) {
            InputOutput.line_warning("pfqn_explicit_ld",
                    "Cancellation has consumed about %.1f decimal digits, more than double "
                            + "precision carries. The result is unreliable, multiprecision "
                            + "arithmetic is required.\n", total.lossDigits);
        }
        return new Result(total.lS, FastMath.exp(total.lS), expr, total.lossDigits);
    }

    /** Overload with the documented defaults: machine precision, "auto", no budget. */
    public static Result pfqn_explicit_ld(Matrix L, Matrix N, Matrix mu) {
        return pfqn_explicit_ld(L, N, mu, Double.NaN, "auto", Double.POSITIVE_INFINITY);
    }

    /** Scaled induced demands sigma_k(t) = (sum_r t_r L(k,r)) / alpha_k(s_k). */
    private static double[] scaled(Matrix L, int[] t, int M, int R, double[] alphaS) {
        double[] th = Pfqn_explicit.induced(L, t, M, R);
        for (int i = 0; i < M; i++) {
            th[i] = th[i] / alphaS[i];
        }
        return th;
    }

    /**
     * Theorem 1 of Casale-Harrison-Ong (2021), Eq. (8): the single-class limited load-dependent
     * constant at induced demands th and total population Nt, as the finite sum over
     * 0 &lt;= v &lt; s of the fixed-rate constant at the scaled demands th/alpha(s), one
     * population level lower for every job held back by v.
     */
    private static Pfqn_explicit.SignedLse hlld(double[] th, int M, int Nt, double[] alphaS,
                                                int[] vcap, double[][] lcum, double[][] lbr,
                                                double[][] sbr, String expr, double tol) {
        double[] sigma = new double[M];
        double[] lth = new double[M];
        for (int i = 0; i < M; i++) {
            sigma[i] = th[i] / alphaS[i];
            lth[i] = (th[i] > 0) ? FastMath.log(th[i]) : Double.NEGATIVE_INFINITY;
        }
        List<Double> lterm = new ArrayList<Double>();
        List<Double> sterm = new ArrayList<Double>();
        double lossDigits = 0.0;
        int[] v = new int[M];
        while (true) {
            int nv = 0;
            for (int i = 0; i < M; i++) {
                nv += v[i];
            }
            if (nv <= Nt) {
                double lval = 0.0;
                double sval = 1.0;
                for (int i = 0; i < M; i++) {
                    int vi = v[i];
                    if (vi > 0) {
                        if (Double.isInfinite(lth[i])) {
                            // theta_k = 0 kills every v_k>0, and 0^0=1 keeps v_k=0
                            sval = 0.0;
                            break;
                        }
                        // kept inside the guard because 0*(-Inf) is NaN, not 0
                        lval += vi * lth[i];
                    }
                    lval += -lcum[i][vi] + lbr[i][vi];
                    sval *= sbr[i][vi];
                }
                if (sval != 0 && !Double.isInfinite(lval) && !Double.isNaN(lval)) {
                    double lg;
                    int sg;
                    double dl;
                    if (Nt - nv == 0) {
                        // g_sigma(0) = 1 by definition. Reading it off the partial fraction
                        // instead would spend digits on an alternating sum whose value is
                        // known exactly.
                        lg = 0.0;
                        sg = 1;
                        dl = 0.0;
                    } else {
                        Pfqn_explicit.SignedLse g = "distinct".equals(expr)
                                ? Pfqn_explicit.gdistinct(sigma, Nt - nv, M)
                                : Pfqn_explicit.grepeated(sigma, Nt - nv, M, tol);
                        lg = g.lS;
                        sg = g.sgn;
                        dl = g.lossDigits;
                    }
                    lossDigits = Math.max(lossDigits, dl);
                    if (sg != 0) {
                        lterm.add(lval + lg);
                        sterm.add(sval * sg);
                    }
                }
            }
            if (!Pfqn_explicit.nextLattice(v, vcap)) {
                break;
            }
        }
        Pfqn_explicit.SignedLse h = Pfqn_explicit.signedLogSumExp(Pfqn_explicit.toArray(lterm),
                Pfqn_explicit.toArray(sterm));
        return new Pfqn_explicit.SignedLse(h.lS, h.sgn, Math.max(lossDigits, h.lossDigits));
    }
}
