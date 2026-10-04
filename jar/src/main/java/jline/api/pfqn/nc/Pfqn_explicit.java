/**
 * @file Explicit closed-form normalizing constant of a multiclass closed network
 *
 * Port at parity of matlab/src/api/pfqn/pfqn_explicit.m, i.e. Eqs. (15) and (16) of
 * G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic
 * Methods", ACM SIGMETRICS 2017.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Pfqn_explicit {
    private Pfqn_explicit() {}

    /** Result of the closed form, mirroring [lG, G, method, lossDigits]. */
    public static final class Result {
        /** Logarithm of the normalizing constant. */
        public final double lG;
        /** The normalizing constant. */
        public final double G;
        /** Expression used, "distinct" (Eq. 15) or "repeated" (Eq. 16). */
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
     * Signed log-sum-exp of S = sum_i s_i exp(l_i): log|S|, sign(S), digits lost.
     *
     * <p>Public because Pfqn_explicit_ld carries the same closed form as the inner kernel of
     * its own two sums and needs the sign alongside the magnitude.</p>
     */
    public static final class SignedLse {
        /** Logarithm of |S|. */
        public final double lS;
        /** Sign of S, zero when every term dropped out. */
        public final int sgn;
        /** Decimal digits lost to cancellation. */
        public final double lossDigits;

        public SignedLse(double lS, int sgn, double lossDigits) {
            this.lS = lS;
            this.sgn = sgn;
            this.lossDigits = lossDigits;
        }
    }

    /**
     * Explicit closed-form normalizing constant of a multiclass closed network.
     *
     * <p>Evaluates the two explicit expressions of Casale (SIGMETRICS 2017), Eqs. (15) and
     * (16). Both instantiate the divided-difference form of Corollary 3.2,
     *
     * <pre>
     *   G(N) = sum_{0&lt;=t&lt;=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r C(N_r,t_r) g_t(|N|)
     * </pre>
     *
     * by substituting a closed form for the single-class constant g_t(|N|) at the induced
     * demands theta_k(t) = sum_r t_r L(k,r). Eq. (15) is Gordon's partial fraction and needs
     * the induced demands PAIRWISE DISTINCT; Eq. (16) is the general partial-fraction
     * expansion over the distinct values and their multiplicities, and reduces term by term
     * to Eq. (15) when every multiplicity is one. The choice is automatic: Eq. (16) is used
     * as soon as two induced demands are closer than tol relative to the largest one.</p>
     *
     * <p>SINGLE CLASS. At R=1 the multiclass constant IS the single-class constant at demands
     * L, so the outer sum is skipped: g_t(N) = t^N g_1(N) and
     * sum_t (-1)^(N-t) t^N/(t!(N-t)!) = S(N,N) = 1. Running the difference anyway would add N
     * alternating terms, and their cancellation, to a closed form that carries none of them.
     * What is left is O(K^2) work at any population.</p>
     *
     * <p>Only single-server load-independent queues are admissible: infinite servers need the
     * integral form of Corollary 3.4 and load-dependent rates need Pfqn_explicit_ld, which keeps
     * this closed form as its inner kernel.</p>
     *
     * <p>NUMERICS. Both expressions alternate in sign with terms far larger than the result,
     * so they are evaluated as signed log-sum-exps: this removes the floating-point RANGE
     * problem but not the cancellation, which is what makes multiprecision arithmetic
     * necessary on all but small models.</p>
     *
     * @param L       service demand matrix (K x R) of single-server load-independent queues
     * @param N       population vector (1 x R)
     * @param tol     relative tolerance declaring two induced demands redundant
     * @param method  "auto", "distinct" (force Eq. 15) or "repeated" (force Eq. 16)
     * @param maxloss cancellation budget in decimal digits; a finite value turns the warnings
     *                into a silent REFUSAL (lG = NaN) once the budget is exceeded, for callers
     *                that hold a fallback; Double.POSITIVE_INFINITY keeps the warnings
     * @return the constant, its logarithm, the expression used and the digits lost
     */
    public static Result pfqn_explicit(Matrix L, Matrix N, double tol, String method,
                                       double maxloss) {
        int R = N.length();
        if (Double.isNaN(tol)) {
            tol = Math.ulp(1.0); // machine precision, the tolerance is relative to max(theta)
        }
        if (method == null || method.isEmpty()) {
            method = "auto";
        }
        if (!"auto".equals(method) && !"distinct".equals(method) && !"repeated".equals(method)) {
            throw new IllegalArgumentException(
                    "pfqn_explicit: unrecognized method, use 'auto', 'distinct' (Eq. 15) or "
                            + "'repeated' (Eq. 16).");
        }
        double Nt = 0.0;
        for (int r = 0; r < R; r++) {
            Nt += N.get(r);
        }
        if (Nt < 0) {
            return new Result(Double.NEGATIVE_INFINITY, 0.0, "distinct", 0.0);
        }
        if (Nt == 0) {
            return new Result(0.0, 1.0, "distinct", 0.0);
        }
        if (L == null || L.isEmpty()) {
            return new Result(Double.NEGATIVE_INFINITY, 0.0, "distinct", 0.0);
        }
        if (L.getNumCols() != R) {
            throw new IllegalArgumentException(
                    "pfqn_explicit: the demand matrix must have one column per class of N.");
        }
        int K = L.getNumRows();
        for (int i = 0; i < K; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) < 0) {
                    throw new IllegalArgumentException(
                            "pfqn_explicit: the demand matrix must be nonnegative.");
                }
            }
        }

        // ---- redundancy scan: are the induced demands pairwise distinct at every t? ----
        boolean isRedundant;
        if (R == 1) {
            // the induced demands at t are t*L, so both the tie structure and the relative
            // tolerance are those of L itself, at every t at once
            double[] th = new double[K];
            for (int i = 0; i < K; i++) {
                th[i] = L.get(i, 0);
            }
            isRedundant = redundantAt(th, tol);
        } else {
            isRedundant = false;
            int[] Nv = new int[R];
            for (int r = 0; r < R; r++) {
                Nv[r] = (int) FastMath.rint(N.get(r));
            }
            int[] t = new int[R];
            while (true) {
                int ts = 0;
                for (int r = 0; r < R; r++) {
                    ts += t[r];
                }
                if (ts > 0 && redundantAt(induced(L, t, K, R), tol)) {
                    isRedundant = true;
                    break;
                }
                if (!nextLattice(t, Nv)) {
                    break;
                }
            }
        }
        String expr = method;
        if ("auto".equals(expr)) {
            expr = isRedundant ? "repeated" : "distinct";
        } else if ("distinct".equals(expr) && isRedundant) {
            throw new IllegalArgumentException(
                    "pfqn_explicit: Eq. (15) requires pairwise distinct induced demands, but two "
                            + "of them agree to within tol. Use 'auto' or 'repeated'.");
        }

        SignedLse total;
        if (R == 1) {
            // ---- single class: the divided difference is the identity, evaluate g ----
            double[] th = new double[K];
            for (int i = 0; i < K; i++) {
                th[i] = L.get(i, 0);
            }
            total = "distinct".equals(expr) ? gdistinct(th, Nt, K) : grepeated(th, Nt, K, tol);
        } else {
            // ---- outer divided-difference sum over 0 <= t <= N ----
            List<Double> lterm = new ArrayList<Double>();
            List<Double> sterm = new ArrayList<Double>();
            double innerLoss = 0.0;
            int[] Nv = new int[R];
            for (int r = 0; r < R; r++) {
                Nv[r] = (int) FastMath.rint(N.get(r));
            }
            int[] t = new int[R];
            while (true) {
                int ts = 0;
                for (int r = 0; r < R; r++) {
                    ts += t[r];
                }
                if (ts > 0) {
                    double[] th = induced(L, t, K, R);
                    double thmax = 0.0;
                    for (int i = 0; i < K; i++) {
                        thmax = Math.max(thmax, th[i]);
                    }
                    if (thmax > 0) {
                        SignedLse g = "distinct".equals(expr) ? gdistinct(th, Nt, K)
                                : grepeated(th, Nt, K, tol);
                        innerLoss = Math.max(innerLoss, g.lossDigits);
                        if (g.sgn != 0) {
                            double l = g.lS;
                            for (int r = 0; r < R; r++) {
                                l -= factln(t[r]);
                                l -= factln(Nv[r] - t[r]);
                            }
                            lterm.add(l);
                            sterm.add(g.sgn * ((((int) Nt - ts) % 2 == 0) ? 1.0 : -1.0));
                        }
                    }
                }
                if (!nextLattice(t, Nv)) {
                    break;
                }
            }
            total = signedLogSumExp(toArray(lterm), toArray(sterm));
            total = new SignedLse(total.lS, total.sgn, Math.max(total.lossDigits, innerLoss));
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
            InputOutput.line_warning("pfqn_explicit",
                    "The explicit expression returned a negative value, double precision is "
                            + "exhausted by cancellation (%.1f digits lost). Multiprecision "
                            + "arithmetic is required.\n", total.lossDigits);
            return new Result(Double.NaN, Double.NaN, expr, total.lossDigits);
        }
        if (total.lossDigits > 15) {
            InputOutput.line_warning("pfqn_explicit",
                    "Cancellation has consumed about %.1f decimal digits, more than double "
                            + "precision carries. The result is unreliable, multiprecision "
                            + "arithmetic is required.\n", total.lossDigits);
        }
        return new Result(total.lS, FastMath.exp(total.lS), expr, total.lossDigits);
    }

    /** Overload with the documented defaults: machine precision, 'auto', no budget. */
    public static Result pfqn_explicit(Matrix L, Matrix N) {
        return pfqn_explicit(L, N, Double.NaN, "auto", Double.POSITIVE_INFINITY);
    }

    /** Induced demands theta_k(t) = sum_r t_r L(k,r). Shared with Pfqn_explicit_ld. */
    public static double[] induced(Matrix L, int[] t, int K, int R) {
        double[] th = new double[K];
        for (int i = 0; i < K; i++) {
            double v = 0.0;
            for (int r = 0; r < R; r++) {
                v += L.get(i, r) * t[r];
            }
            th[i] = v;
        }
        return th;
    }

    /**
     * Whether two of the induced demands agree to within tol, relatively to the largest.
     * A scale of zero leaves every induced demand at zero, so g_t(|N|) = 0 at sum(N) &gt; 0
     * and the term takes no part in the sum.
     */
    public static boolean redundantAt(double[] th, double tol) {
        double[] s = th.clone();
        Arrays.sort(s);
        double scale = s[s.length - 1];
        if (!(scale > 0)) {
            return false;
        }
        for (int i = 1; i < s.length; i++) {
            if (s[i] - s[i - 1] <= tol * scale) {
                return true;
            }
        }
        return false;
    }

    /** Next vector of the lattice 0 &lt;= t &lt;= N; false once it is exhausted. Shared with Pfqn_explicit_ld. */
    public static boolean nextLattice(int[] t, int[] N) {
        int r = t.length;
        while (r > 0 && t[r - 1] == N[r - 1]) {
            t[--r] = 0;
        }
        if (r == 0) {
            return false;
        }
        t[r - 1]++;
        return true;
    }

    /**
     * Eq. (14): single-class constant at pairwise distinct demands th, population Nt over K
     * queues. A zero demand contributes nothing, which also realizes the 0/0 = 0 convention
     * of Eq. (15) when the zero is repeated.
     */
    public static SignedLse gdistinct(double[] th, double Nt, int K) {
        double[] lin = new double[K];
        double[] sgv = new double[K];
        Arrays.fill(lin, Double.NEGATIVE_INFINITY);
        for (int k = 0; k < K; k++) {
            if (th[k] <= 0) {
                continue;
            }
            double acc = (Nt + K - 1) * FastMath.log(th[k]);
            double sign = 1.0;
            for (int i = 0; i < K; i++) {
                if (i == k) {
                    continue;
                }
                double d = th[k] - th[i];
                acc -= FastMath.log(Math.abs(d));
                sign *= Math.signum(d);
            }
            lin[k] = acc;
            sgv[k] = sign;
        }
        return signedLogSumExp(lin, sgv);
    }

    /**
     * Eq. (16): single-class constant at demands th of arbitrary multiplicity, population Nt
     * over K queues. Demands within tol of each other, relatively to the largest one, are
     * merged into one distinct value carrying their count.
     */
    public static SignedLse grepeated(double[] th, double Nt, int K, double tol) {
        double[] ths = th.clone();
        Arrays.sort(ths);
        double scale = ths[ths.length - 1];
        if (scale <= 0) {
            scale = 1.0;
        }
        List<Double> thdList = new ArrayList<Double>();
        List<Integer> mList = new ArrayList<Integer>();
        int i = 0;
        while (i < ths.length) {
            int j = i + 1;
            while (j < ths.length && !(ths[j] - ths[j - 1] > tol * scale)) {
                j++;
            }
            double sum = 0.0;
            for (int q = i; q < j; q++) {
                sum += ths[q];
            }
            thdList.add(sum / (j - i)); // the centroid represents a cluster of near-ties
            mList.add(j - i);
            i = j;
        }
        int Kp = thdList.size();
        double[] thd = new double[Kp];
        int[] m = new int[Kp];
        for (int j = 0; j < Kp; j++) {
            thd[j] = thdList.get(j);
            m[j] = mList.get(j);
        }
        List<Double> lin = new ArrayList<Double>();
        List<Double> sgv = new ArrayList<Double>();
        for (int j = 0; j < Kp; j++) {
            if (thd[j] <= 0) {
                // the exponent Nt+K-m_j is at least Nt>=1, so a zero cluster contributes nothing
                continue;
            }
            double louter = (Nt + K - m[j]) * FastMath.log(thd[j]);
            double souter = ((m[j] - 1) % 2 == 0) ? 1.0 : -1.0;
            List<int[]> rs = multichoose(Kp, m[j] - 1); // every K'-vector r>=0 with sum(r)=m_j-1
            for (int a = 0; a < rs.size(); a++) {
                int[] r = rs.get(a);
                double lval = louter + nchoosekln(Nt + r[j], r[j]);
                double sval = souter * ((r[j] % 2 == 0) ? 1.0 : -1.0);
                boolean vanished = false;
                for (int k = 0; k < Kp; k++) {
                    if (k == j) {
                        continue;
                    }
                    lval += nchoosekln(m[k] + r[k] - 1, r[k]);
                    if (r[k] > 0) {
                        if (thd[k] <= 0) {
                            // theta_k^r_k vanishes, 0^0 = 1 is the r_k = 0 case
                            vanished = true;
                            break;
                        }
                        lval += r[k] * FastMath.log(thd[k]);
                    }
                    double dd = thd[j] - thd[k];
                    lval -= (m[k] + r[k]) * FastMath.log(Math.abs(dd));
                    if (dd < 0 && (m[k] + r[k]) % 2 != 0) {
                        sval = -sval;
                    }
                }
                if (vanished) {
                    lin.add(Double.NEGATIVE_INFINITY);
                    sgv.add(0.0);
                } else {
                    lin.add(lval);
                    sgv.add(sval);
                }
            }
        }
        return signedLogSumExp(toArray(lin), toArray(sgv));
    }

    /** Every Kp-vector r &gt;= 0 with sum(r) = k, i.e. MATLAB multichoose(Kp,k). */
    private static List<int[]> multichoose(int Kp, int k) {
        List<int[]> out = new ArrayList<int[]>();
        if (Kp <= 0 || k < 0) {
            return out;
        }
        multichooseRec(Kp, k, new int[Kp], 0, out);
        return out;
    }

    private static void multichooseRec(int Kp, int k, int[] current, int idx, List<int[]> out) {
        if (idx == Kp - 1) {
            current[idx] = k;
            out.add(current.clone());
            return;
        }
        for (int i = 0; i <= k; i++) {
            current[idx] = i;
            multichooseRec(Kp, k - i, current, idx + 1, out);
        }
    }

    /**
     * Signed log-sum-exp of S = sum_i s_i exp(l_i). Returns log|S|, sign(S) and the decimal
     * digits lost to cancellation.
     */
    public static SignedLse signedLogSumExp(double[] lterm, double[] sterm) {
        int n = 0;
        double a = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < lterm.length; i++) {
            if (!Double.isInfinite(lterm[i]) && !Double.isNaN(lterm[i]) && sterm[i] != 0) {
                n++;
                if (lterm[i] > a) {
                    a = lterm[i];
                }
            }
        }
        if (n == 0) {
            return new SignedLse(Double.NEGATIVE_INFINITY, 0, 0.0);
        }
        double s = 0.0;
        for (int i = 0; i < lterm.length; i++) {
            if (!Double.isInfinite(lterm[i]) && !Double.isNaN(lterm[i]) && sterm[i] != 0) {
                s += sterm[i] * FastMath.exp(lterm[i] - a);
            }
        }
        if (s == 0) {
            return new SignedLse(Double.NEGATIVE_INFINITY, 0, Double.POSITIVE_INFINITY);
        }
        // max(exp(l_i - a)) is 1, so -log10|s| is the shortfall of the sum against its largest
        // term. Every one of the n terms carries a rounding error of order eps*max_term, so the
        // digits actually lost are that shortfall PLUS log10(n); dropping the count understates
        // the loss and lets a wrong answer past the guard.
        double lossDigits = Math.max(0.0, FastMath.log10(n / Math.abs(s)));
        return new SignedLse(a + FastMath.log(Math.abs(s)), s > 0 ? 1 : -1, lossDigits);
    }

    public static double factln(double n) {
        return Gamma.logGamma(1.0 + n);
    }

    /** Logarithm of the binomial coefficient C(n,m); MATLAB nchoosekln.m. */
    private static double nchoosekln(double n, double m) {
        return Gamma.logGamma(1.0 + n) - Gamma.logGamma(1.0 + n - m) - Gamma.logGamma(1.0 + m);
    }

    public static double[] toArray(List<Double> v) {
        double[] out = new double[v.size()];
        for (int i = 0; i < v.size(); i++) {
            out[i] = v.get(i);
        }
        return out;
    }
}
