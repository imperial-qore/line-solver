/**
 * @file Multiclass Recursion by Generating Functions (RGF), with think times
 *
 * Exact normalizing constant of a closed MULTICLASS product-form network by
 * eliminating one class at a time by residues, finishing in the single-class
 * convolution of {@link Pfqn_rgf}. Ported at parity from MATLAB pfqn_rgfmc.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;

/**
 * Harrison and Coury, "On the asymptotic behaviour of closed multiclass queueing
 * networks", Perf. Eval. 47:131-138, 2002, Thm 1, as algorithmised by Harrison
 * and Lee, "A new recursive algorithm for computing generating functions in
 * closed multi-class queueing networks", IEEE MASCOTS 2004, eqs. (4)-(5).
 *
 * <p>THINK TIMES ARE NOT IN EITHER PAPER. Both write the generating function as
 * the RATIONAL prod_i (1 - rho_i z)^-m_i, every node a load-independent single
 * server. An infinite server multiplies it by the ENTIRE exp(sum_r Z_r z_r),
 * which breaks the step Thm 1 rests on: G_n(z') = -sum_i r_i holds only because
 * the residues of n(z)/d(z) sum to zero when deg d &gt;= deg n + 2 (Bertozzi and
 * McKenna, SIAM Review 35(2):239-268, 1993, fact (IV), p. 246), and an
 * exponential numerator does not decay at infinity. The delay is therefore
 * carried by their own repair, eqs. (3.19)-(3.21): only the first k_r+1 Taylor
 * coefficients of exp(Z_r z_r) can reach the coefficient of z_r^k_r, so
 * replacing the exponential by that polynomial is EXACT and leaves a rational
 * integrand. The price is that the eliminated class's population re-enters the
 * term count, which is the population-insensitivity Harrison-Lee Sec. 4
 * advertises; the class kept for the base case pays nothing.</p>
 *
 * <p>The elimination is exact in exact arithmetic but is an ALTERNATING sum over
 * residues, so near-coincident loads over an eliminated class destroy
 * significance. The worst cancellation ratio is tracked and the routine REFUSES
 * past maxcancel rather than returning a confidently wrong lG.</p>
 */
public final class Pfqn_rgfmc {
    private Pfqn_rgfmc() {}

    /** Default relative tolerance for calling two affine forms proportional. */
    public static final double DEF_TOL = 1e-12;
    /** Default cap on residue terms carried between eliminations. */
    public static final int DEF_MAXTERMS = 1000000;
    /** Default nats of cancellation tolerated before refusing. */
    public static final double DEF_MAXCANCEL = 15.0;

    /** Normalizing constant and its logarithm. */
    public static final class Result {
        /** Normalizing constant G(N). */
        public final double G;
        /** Logarithm of G(N). */
        public final double lG;

        public Result(double G, double lG) {
            this.G = G;
            this.lG = lG;
        }
    }

    /** coef * prod_t ( F[t][0] + sum_{p>=1} F[t][p] z_p ) ^ (-m[t]). */
    private static final class Term {
        double lc;
        double sc;
        double[][] F;
        int[] m;

        Term(double lc, double sc, double[][] F, int[] m) {
            this.lc = lc;
            this.sc = sc;
            this.F = F;
            this.m = m;
        }
    }

    /**
     * Exact multiclass normalizing constant by recursion on generating functions.
     *
     * @param L service demand matrix (M x R)
     * @param N population vector (R), nonnegative integers
     * @param Z think time vector (R)
     * @return the normalizing constant and its logarithm
     */
    public static Result pfqn_rgfmc(Matrix L, double[] N, double[] Z) {
        return pfqn_rgfmc(L, N, Z, DEF_TOL, DEF_MAXTERMS, DEF_MAXCANCEL);
    }

    /**
     * @param L service demand matrix (M x R)
     * @param N population vector (R), nonnegative integers
     * @param Z think time vector (R)
     * @param tol relative tolerance for calling two affine forms proportional
     * @param maxterms cap on residue terms carried between eliminations
     * @param maxcancel nats of cancellation tolerated before refusing
     * @return the normalizing constant and its logarithm
     */
    public static Result pfqn_rgfmc(Matrix L, double[] N, double[] Z,
                                    double tol, int maxterms, double maxcancel) {
        int R = L.getNumCols();
        if (N.length != R || Z.length != R) {
            throw new IllegalArgumentException(
                    "pfqn_rgfmc requires numel(N) and numel(Z) to match the number of columns of L.");
        }
        for (int r = 0; r < R; r++) {
            if (Z[r] < 0 || N[r] < 0) {
                throw new IllegalArgumentException("pfqn_rgfmc requires nonnegative L, N and Z.");
            }
            if (N[r] != FastMath.rint(N[r])) {
                throw new IllegalArgumentException("pfqn_rgfmc requires integer populations.");
            }
        }
        // Drop the classes with no jobs, then the stations with no demand at all.
        List<Integer> kc = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (N[r] > 0) {
                kc.add(Integer.valueOf(r));
            }
        }
        int Rk = kc.size();
        if (Rk == 0) {
            return new Result(1.0, 0.0);
        }
        List<Integer> kr = new ArrayList<Integer>();
        for (int i = 0; i < L.getNumRows(); i++) {
            boolean any = false;
            for (int c = 0; c < Rk; c++) {
                double d = L.get(i, kc.get(c).intValue());
                if (d < 0) {
                    throw new IllegalArgumentException("pfqn_rgfmc requires nonnegative L, N and Z.");
                }
                if (d > 0) {
                    any = true;
                }
            }
            if (any) {
                kr.add(Integer.valueOf(i));
            }
        }
        int M = kr.size();
        double[] Nk = new double[Rk];
        double[] Zk = new double[Rk];
        for (int c = 0; c < Rk; c++) {
            Nk[c] = N[kc.get(c).intValue()];
            Zk[c] = Z[kc.get(c).intValue()];
        }
        if (M == 0) {
            double lG = 0.0;
            for (int c = 0; c < Rk; c++) {
                lG += Nk[c] * FastMath.log(Zk[c]) - Gamma.logGamma(Nk[c] + 1.0);
            }
            return new Result(FastMath.exp(lG), lG);
        }
        if (Rk == 1) {
            Matrix col = new Matrix(M, 1);
            for (int i = 0; i < M; i++) {
                col.set(i, 0, L.get(kr.get(i).intValue(), kc.get(0).intValue()));
            }
            Pfqn_rgf.Result r1 = Pfqn_rgf.pfqn_rgf(col, Nk[0], Zk[0]);
            return new Result(r1.G, r1.lG);
        }
        // Base class = smallest population: that population is the degree the
        // sign-indefinite base series is carried to, so it drives cancellation.
        Integer[] ord = new Integer[Rk];
        for (int c = 0; c < Rk; c++) {
            ord[c] = Integer.valueOf(c);
        }
        final double[] Nsort = Nk;
        Arrays.sort(ord, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(Nsort[a.intValue()], Nsort[b.intValue()]);
            }
        });
        double[][] Ls = new double[M][Rk];
        double[] Ns = new double[Rk];
        double[] Zs = new double[Rk];
        for (int c = 0; c < Rk; c++) {
            int src = ord[c].intValue();
            Ns[c] = Nk[src];
            Zs[c] = Zk[src];
            for (int i = 0; i < M; i++) {
                Ls[i][c] = L.get(kr.get(i).intValue(), kc.get(src).intValue());
            }
        }
        double lGscale = 0.0;
        for (int c = 0; c < Rk; c++) {
            double cs = Zs[c];
            for (int i = 0; i < M; i++) {
                if (Ls[i][c] > cs) {
                    cs = Ls[i][c];
                }
            }
            if (cs <= 0) {
                cs = 1.0;
            }
            for (int i = 0; i < M; i++) {
                Ls[i][c] /= cs;
            }
            Zs[c] /= cs;
            lGscale += Ns[c] * FastMath.log(cs);
        }
        double[][] F0 = new double[M][Rk + 1];
        int[] m0 = new int[M];
        for (int i = 0; i < M; i++) {
            F0[i][0] = 1.0;
            for (int c = 0; c < Rk; c++) {
                F0[i][c + 1] = -Ls[i][c];
            }
            m0[i] = 1;
        }
        List<Term> terms = new ArrayList<Term>();
        terms.add(new Term(0.0, 1.0, F0, m0));
        double[] cond = new double[1];
        // F column p+1 carries class p, so eliminating class p means column p+1.
        for (int col = Rk; col >= 2; col--) {
            terms = step(terms, col, (int) FastMath.rint(Ns[col - 1]), Zs[col - 1],
                    tol, maxterms, cond);
            if (terms.isEmpty()) {
                return new Result(0.0, Double.NEGATIVE_INFINITY);
            }
        }
        double[] br = base(terms, (int) FastMath.rint(Ns[0]), Zs[0], tol, cond);
        if (br[1] == 0.0) {
            return new Result(0.0, Double.NEGATIVE_INFINITY);
        }
        if (br[1] < 0 || cond[0] > maxcancel) {
            throw new IllegalStateException(String.format(
                    "pfqn_rgfmc: the residue sum cancelled %.1f nats, past the %.1f allowed, so lG "
                    + "carries no significant digits. The eliminated classes have near-coincident "
                    + "loads over the stations. Use method 'ca' for the exact convolution.",
                    Double.valueOf(cond[0]), Double.valueOf(maxcancel)));
        }
        double lG = br[0] + lGscale;
        return new Result(FastMath.exp(lG), lG);
    }

    /** Signed log-domain sum; returns {log|sum|, sign} and tracks cancellation. */
    private static double[] slogsum(double[] lv, double[] sv, int n, double[] cond) {
        double mx = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            if (sv[i] != 0.0 && !Double.isInfinite(lv[i]) && lv[i] > mx) {
                mx = lv[i];
            }
        }
        if (Double.isInfinite(mx)) {
            return new double[]{Double.NEGATIVE_INFINITY, 0.0};
        }
        double tot = 0.0;
        double abs = 0.0;
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            if (sv[i] == 0.0 || Double.isInfinite(lv[i])) {
                continue;
            }
            double e = FastMath.exp(lv[i] - mx);
            tot += sv[i] * e;
            abs += e;
            cnt++;
        }
        if (tot == 0.0) {
            return new double[]{Double.NEGATIVE_INFINITY, 0.0};
        }
        double lsum = mx + FastMath.log(FastMath.abs(tot));
        if (cnt > 1) {
            double c = mx + FastMath.log(abs) - lsum;
            if (c > cond[0]) {
                cond[0] = c;
            }
        }
        return new double[]{lsum, tot > 0 ? 1.0 : -1.0};
    }

    /** Signed log-domain linear convolution truncated at the common length. */
    private static void slogconv(double[] lu, double[] su, double[] lv, double[] sv,
                                 double[] cond) {
        int n = lu.length;
        double[] lo = new double[n];
        double[] so = new double[n];
        double[] tl = new double[n];
        double[] ts = new double[n];
        for (int k = 0; k < n; k++) {
            for (int j = 0; j <= k; j++) {
                tl[j] = lu[j] + lv[k - j];
                ts[j] = su[j] * sv[k - j];
            }
            double[] r = slogsum(tl, ts, k + 1, cond);
            lo[k] = r[0];
            so[k] = r[1];
        }
        System.arraycopy(lo, 0, lu, 0, n);
        System.arraycopy(so, 0, su, 0, n);
    }

    /** log C(n,r); never a factorial quotient. */
    private static double lbinom(double n, double r) {
        return Gamma.logGamma(n + 1.0) - Gamma.logGamma(r + 1.0) - Gamma.logGamma(n - r + 1.0);
    }

    /** Nonnegative integer rows of length parts summing to total. */
    private static List<int[]> compositions(int total, int parts) {
        List<int[]> out = new ArrayList<int[]>();
        if (parts == 0) {
            if (total == 0) {
                out.add(new int[0]);
            }
            return out;
        }
        if (parts == 1) {
            out.add(new int[]{total});
            return out;
        }
        for (int first = 0; first <= total; first++) {
            List<int[]> sub = compositions(total - first, parts - 1);
            for (int i = 0; i < sub.size(); i++) {
                int[] s = sub.get(i);
                int[] row = new int[parts];
                row[0] = first;
                System.arraycopy(s, 0, row, 1, parts - 1);
                out.add(row);
            }
        }
        return out;
    }

    /**
     * Fuse PROPORTIONAL affine forms into one factor of summed multiplicity.
     *
     * <p>Two forms name the same pole exactly when they are proportional, and
     * that is the degeneracy Harrison-Coury Thm 1 excludes and its Conclusion
     * leaves open; merging disposes of it with nothing else changed.</p>
     */
    private static Object[] merge(double[][] F, int[] m, double tol) {
        int T = F.length;
        boolean[] used = new boolean[T];
        List<double[]> outF = new ArrayList<double[]>();
        List<Integer> outm = new ArrayList<Integer>();
        double lc = 0.0;
        double sc = 1.0;
        for (int i = 0; i < T; i++) {
            if (used[i]) {
                continue;
            }
            used[i] = true;
            double[] Fi = F[i];
            int mi = m[i];
            int pi = 0;
            for (int c = 1; c < Fi.length; c++) {
                if (FastMath.abs(Fi[c]) > FastMath.abs(Fi[pi])) {
                    pi = c;
                }
            }
            if (Fi[pi] == 0.0) {
                throw new IllegalStateException("pfqn_rgfmc met an identically zero factor.");
            }
            for (int j = i + 1; j < T; j++) {
                if (used[j]) {
                    continue;
                }
                double[] Fj = F[j];
                double nj = 0.0;
                for (int c = 0; c < Fj.length; c++) {
                    nj = FastMath.max(nj, FastMath.abs(Fj[c]));
                }
                if (nj <= 0.0) {
                    continue;
                }
                double r = Fj[pi] / Fi[pi];
                if (r == 0.0) {
                    continue;
                }
                boolean prop = true;
                for (int c = 0; c < Fj.length; c++) {
                    if (FastMath.abs(Fj[c] - r * Fi[c]) > tol * nj) {
                        prop = false;
                        break;
                    }
                }
                if (prop) {
                    lc -= m[j] * FastMath.log(FastMath.abs(r));
                    if (r < 0 && (m[j] % 2) == 1) {
                        sc = -sc;
                    }
                    mi += m[j];
                    used[j] = true;
                }
            }
            outF.add(Fi);
            outm.add(Integer.valueOf(mi));
        }
        double[][] nF = new double[outF.size()][];
        int[] nm = new int[outm.size()];
        for (int i = 0; i < outF.size(); i++) {
            nF[i] = outF.get(i);
            nm[i] = outm.get(i).intValue();
        }
        return new Object[]{nF, nm, Double.valueOf(lc), Double.valueOf(sc)};
    }

    private static double signPow(double x, int k) {
        return (x < 0 && (k % 2) == 1) ? -1.0 : 1.0;
    }

    /**
     * Eliminate the class in column col of F: Harrison-Coury Thm 1 as a partial
     * fraction, the poles of z_col sitting at A_j/B_j with order m_j.
     */
    private static List<Term> step(List<Term> terms, int col, int kr, double Zr,
                                   double tol, int maxterms, double[] cond) {
        List<Term> out = new ArrayList<Term>();
        for (int it = 0; it < terms.size(); it++) {
            Term t = terms.get(it);
            Object[] mg = merge(t.F, t.m, tol);
            double[][] F = (double[][]) mg[0];
            int[] m = (int[]) mg[1];
            double lc = t.lc + ((Double) mg[2]).doubleValue();
            double sc = t.sc * ((Double) mg[3]).doubleValue();
            int T = F.length;
            double[][] A = new double[T][col];
            double[] B = new double[T];
            boolean[] isMono = new boolean[T];
            int shift = 0;
            for (int i = 0; i < T; i++) {
                double scale = 0.0;
                for (int c = 0; c < F[i].length; c++) {
                    scale = FastMath.max(scale, FastMath.abs(F[i][c]));
                }
                if (scale == 0.0) {
                    scale = 1.0;
                }
                boolean mono = true;
                for (int c = 0; c < col; c++) {
                    A[i][c] = F[i][c];
                    if (FastMath.abs(F[i][c]) > tol * scale) {
                        mono = false;
                    }
                }
                B[i] = -F[i][col];
                isMono[i] = mono;
                if (mono) {
                    lc -= m[i] * FastMath.log(FastMath.abs(B[i]));
                    sc *= signPow(-B[i], m[i]);
                    shift += m[i];
                }
            }
            List<Integer> S = new ArrayList<Integer>();
            List<Integer> P = new ArrayList<Integer>();
            for (int i = 0; i < T; i++) {
                if (isMono[i]) {
                    continue;
                }
                if (B[i] != 0.0) {
                    S.add(Integer.valueOf(i));
                } else {
                    P.add(Integer.valueOf(i));
                }
            }
            int Ntot = kr + shift;
            int pmax = Zr > 0 ? Ntot : 0;
            for (int p = 0; p <= pmax; p++) {
                // Bertozzi-McKenna truncation, in logs: naive Z^p/p! overflows.
                double lcz = lc;
                if (p > 0) {
                    lcz += p * FastMath.log(Zr) - Gamma.logGamma(p + 1.0);
                }
                int n = Ntot - p;
                if (S.isEmpty()) {
                    if (n == 0) {
                        double[][] nF = new double[P.size()][];
                        int[] nm = new int[P.size()];
                        for (int a = 0; a < P.size(); a++) {
                            nF[a] = A[P.get(a).intValue()];
                            nm[a] = m[P.get(a).intValue()];
                        }
                        out.add(new Term(lcz, sc, nF, nm));
                    }
                    continue;
                }
                for (int jj = 0; jj < S.size(); jj++) {
                    int j = S.get(jj).intValue();
                    List<Integer> oth = new ArrayList<Integer>();
                    for (int a = 0; a < S.size(); a++) {
                        if (S.get(a).intValue() != j) {
                            oth.add(S.get(a));
                        }
                    }
                    int no = oth.size();
                    double Bj = B[j];
                    double[] Aj = A[j];
                    double[][] Cjl = new double[no][col];
                    for (int a = 0; a < no; a++) {
                        int l = oth.get(a).intValue();
                        for (int c = 0; c < col; c++) {
                            Cjl[a][c] = (A[l][c] * Bj - B[l] * Aj[c]) / Bj;
                        }
                    }
                    for (int k = 0; k < m[j]; k++) {
                        double lbase = lcz - k * FastMath.log(FastMath.abs(Bj))
                                + lbinom(n + m[j] - k - 1, n)
                                + n * FastMath.log(FastMath.abs(Bj));
                        double sbase = sc * signPow(-Bj, k) * signPow(Bj, n);
                        List<int[]> comps = compositions(k, no);
                        for (int cc = 0; cc < comps.size(); cc++) {
                            int[] jl = comps.get(cc);
                            double lt = lbase;
                            double st = sbase;
                            for (int a = 0; a < no; a++) {
                                if (jl[a] > 0) {
                                    int l = oth.get(a).intValue();
                                    lt += lbinom(m[l] + jl[a] - 1, jl[a])
                                            + jl[a] * FastMath.log(FastMath.abs(B[l]));
                                    st *= signPow(B[l], jl[a]);
                                }
                            }
                            int nn = 1 + no + P.size();
                            double[][] nF = new double[nn][];
                            int[] nm = new int[nn];
                            nF[0] = Aj;
                            nm[0] = n + m[j] - k;
                            for (int a = 0; a < no; a++) {
                                nF[1 + a] = Cjl[a];
                                nm[1 + a] = m[oth.get(a).intValue()] + jl[a];
                            }
                            for (int a = 0; a < P.size(); a++) {
                                nF[1 + no + a] = A[P.get(a).intValue()];
                                nm[1 + no + a] = m[P.get(a).intValue()];
                            }
                            out.add(new Term(lt, st, nF, nm));
                        }
                    }
                }
            }
            if (out.size() > maxterms) {
                throw new IllegalStateException(String.format(
                        "pfqn_rgfmc exceeded maxterms=%d. The residue term count grows as "
                        + "C(S+M-1,M-1) per further elimination; use method 'ca'.",
                        Integer.valueOf(maxterms)));
            }
        }
        return out;
    }

    /** Single-class base case, memoised on (loads, multiplicities) per Sec. 3.4. */
    private static double[] base(List<Term> terms, int k1, double Z1, double tol,
                                 double[] cond) {
        Map<String, double[]> cache = new HashMap<String, double[]>();
        double[] lv = new double[terms.size()];
        double[] sv = new double[terms.size()];
        int cnt = 0;
        for (int it = 0; it < terms.size(); it++) {
            Term t = terms.get(it);
            Object[] mg = merge(t.F, t.m, tol);
            double[][] F = (double[][]) mg[0];
            int[] m = (int[]) mg[1];
            double lc = t.lc + ((Double) mg[2]).doubleValue();
            double sc = t.sc * ((Double) mg[3]).doubleValue();
            int shift = 0;
            List<Double> loads = new ArrayList<Double>();
            List<Integer> mults = new ArrayList<Integer>();
            for (int a = 0; a < F.length; a++) {
                double a0 = F[a][0];
                double a1 = F[a][1];
                int ma = m[a];
                double sca = FastMath.max(FastMath.abs(a0), FastMath.abs(a1));
                if (sca == 0.0) {
                    throw new IllegalStateException("pfqn_rgfmc met an identically zero factor.");
                }
                if (FastMath.abs(a0) <= tol * sca) {
                    lc -= ma * FastMath.log(FastMath.abs(a1));
                    sc *= signPow(a1, ma);
                    shift += ma;
                } else {
                    lc -= ma * FastMath.log(FastMath.abs(a0));
                    sc *= signPow(a0, ma);
                    if (FastMath.abs(a1) > tol * sca) {
                        loads.add(Double.valueOf(-a1 / a0));
                        mults.add(Integer.valueOf(ma));
                    }
                }
            }
            int Ntot = k1 + shift;
            StringBuilder key = new StringBuilder();
            key.append(Ntot).append('|').append(Z1);
            for (int a = 0; a < loads.size(); a++) {
                key.append('|').append(loads.get(a)).append(':').append(mults.get(a));
            }
            double[] hit = cache.get(key.toString());
            if (hit == null) {
                double[] ld = new double[loads.size()];
                int[] mu = new int[mults.size()];
                for (int a = 0; a < ld.length; a++) {
                    ld[a] = loads.get(a).doubleValue();
                    mu[a] = mults.get(a).intValue();
                }
                hit = baseKernel(ld, mu, Ntot, Z1, cond);
                cache.put(key.toString(), hit);
            }
            if (hit[1] != 0.0) {
                lv[cnt] = lc + hit[0];
                sv[cnt] = sc * hit[1];
                cnt++;
            }
        }
        return slogsum(lv, sv, cnt, cond);
    }

    /**
     * [z^N] exp(Z z) prod_t (1 - p_t z)^-m_t, signed and in the log domain.
     *
     * <p>Coury-Harrison (1997) Property 1 with the loads allowed to be negative:
     * an eliminated class leaves pole differences that are not sign-definite.</p>
     */
    private static double[] baseKernel(double[] loads, int[] mults, int N, double Z,
                                       double[] cond) {
        double[] lg = new double[N + 1];
        double[] sg = new double[N + 1];
        Arrays.fill(lg, Double.NEGATIVE_INFINITY);
        lg[0] = 0.0;
        sg[0] = 1.0;
        double[] lr = new double[N + 1];
        double[] sr = new double[N + 1];
        if (Z > 0) {
            for (int k = 0; k <= N; k++) {
                lr[k] = k * FastMath.log(Z) - Gamma.logGamma(k + 1.0);
                sr[k] = 1.0;
            }
            slogconv(lg, sg, lr, sr, cond);
        }
        for (int a = 0; a < loads.length; a++) {
            double p = loads[a];
            int mm = mults[a];
            if (p == 0.0) {
                continue;
            }
            for (int k = 0; k <= N; k++) {
                if (mm == 1) {
                    lr[k] = k * FastMath.log(FastMath.abs(p));
                } else {
                    lr[k] = Gamma.logGamma(k + (double) mm) - Gamma.logGamma(k + 1.0)
                            - Gamma.logGamma((double) mm) + k * FastMath.log(FastMath.abs(p));
                }
                sr[k] = (p > 0 || (k % 2) == 0) ? 1.0 : -1.0;
            }
            slogconv(lg, sg, lr, sr, cond);
        }
        return new double[]{lg[N], sg[N]};
    }
}
