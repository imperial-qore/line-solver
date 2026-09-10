/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Gerasimov's residue (closed-form) normalizing constant, generalized to R classes
 *
 * A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing Networks",
 * Operations Research 43(4):704-711, 1995. Ported at parity from
 * matlab/src/api/pfqn/pfqn_gerasimov.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_gerasimov {

    private Pfqn_gerasimov() {}

    /** A product of powers of affine forms, carrying a scalar coefficient. */
    private static final class Term {
        double c;
        double[][] F;   // form j is F[j][0] + sum_s F[j][s] u_s
        int[] m;        // multiplicity of form j

        Term(double c, double[][] F, int[] m) {
            this.c = c;
            this.F = F;
            this.m = m;
        }
    }

    /**
     * Exact normalizing constant of a closed multiclass product-form network by ITERATED
     * RESIDUES of its rational generating function, one class at a time.
     *
     * <p>Gerasimov (1995) evaluates
     *
     * <pre>
     *   G(N_1,...,N_R) = (2 pi i)^-R int_G1 ... int_GR
     *                      prod_s z_s^(N_s-1) prod_i (1 - sum_s x_is/z_s)^-1
     * </pre>
     *
     * by residues, and gives the resulting CLOSED FORM only for R = 1 (Thm 1-2) and R = 2
     * (Thm 3 for simple poles, Thm 4 for multiple ones), stating that "for three or more
     * classes of customers, the normalizing constants can be found by numerical methods".
     * This routine implements the residue elimination itself, so the closed form is
     * produced for ANY R; at R = 2 it reproduces Thm 3/4 term by term.
     *
     * <p>Written as a coefficient of the u_s = 1/z_s series,
     *
     * <pre>
     *   G(N) = [prod_s u_s^(N_s)] exp(sum_s Z_s u_s) prod_i (1 - sum_s x_is u_s)^-1,
     * </pre>
     *
     * every factor is AFFINE in u, so singling out u_r gives f = A - B u_r with A affine in
     * the surviving variables. Partial fractions in u_r,
     *
     * <pre>
     *   [u_r^n] prod_j (A_j - B_j u_r)^-m_j
     *     = sum_j sum_{k=0}^{m_j-1} (-B_j)^-k C(n+m_j-k-1,n) B_j^n A_j^-(n+m_j-k)
     *       [t^k] prod_{l!=j} (C_jl - B_l t)^-m_l,   C_jl = (A_l B_j - B_l A_j)/B_j,
     * </pre>
     *
     * map a sum of products of affine powers into another one with one variable fewer, and
     * R-1 such steps leave a univariate coefficient extraction. At R = 2 the single step
     * returns one term per station i, with outer factor x_i2^(N_2+M-1)/prod_(k!=i)(x_i2-x_k2),
     * a pole of order N_2+1 at x_i1 and simple poles at the paper's
     * z_1ik = (x_k1 x_i2 - x_i1 x_k2)/(x_i2 - x_k2): exactly Thm 3, with the multiple poles
     * of Thm 4 (his xi_i &lt; M) handled by the same step. Tied x_i2, vanishing x_i2 and
     * identical station rows, all outside the paper's hypotheses, are ordinary cases here.
     *
     * <p>Cost. Let M be the number of stations and order the populations
     * N_(1) &lt;= ... &lt;= N_(R). The first elimination turns the single input term into M,
     * and every later one multiplies the count by C(S+M-1,M-1) + M-1, where S is the total
     * population already eliminated: a pole of order S+1 has to be differentiated against
     * the M-1 remaining ones. The innermost extraction then convolves M series of length
     * N_(1). Hence R = 1 costs O(M N), Buzen's own cost; R = 2 costs O(M^2 N_(1)^2),
     * INDEPENDENT OF N_(2); and R &gt;= 3 costs the same times prod_(r=3..R) C(N_(r)+M-1,M-1).
     * The R = 2 line is the reason to reach for this method: a population removed by
     * residues enters only as a pole ORDER, i.e. through binomial coefficients, so it costs
     * nothing at all. On a 4-station two-class model at N = [6, 20000] this returns lG in
     * 0.1 ms where pfqn_ca needs 79 s, to the same 1.3e-16. For R &gt;= 3 the term count is
     * polynomial in the populations of degree (M-1)(R-2) and exponential in R, which is why
     * the paper stops at two classes and why maxterms exists.
     *
     * <p>Conditioning. The sum is alternating, exactly as the paper writes it, and two
     * decisions keep it usable: near-coincident poles are merged under a RELATIVE tolerance,
     * so they are one multiple pole rather than two nearly cancelling simple ones, and the
     * class left for the innermost extraction is the one with the SMALLEST population,
     * because that population is the degree the final, sign-indefinite series is carried to.
     * Measured on 372 random models against pfqn_ca: median 2.0e-16, p90 4.6e-15, p99
     * 1.1e-11, worst 1.5e-09. On an ill-conditioned demand matrix pfqn_ca or pfqn_nc are
     * still the safer routes to the same number.
     *
     * @param L service demand matrix (M x R), L(i,r) = demand of class r at station i
     * @param N population vector (1 x R), nonnegative integers
     * @param Z think time vector (1 x R), may be null. A delay contributes the entire factor
     *          exp(sum_s Z_s u_s), handled exactly by convolving its Poisson coefficients
     *          into each elimination.
     * @param tol relative tolerance for declaring two affine forms proportional, hence one
     *          pole rather than two
     * @param maxterms cap on the number of residue terms carried between eliminations.
     *          Exceeding it is an error, not a truncation: a truncated residue sum is not a
     *          bound or an approximation of G, it is a wrong number.
     * @return the normalizing constant and its logarithm
     */
    public static Ret.pfqnNc pfqn_gerasimov(Matrix L, Matrix N, Matrix Z, double tol, int maxterms) {
        int R = L.getNumCols();
        double[] nn = new double[R];
        double[] zz = new double[R];
        if (N.length() != R) {
            throw new IllegalArgumentException("pfqn_gerasimov requires numel(N) to match the number of columns of L.");
        }
        for (int r = 0; r < R; r++) {
            nn[r] = N.get(r);
            zz[r] = (Z == null || Z.isEmpty()) ? 0.0 : Z.get(r);
        }
        for (int r = 0; r < R; r++) {
            if (nn[r] < 0 || zz[r] < 0) {
                throw new IllegalArgumentException("pfqn_gerasimov requires nonnegative L, N and Z.");
            }
            if (nn[r] != Math.rint(nn[r])) {
                throw new IllegalArgumentException("pfqn_gerasimov requires integer populations.");
            }
        }
        int M0 = L.getNumRows();
        for (int i = 0; i < M0; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) < 0) {
                    throw new IllegalArgumentException("pfqn_gerasimov requires nonnegative L, N and Z.");
                }
            }
        }

        // A class with no jobs is eliminated by evaluating the generating function at
        // u_r = 0, i.e. by deleting its column outright.
        List<Integer> cls = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (nn[r] > 0) {
                cls.add(Integer.valueOf(r));
            }
        }
        if (cls.isEmpty()) {
            return new Ret.pfqnNc(Double.valueOf(1.0), Double.valueOf(0.0));
        }
        // Class order, which decides both the cost and the accuracy.
        //  - The class left for the innermost extraction sets the CONDITIONING. Its
        //    population is the degree the final series is carried to, and the poles of the
        //    reduced problem have arbitrary sign, so that series cancels; the populations
        //    eliminated by residues enter only as pole ORDERS, through binomial
        //    coefficients, and cancel nothing. Basing on N = 100 rather than on N = 6 in
        //    one 4-station model cost 39 nats of lG. The SMALLEST population goes to base.
        //  - Eliminating class r leaves a pole of order N_r+1 that every LATER elimination
        //    has to differentiate, so the rest are eliminated smallest-first to keep the
        //    multiplicities low for as long as possible.
        // Eliminations run from index R down to 2, so indices 2..R hold the remaining
        // populations in DECREASING order and index 1 holds the smallest.
        for (int a = 0; a < cls.size(); a++) {
            for (int b = a + 1; b < cls.size(); b++) {
                if (nn[cls.get(b).intValue()] < nn[cls.get(a).intValue()]) {
                    Integer t = cls.get(a);
                    cls.set(a, cls.get(b));
                    cls.set(b, t);
                }
            }
        }
        for (int a = 1, b = cls.size() - 1; a < b; a++, b--) {
            Integer t = cls.get(a);
            cls.set(a, cls.get(b));
            cls.set(b, t);
        }
        int Rk = cls.size();
        int[] Nk = new int[Rk];
        double[] Zk = new double[Rk];
        for (int r = 0; r < Rk; r++) {
            Nk[r] = (int) Math.rint(nn[cls.get(r).intValue()]);
            Zk[r] = zz[cls.get(r).intValue()];
        }
        // Per-class scaling. The residue coefficients carry x_ir^(N_r+M-1), which in double
        // overflows well before G itself does: at x = 4 and N_r = 400 the factor alone is
        // 1e240 while G is finite. Dividing column r by c_r divides G by exactly c_r^N_r
        // (substitute u_r -> u_r/c_r in the generating function), so the scaling is exact
        // and is undone in the log domain at the end.
        double[] cs = new double[Rk];
        for (int r = 0; r < Rk; r++) {
            double c = Zk[r];
            for (int i = 0; i < M0; i++) {
                c = Math.max(c, L.get(i, cls.get(r).intValue()));
            }
            cs[r] = c > 0 ? c : 1.0;
            Zk[r] /= cs[r];
        }
        double lGscale = 0;
        for (int r = 0; r < Rk; r++) {
            lGscale += Nk[r] * Math.log(cs[r]);
        }
        // A station with no demand at all contributes the factor 1.
        List<double[]> rows = new ArrayList<double[]>();
        for (int i = 0; i < M0; i++) {
            double[] row = new double[Rk + 1];
            row[0] = 1.0;
            boolean any = false;
            for (int r = 0; r < Rk; r++) {
                double v = L.get(i, cls.get(r).intValue());
                row[r + 1] = -v / cs[r];
                if (v > 0) {
                    any = true;
                }
            }
            if (any) {
                rows.add(row);
            }
        }
        int M = rows.size();
        double[][] F0 = new double[M][];
        int[] m0 = new int[M];
        for (int i = 0; i < M; i++) {
            F0[i] = rows.get(i);
            m0[i] = 1;
        }

        List<Term> terms = new ArrayList<Term>();
        terms.add(new Term(1.0, F0, m0));
        // class s (1-based) sits in column s of F = [1, -L]; eliminate R, R-1, ..., 2
        for (int r = Rk; r >= 2; r--) {
            terms = step(terms, r, Nk[r - 1], Zk[r - 1], tol, maxterms);
            if (terms.isEmpty()) {
                return new Ret.pfqnNc(Double.valueOf(0.0), Double.valueOf(Double.NEGATIVE_INFINITY));
            }
        }
        double Gs = base(terms, Nk[0], Zk[0], tol);
        if (Gs <= 0) {
            return new Ret.pfqnNc(Double.valueOf(0.0), Double.valueOf(Double.NEGATIVE_INFINITY));
        }
        double lGout = Math.log(Gs) + lGscale;
        // exp(lGout), not Gs*exp(lGscale): the latter overflows forming
        // exp(lGscale) alone even when the product is inside the double range.
        return new Ret.pfqnNc(Double.valueOf(Math.exp(lGout)), Double.valueOf(lGout));
    }

    /** Defaults: tol = 1e-12, maxterms = 200000. */
    public static Ret.pfqnNc pfqn_gerasimov(Matrix L, Matrix N, Matrix Z) {
        return pfqn_gerasimov(L, N, Z, 1e-12, 200000);
    }

    /** Delay-free case. */
    public static Ret.pfqnNc pfqn_gerasimov(Matrix L, Matrix N) {
        return pfqn_gerasimov(L, N, null, 1e-12, 200000);
    }

    // ------------------------------------------------------------------
    // One residue elimination: integrate out u_r and return the surviving sum of products
    // of affine powers, each form narrowed from r+1 to r columns.
    // ------------------------------------------------------------------
    private static List<Term> step(List<Term> terms, int r, int Nr, double Zr, double tol, int maxterms) {
        List<Term> out = new ArrayList<Term>();
        for (int it = 0; it < terms.size(); it++) {
            Term t = merge(terms.get(it), tol);
            int nf = t.F.length;
            double[][] A = new double[nf][r];
            double[] B = new double[nf];
            double[] scale = new double[nf];
            for (int j = 0; j < nf; j++) {
                for (int s = 0; s < r; s++) {
                    A[j][s] = t.F[j][s];
                }
                B[j] = -t.F[j][r];
                double sc = 0;
                for (int s = 0; s <= r; s++) {
                    sc = Math.max(sc, Math.abs(t.F[j][s]));
                }
                scale[j] = sc;
            }
            double c = t.c;
            int shift = 0;
            boolean[] drop = new boolean[nf];
            for (int j = 0; j < nf; j++) {
                boolean mono = true;
                for (int s = 0; s < r; s++) {
                    if (Math.abs(A[j][s]) > tol * scale[j]) {
                        mono = false;
                        break;
                    }
                }
                if (mono && Math.abs(B[j]) <= tol * scale[j]) {
                    throw new IllegalStateException("pfqn_gerasimov met an identically zero factor, which cannot happen after merging proportional ones.");
                }
                if (mono) {
                    // A factor -B u_r carries no finite pole: it only shifts the exponent.
                    c *= Math.pow(-B[j], -t.m[j]);
                    shift += t.m[j];
                    drop[j] = true;
                }
            }
            List<Integer> S = new ArrayList<Integer>();  // factors carrying a pole in u_r
            List<Integer> P = new ArrayList<Integer>();  // factors free of u_r
            for (int j = 0; j < nf; j++) {
                if (drop[j]) {
                    continue;
                }
                if (B[j] != 0.0) {
                    S.add(Integer.valueOf(j));
                } else {
                    P.add(Integer.valueOf(j));
                }
            }
            int Ntot = Nr + shift;
            int pmax = Zr > 0 ? Ntot : 0;
            for (int p = 0; p <= pmax; p++) {
                double cz = c;
                if (p > 0) {
                    // Poisson weight Z_r^p/p! through logs: the naive ratio overflows
                    // for p >~ 171, reachable when the eliminated class has think time.
                    cz = c * Math.exp(p * Math.log(Zr) - Gamma.logGamma(p + 1.0));
                }
                int Neff = Ntot - p;
                if (S.isEmpty()) {
                    if (Neff == 0) {
                        double[][] Fn = new double[P.size()][];
                        int[] mn = new int[P.size()];
                        for (int q = 0; q < P.size(); q++) {
                            Fn[q] = A[P.get(q).intValue()].clone();
                            mn[q] = t.m[P.get(q).intValue()];
                        }
                        out.add(new Term(cz, Fn, mn));
                    }
                    continue;
                }
                for (int jj = 0; jj < S.size(); jj++) {
                    int j = S.get(jj).intValue();
                    List<Integer> oth = new ArrayList<Integer>();
                    for (int q = 0; q < S.size(); q++) {
                        if (q != jj) {
                            oth.add(S.get(q));
                        }
                    }
                    int no = oth.size();
                    double[][] Cjl = new double[no][r];
                    for (int l = 0; l < no; l++) {
                        int k = oth.get(l).intValue();
                        for (int s = 0; s < r; s++) {
                            Cjl[l][s] = (A[k][s] * B[j] - B[k] * A[j][s]) / B[j];
                        }
                    }
                    for (int k = 0; k < t.m[j]; k++) {
                        List<int[]> comps = compositions(k, no);
                        for (int ci = 0; ci < comps.size(); ci++) {
                            int[] nk = comps.get(ci);
                            double coef = cz * Math.pow(-B[j], -k)
                                    * binom(Neff + t.m[j] - k - 1, Neff) * Math.pow(B[j], Neff);
                            for (int l = 0; l < no; l++) {
                                int kk = oth.get(l).intValue();
                                coef *= binom(t.m[kk] + nk[l] - 1, nk[l]) * Math.pow(B[kk], nk[l]);
                            }
                            if (coef == 0.0) {
                                continue;
                            }
                            int nn2 = 1 + no + P.size();
                            double[][] Fn = new double[nn2][];
                            int[] mn = new int[nn2];
                            Fn[0] = A[j].clone();
                            mn[0] = Neff + t.m[j] - k;
                            for (int l = 0; l < no; l++) {
                                Fn[1 + l] = Cjl[l].clone();
                                mn[1 + l] = t.m[oth.get(l).intValue()] + nk[l];
                            }
                            for (int q = 0; q < P.size(); q++) {
                                Fn[1 + no + q] = A[P.get(q).intValue()].clone();
                                mn[1 + no + q] = t.m[P.get(q).intValue()];
                            }
                            out.add(new Term(coef, Fn, mn));
                        }
                    }
                }
            }
            if (out.size() > maxterms) {
                throw new IllegalStateException("pfqn_gerasimov exceeded maxterms (" + maxterms
                        + ") while eliminating class " + r
                        + "; the residue expansion of this model is too large. Use pfqn_ca or pfqn_nc.");
            }
        }
        return out;
    }

    // ------------------------------------------------------------------
    // Last class: a univariate coefficient extraction. Summing the residues here too would
    // repeat the step above, but convolving the series of each factor returns the same
    // number without expanding the multiple poles.
    // ------------------------------------------------------------------
    private static double base(List<Term> terms, int N1, double Z1, double tol) {
        double G = 0;
        for (int it = 0; it < terms.size(); it++) {
            Term t = merge(terms.get(it), tol);
            int nf = t.F.length;
            double c = t.c;
            int shift = 0;
            boolean[] drop = new boolean[nf];
            for (int j = 0; j < nf; j++) {
                double scale = Math.max(Math.abs(t.F[j][0]), Math.abs(t.F[j][1]));
                double Aj = t.F[j][0];
                double Bj = -t.F[j][1];
                if (Math.abs(Aj) <= tol * scale) {
                    if (Math.abs(Bj) <= tol * scale) {
                        throw new IllegalStateException("pfqn_gerasimov met an identically zero factor at the innermost coefficient extraction.");
                    }
                    c *= Math.pow(-Bj, -t.m[j]);
                    shift += t.m[j];
                    drop[j] = true;
                }
            }
            int Ntot = N1 + shift;
            double[] s = new double[Ntot + 1];
            s[0] = 1.0;
            if (Z1 > 0) {
                double[] pois = new double[Ntot + 1];
                for (int n = 0; n <= Ntot; n++) {
                    pois[n] = Math.exp(n * Math.log(Z1) - Gamma.logGamma(n + 1.0));
                }
                s = conv(s, pois, Ntot);
            }
            for (int j = 0; j < nf; j++) {
                if (drop[j]) {
                    continue;
                }
                double Aj = t.F[j][0];
                double Bj = -t.F[j][1];
                c *= Math.pow(Aj, -t.m[j]);
                if (Bj == 0.0) {
                    continue;
                }
                double ratio = Bj / Aj;
                double[] seq = new double[Ntot + 1];
                for (int n = 0; n <= Ntot; n++) {
                    seq[n] = binom(t.m[j] + n - 1, n) * Math.pow(ratio, n);
                }
                s = conv(s, seq, Ntot);
            }
            G += c * s[Ntot];
        }
        return G;
    }

    /** Truncated convolution of two coefficient sequences, kept to degree n. */
    private static double[] conv(double[] a, double[] b, int n) {
        double[] y = new double[n + 1];
        for (int i = 0; i < a.length && i <= n; i++) {
            if (a[i] == 0.0) {
                continue;
            }
            for (int j = 0; j < b.length && i + j <= n; j++) {
                y[i + j] += a[i] * b[j];
            }
        }
        return y;
    }

    /**
     * Merge proportional affine forms: f_k = lambda f_j is one pole of order m_j+m_k, not
     * two nearby simple ones, and lambda^-m_k moves into the scalar.
     */
    private static Term merge(Term t, double tol) {
        int nf = t.F.length;
        if (nf <= 1) {
            return t;
        }
        int w = t.F[0].length;
        boolean[] keep = new boolean[nf];
        int[] m = t.m.clone();
        double c = t.c;
        for (int j = 0; j < nf; j++) {
            keep[j] = true;
        }
        for (int j = 0; j < nf; j++) {
            if (!keep[j]) {
                continue;
            }
            int pj = 0;
            for (int s = 1; s < w; s++) {
                if (Math.abs(t.F[j][s]) > Math.abs(t.F[j][pj])) {
                    pj = s;
                }
            }
            if (t.F[j][pj] == 0.0) {
                continue;
            }
            for (int k = j + 1; k < nf; k++) {
                if (!keep[k]) {
                    continue;
                }
                double lam = t.F[k][pj] / t.F[j][pj];
                if (lam == 0.0) {
                    continue;
                }
                double dev = 0;
                double scale = 0;
                for (int s = 0; s < w; s++) {
                    dev = Math.max(dev, Math.abs(t.F[k][s] - lam * t.F[j][s]));
                    scale = Math.max(scale, Math.max(Math.abs(t.F[k][s]), Math.abs(t.F[j][s])));
                }
                if (dev <= tol * scale) {
                    c *= Math.pow(lam, -m[k]);
                    m[j] += m[k];
                    keep[k] = false;
                }
            }
        }
        int nk = 0;
        for (int j = 0; j < nf; j++) {
            if (keep[j]) {
                nk++;
            }
        }
        if (nk == nf) {
            return new Term(c, t.F, m);
        }
        double[][] Fn = new double[nk][];
        int[] mn = new int[nk];
        int q = 0;
        for (int j = 0; j < nf; j++) {
            if (keep[j]) {
                Fn[q] = t.F[j];
                mn[q] = m[j];
                q++;
            }
        }
        return new Term(c, Fn, mn);
    }

    /** All k-tuples of nonnegative integers summing to n. */
    private static List<int[]> compositions(int n, int k) {
        List<int[]> out = new ArrayList<int[]>();
        if (k == 0) {
            if (n == 0) {
                out.add(new int[0]);
            }
            return out;
        }
        if (k == 1) {
            out.add(new int[] { n });
            return out;
        }
        for (int a = 0; a <= n; a++) {
            List<int[]> sub = compositions(n - a, k - 1);
            for (int i = 0; i < sub.size(); i++) {
                int[] row = new int[k];
                row[0] = a;
                System.arraycopy(sub.get(i), 0, row, 1, k - 1);
                out.add(row);
            }
        }
        return out;
    }

    /** Binomial coefficient, exact while the value stays representable. */
    private static double binom(int n, int k) {
        if (k < 0 || n < 0 || k > n) {
            return 0.0;
        }
        int kk = Math.min(k, n - k);
        double b = 1.0;
        for (int i = 1; i <= kk; i++) {
            b = b * (n - kk + i) / i;
        }
        if (b < 9007199254740992.0) {
            b = Math.rint(b);
        }
        return b;
    }

    private static double factorial(int n) {
        double f = 1.0;
        for (int i = 2; i <= n; i++) {
            f *= i;
        }
        return f;
    }
}
