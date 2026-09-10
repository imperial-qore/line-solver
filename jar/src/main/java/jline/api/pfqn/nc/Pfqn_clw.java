/**
 * Choudhury-Leung-Whitt normalization constant by numerical inversion of the
 * generating function (JACM 42(5):935-970, 1995).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Computes g(K) of a multichain closed product-form network with single-server
 * and (optionally) infinite-server queues by numerically inverting its
 * p-dimensional generating function (Choudhury, Leung and Whitt, 1995, eq. 4.5)
 *
 *   G(z) = exp(sum_j rho_{j0} z_j) / prod_i (1 - sum_j rho_{ji} z_j)^{m_i}
 *
 * where j=1..p indexes chains and i=1..q' the distinct single-server queues with
 * multiplicity m_i. g(K) is recovered by nested one-dimensional lattice-Poisson
 * inversions (eq. 2.3) with restrictive static scaling (eqs. 5.41-5.46) and
 * log-domain recovery (eq. 7.1).
 *
 * Both of the paper's accelerations are applied. Dimension reduction by
 * decomposition (Sec. 3, Sec. 5.4) reads the interdependence graph of the
 * factors of (4.5), removes the subset D minimizing |D| + max_i |S_i(D)|
 * (eq. 3.3), and inverts each connected component of the remainder separately,
 * multiplying the results. Euler summation (Sec. 2.4, eq. 2.22) replaces the
 * nearly alternating inner sum of (2.3) by the Euler sum of its first n+m+1
 * terms wherever K_j exceeds n+m, so prod_j K_j becomes prod_j min(n+m+1, K_j)
 * in the cost (eq. 2.26).
 */
public final class Pfqn_clw {
    private Pfqn_clw() {}

    /** Accelerations and their parameters; the defaults are the paper's. */
    public static final class Options {
        /** Apply Euler summation where K_j > eulerN + eulerM. */
        public boolean euler = true;
        /** Terms summed exactly before averaging (n in eq. 2.22). */
        public int eulerN = 11;
        /** Starting order of the Euler averaging (m in eq. 2.22); doubled on demand. */
        public int eulerM = 20;
        /** Relative tolerance on |E(m,n) - E(m,n+1)| for accepting the Euler sum. */
        public double eulerTol = 1e-10;
        /** Largest Euler order reached by doubling. */
        public int eulerMaxM = 160;
        /** Apply dimension reduction by decomposition (Section 3). */
        public boolean dimred = true;
        /** Largest |D| examined when minimizing (3.3). */
        public int dimredMaxD = 4;
        /** Multipliers on the scale parameters alpha_j, the manual tuning of page 956. */
        public double[] beta = null;
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N) {
        return pfqn_clw(L, N, null, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z) {
        return pfqn_clw(L, N, Z, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z, Matrix m) {
        return pfqn_clw(L, N, Z, m, null, null);
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z, Matrix m, Matrix lpar, Matrix gampar) {
        return pfqn_clw(L, N, Z, m, lpar, gampar, new Options());
    }

    /**
     * @param L     (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
     * @param N     (1 x p or p x 1) closed-chain population vector K.
     * @param Z     (1 x p) aggregate infinite-server relative intensities rho_{j0}; null = 0.
     * @param m     (q' x 1) queue multiplicities m_i; null = ones.
     * @param lpar  (1 x p) inner lattice parameters l_j indexed by chain; null = by
     *              inversion depth, 1 at depth 1, 2 at depths 2-3, 3 deeper.
     * @param gampar(1 x p) aliasing parameters gamma_j indexed by chain; null = by
     *              depth, 11, 13, 13, 15, ...
     * @param opt   accelerations; null = defaults.
     * @return normalization constant g(K) (Inf on overflow) and its natural log.
     */
    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z, Matrix m, Matrix lpar, Matrix gampar,
                                      Options opt) {
        int qd = L.getNumRows();
        int p = L.getNumCols();
        if (opt == null) {
            opt = new Options();
        }

        int[] Nv = new int[p];
        for (int j = 0; j < p; j++) {
            Nv[j] = (int) Math.round(N.get(j));
        }
        double[] Zv = new double[p];
        if (Z != null && !Z.isEmpty()) {
            for (int j = 0; j < p; j++) {
                Zv[j] = Z.get(j);
            }
        }
        double[] mv = new double[qd];
        if (m != null && !m.isEmpty()) {
            for (int i = 0; i < qd; i++) {
                mv[i] = m.get(i);
            }
        } else {
            Arrays.fill(mv, 1.0);
        }

        boolean anyNeg = false;
        long sumN = 0;
        for (int j = 0; j < p; j++) {
            if (Nv[j] < 0) anyNeg = true;
            sumN += Nv[j];
        }
        if (anyNeg) {
            return new Ret.pfqnNc(0.0, Double.NEGATIVE_INFINITY);
        }
        if (sumN == 0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        // read L into a dense array
        double[][] Lv = new double[qd][p];
        for (int i = 0; i < qd; i++) {
            for (int j = 0; j < p; j++) {
                Lv[i][j] = L.get(i, j);
            }
        }

        // dimension reduction (Section 3): D is inverted first, then each connected
        // component of the interdependence graph minus D, independently
        Decomposition dec = dimred(Lv, qd, p, opt);
        int[] ordD = dec.D;
        List<int[]> comps = dec.comps;
        int d = ordD.length;

        int[] order = new int[p];
        int pos = 0;
        for (int t = 0; t < d; t++) {
            order[pos++] = ordD[t];
        }
        for (int c = 0; c < comps.size(); c++) {
            int[] vc = comps.get(c);
            for (int s = 0; s < vc.length; s++) {
                order[pos++] = vc[s];
            }
        }

        // depth of each chain: D occupies depths 1..d and every component restarts
        // at depth d+1, since components are inverted in parallel
        int[] depth = new int[p];
        for (int t = 0; t < d; t++) {
            depth[ordD[t]] = t + 1;
        }
        for (int c = 0; c < comps.size(); c++) {
            int[] vc = comps.get(c);
            for (int s = 0; s < vc.length; s++) {
                depth[vc[s]] = d + s + 1;
            }
        }

        int[] lv = new int[p];
        if (lpar != null && !lpar.isEmpty()) {
            for (int j = 0; j < p; j++) {
                lv[j] = (int) Math.round(lpar.get(j));
            }
        } else {
            for (int j = 0; j < p; j++) {
                lv[j] = (depth[j] == 1) ? 1 : ((depth[j] <= 3) ? 2 : 3);
            }
        }
        double[] gamv = new double[p];
        if (gampar != null && !gampar.isEmpty()) {
            for (int j = 0; j < p; j++) {
                gamv[j] = gampar.get(j);
            }
        } else {
            for (int j = 0; j < p; j++) {
                gamv[j] = (depth[j] == 1) ? 11.0 : ((depth[j] <= 3) ? 13.0 : 15.0);
            }
        }

        // contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (eq. 2.7)
        double[] r = new double[p];
        for (int j = 0; j < p; j++) {
            r[j] = (Nv[j] == 0) ? 1.0
                    : Math.pow(10.0, -gamv[j] / (2.0 * lv[j] * Nv[j]));
        }

        // restrictive static scaling (eqs. 5.41-5.46), outer vars at |z_k| = r_k.
        // The loop runs in inversion order, which is what dimension reduction
        // changes (Section 5.4); chains of distinct components never deflate one
        // another, because they share no queue.
        double[] alpha = new double[p];
        Arrays.fill(alpha, 1.0);
        double[] used = new double[qd];
        for (int t = 0; t < p; t++) {
            int j = order[t];
            int Kj = Nv[j];
            int lj = lv[j];
            if (Kj == 0) { continue; }  // empty chain: 2*lj*Kj = 0 below
            int nc = 0;
            for (int i = 0; i < qd; i++) {
                if (Lv[i][j] > 0) nc++;
            }
            double aj = Double.POSITIVE_INFINITY;
            if (nc > 0) {
                Integer[] idx = new Integer[nc];
                double[] eAll = new double[qd];
                int c = 0;
                for (int i = 0; i < qd; i++) {
                    double denom = 1.0 - used[i];
                    if (denom <= 0) denom = Double.MIN_VALUE;
                    eAll[i] = Lv[i][j] / denom;
                    if (Lv[i][j] > 0) {
                        idx[c++] = i;
                    }
                }
                final double[] eRef = eAll;
                Arrays.sort(idx, (x, y) -> Double.compare(eRef[y], eRef[x]));  // descending
                double sumE = 0.0;
                double cummb = 0.0;
                for (int n = 0; n < nc; n++) {
                    int qi = idx[n];
                    sumE += eAll[qi];
                    cummb += mv[qi];
                    double cumrho = sumE / (n + 1);
                    // N_{ij} = mbar_n - 1 + sum_{k inside j} K_k eta_{k,qi} (eq. 5.43)
                    double Nnd = cummb - 1.0;
                    for (int u = t + 1; u < p; u++) {
                        int k = order[u];
                        if (Lv[qi][k] != 0) Nnd += Nv[k];
                    }
                    int Nn = (int) Math.round(Nnd);
                    double an;
                    if (Nn <= 0) {
                        an = 1.0;
                    } else {
                        // in the log domain: the product runs over N_{ij} factors
                        // below one and underflows to zero at a few hundred of
                        // them, which would silently set alpha_j = 0 and lG = NaN
                        double lp = 0.0;
                        double twolK = 2.0 * lj * Kj;
                        for (int ll = 1; ll <= Nn; ll++) {
                            lp += Math.log((Kj + ll) / (Kj + twolK + ll));
                        }
                        an = Math.exp(lp / twolK);
                    }
                    aj = Math.min(aj, an / cumrho);
                }
            }
            if (Zv[j] > 0) {
                aj = Math.min(aj, Kj / Zv[j]);
            }
            if (!Double.isFinite(aj)) {
                aj = 1.0;
            }
            alpha[j] = ((opt.beta != null) ? opt.beta[j] : 1.0) * aj;   // page 956
            for (int i = 0; i < qd; i++) {
                used[i] += alpha[j] * Lv[i][j] * r[j];
            }
        }

        // scaled process parameters
        double[] arho0 = new double[p];      // alpha_j rho_{j0}
        double[][] rhoS = new double[qd][p]; // alpha_j rho_{ji}
        for (int j = 0; j < p; j++) {
            arho0[j] = alpha[j] * Zv[j];
            for (int i = 0; i < qd; i++) {
                rhoS[i][j] = Lv[i][j] * alpha[j];
            }
        }

        // split the factors of (4.5) over D and the components: a queue whose
        // chains all lie in D is constant during the component inversions, and
        // every other queue has all of its non-D chains inside a single component
        boolean[] inD = new boolean[p];
        for (int t = 0; t < d; t++) {
            inD[ordD[t]] = true;
        }
        int[] compOf = new int[p];
        for (int c = 0; c < comps.size(); c++) {
            int[] vc = comps.get(c);
            for (int s = 0; s < vc.length; s++) {
                compOf[vc[s]] = c;
            }
        }
        int[] qBucket = new int[qd];
        for (int i = 0; i < qd; i++) {
            qBucket[i] = -1;
            for (int j = 0; j < p; j++) {
                if (Lv[i][j] != 0 && !inD[j]) {
                    qBucket[i] = compOf[j];
                    break;
                }
            }
        }
        int[] qD = collect(qBucket, -1);
        int[][] qC = new int[comps.size()][];
        for (int c = 0; c < comps.size(); c++) {
            qC[c] = collect(qBucket, c);
        }

        // per-group normalization (Section 2.2, page 944): the scaling normalizes
        // the whole generating function, not each group, and a decomposition
        // multiplies the groups together. Every factor has nonnegative
        // coefficients, so a group's modulus is maximized at w = r; that constant
        // cancels in the recovery and is left at zero on the undecomposed path.
        double[] off = new double[comps.size() + 1];   // [0] is the D group
        if (d > 0 || comps.size() > 1) {
            for (int b = 0; b < off.length; b++) {
                int tag = b - 1;
                int[] ch = (tag < 0) ? ordD : comps.get(tag);
                double o = 0.0;
                for (int t = 0; t < ch.length; t++) {
                    o += arho0[ch[t]] * (r[ch[t]] - 1.0);
                }
                for (int i = 0; i < qd; i++) {
                    if (qBucket[i] != tag) continue;
                    double x = 0.0;
                    for (int j = 0; j < p; j++) {
                        x += rhoS[i][j] * r[j];
                    }
                    double pole = 1.0 - x;
                    if (!(pole > 0)) pole = Double.MIN_VALUE;
                    o -= mv[i] * Math.log(pole);
                }
                off[b] = o;
            }
        }

        Ctx ctx = new Ctx(qd, p, Nv, lv, r, arho0, rhoS, mv);
        ctx.off = off;
        ctx.d = d;
        ctx.ordD = ordD;
        ctx.comps = comps;
        ctx.qD = qD;
        ctx.qC = qC;
        ctx.euler = opt.euler;
        ctx.eulerN = opt.eulerN;
        ctx.eulerM = opt.eulerM;
        ctx.eulerTol = opt.eulerTol;
        ctx.eulerMaxM = opt.eulerMaxM;

        double[] wRe = new double[p];
        double[] wIm = new double[p];
        double[] gbar = invertD(0, ctx, wRe, wIm);   // real part carries gbar(K)

        double lsum = 0.0;
        for (int b = 0; b < off.length; b++) {
            lsum += off[b];
        }
        for (int j = 0; j < p; j++) {
            lsum += arho0[j] - Nv[j] * Math.log(alpha[j]);
        }
        double lG = Math.log(gbar[0]) + lsum;
        double G = (lG > 709) ? Double.POSITIVE_INFINITY : Math.exp(lG);
        return new Ret.pfqnNc(G, lG);
    }

    private static int[] collect(int[] bucket, int want) {
        int n = 0;
        for (int i = 0; i < bucket.length; i++) {
            if (bucket[i] == want) n++;
        }
        int[] out = new int[n];
        int c = 0;
        for (int i = 0; i < bucket.length; i++) {
            if (bucket[i] == want) out[c++] = i;
        }
        return out;
    }

    // immutable context shared (read-only) across the recursion
    private static final class Ctx {
        final int qd, p;
        final int[] N, l;
        final double[] r, arho0;
        final double[][] rhoS;
        final double[] m;
        int d;
        int[] ordD;
        List<int[]> comps;
        int[] qD;
        int[][] qC;
        double[] off;
        boolean euler;
        int eulerN, eulerM, eulerMaxM;
        double eulerTol;
        Ctx(int qd, int p, int[] N, int[] l, double[] r, double[] arho0, double[][] rhoS, double[] m) {
            this.qd = qd; this.p = p; this.N = N; this.l = l;
            this.r = r; this.arho0 = arho0; this.rhoS = rhoS; this.m = m;
        }
    }

    // ---- outer inversion over the committed variables D (Section 3) ----
    private static double[] invertD(int t, Ctx ctx, double[] wRe, double[] wIm) {
        if (t >= ctx.d) {
            // D fixed: the remaining factors have no variable in common, so the
            // coefficient of the inner monomial is the product of the components'
            double[] val = gbarSub(ctx, ctx.qD, ctx.ordD, ctx.off[0], wRe, wIm);
            for (int c = 0; c < ctx.comps.size(); c++) {
                double[] cv = invertC(c, 0, ctx, wRe, wIm);
                double re = val[0] * cv[0] - val[1] * cv[1];
                double im = val[0] * cv[1] + val[1] * cv[0];
                val = new double[]{re, im};
            }
            return val;
        }
        int j = ctx.ordD[t];
        double[] val = lattice(j, ctx, wRe, wIm, t + 1, -1, -1);
        if (t == 0) {
            val[1] = 0.0;
        }
        return val;
    }

    // ---- inversion of one component of the interdependence graph minus D ----
    private static double[] invertC(int c, int s, Ctx ctx, double[] wRe, double[] wIm) {
        int[] vars = ctx.comps.get(c);
        if (s >= vars.length) {
            return gbarSub(ctx, ctx.qC[c], vars, ctx.off[c + 1], wRe, wIm);
        }
        double[] val = lattice(vars[s], ctx, wRe, wIm, -1, c, s + 1);
        if (ctx.d == 0 && s == 0) {
            // with D empty every component is an independent subnetwork, so its
            // coefficient is real
            val[1] = 0.0;
        }
        return val;
    }

    // one-dimensional lattice-Poisson inversion (eq. 2.3), scaled; returns {re,im}.
    // Exactly one of (nextD >= 0) and (comp >= 0) selects what is evaluated at the
    // contour points: the next D level, or the next level of component comp.
    private static double[] lattice(int j, Ctx ctx, double[] wRe, double[] wIm,
                                    int nextD, int comp, int nextS) {
        int Kj = ctx.N[j];
        int lj = ctx.l[j];
        double rj = ctx.r[j];
        if (Kj == 0) {
            // [w_j^0] Gbar = Gbar(w_j=0): the K=0 lattice is the single point 0, and
            // exp(-arho0_j) there cancels the +arho0_j added back into lG
            wRe[j] = 0.0;
            wIm[j] = 0.0;
            return (nextD >= 0) ? invertD(nextD, ctx, wRe, wIm)
                    : invertC(comp, nextS, ctx, wRe, wIm);
        }
        double accRe = 0.0, accIm = 0.0;
        for (int k1 = 0; k1 < lj; k1++) {
            double phRe = Math.cos(-Math.PI * k1 / lj);
            double phIm = Math.sin(-Math.PI * k1 / lj);
            double[] inner = innerSum(j, Kj, lj, k1, rj, ctx, wRe, wIm, nextD, comp, nextS);
            accRe += phRe * inner[0] - phIm * inner[1];
            accIm += phRe * inner[1] + phIm * inner[0];
        }
        double denom = 2.0 * lj * Kj * Math.pow(rj, Kj);
        return new double[]{accRe / denom, accIm / denom};
    }

    // ---- inner sum of (2.3) over the lattice index k, with Euler summation ----
    // The sum splits at k = 0 into two nearly alternating series (Section 2.4);
    // each is replaced by its Euler sum (eq. 2.22). The order m is doubled until
    // the paper's own estimate |E(m,n) - E(m,n+1)| falls under the tolerance, and
    // the exact sum is taken once n+m reaches K_j, so accuracy is not traded away.
    private static double[] innerSum(int j, int Kj, int lj, int k1, double rj, Ctx ctx,
                                     double[] wRe, double[] wIm, int nextD, int comp, int nextS) {
        int mCur = ctx.eulerM;
        while (true) {
            int T = ctx.eulerN + mCur;
            if (!ctx.euler || Kj <= T + 1) {
                double sRe = 0.0, sIm = 0.0;
                for (int k = -Kj; k < Kj; k++) {
                    double[] v = evalAt(j, k, Kj, lj, k1, rj, ctx, wRe, wIm, nextD, comp, nextS);
                    double sg = ((k % 2 == 0) ? 1.0 : -1.0);
                    sRe += sg * v[0];
                    sIm += sg * v[1];
                }
                return new double[]{sRe, sIm};
            }
            double[] w1 = eulerWeights(ctx.eulerN, mCur);        // s = 0..T
            double[] w2 = eulerWeights(ctx.eulerN + 1, mCur);    // s = 0..T+1
            double e1Re = 0.0, e1Im = 0.0, e2Re = 0.0, e2Im = 0.0;
            for (int s = 0; s <= T + 1; s++) {
                double[] vp = evalAt(j, s, Kj, lj, k1, rj, ctx, wRe, wIm, nextD, comp, nextS);
                double[] vn = evalAt(j, -(s + 1), Kj, lj, k1, rj, ctx, wRe, wIm, nextD, comp, nextS);
                double dRe = vp[0] - vn[0];
                double dIm = vp[1] - vn[1];
                double sg = ((s % 2 == 0) ? 1.0 : -1.0);
                if (s <= T) {
                    e1Re += sg * w1[s] * dRe;
                    e1Im += sg * w1[s] * dIm;
                }
                e2Re += sg * w2[s] * dRe;
                e2Im += sg * w2[s] * dIm;
            }
            if (Math.hypot(e1Re - e2Re, e1Im - e2Im) <= ctx.eulerTol * Math.hypot(e2Re, e2Im)
                    || mCur >= ctx.eulerMaxM) {
                return new double[]{e2Re, e2Im};
            }
            mCur *= 2;
        }
    }

    private static double[] evalAt(int j, int k, int Kj, int lj, int k1, double rj, Ctx ctx,
                                   double[] wRe, double[] wIm, int nextD, int comp, int nextS) {
        double theta = Math.PI * (k1 + (double) lj * k) / ((double) lj * Kj);
        wRe[j] = rj * Math.cos(theta);
        wIm[j] = rj * Math.sin(theta);
        return (nextD >= 0) ? invertD(nextD, ctx, wRe, wIm)
                : invertC(comp, nextS, ctx, wRe, wIm);
    }

    // ---- Euler weights: E(m,n) = sum_i w_i (-1)^i a_i, from eq. (2.22) ----
    // E(m,n) = 2^-m sum_{k=0}^{m} C(m,k) S_{n+k} with S_t = sum_{i<=t} (-1)^i a_i,
    // so a_i carries the mass of every partial sum that contains it.
    static double[] eulerWeights(int n, int mm) {
        double[] b = new double[mm + 1];
        b[0] = 1.0;
        for (int k = 1; k <= mm; k++) {
            b[k] = b[k - 1] * (mm - k + 1) / k;      // C(mm,k)
        }
        double scale = Math.pow(2.0, -mm);
        for (int k = 0; k <= mm; k++) {
            b[k] *= scale;
        }
        double[] tail = new double[mm + 1];          // tail[k] = 2^-mm sum_{j>=k} C(mm,j)
        double acc = 0.0;
        for (int k = mm; k >= 0; k--) {
            acc += b[k];
            tail[k] = acc;
        }
        double[] w = new double[n + mm + 1];
        for (int i = 0; i <= n; i++) {
            w[i] = 1.0;
        }
        for (int i = n + 1; i <= n + mm; i++) {
            w[i] = tail[i - n];
        }
        return w;
    }

    // ---- interdependence graph and the minimizing subset D (eqs. 3.1-3.3) ----
    private static final class Decomposition {
        final int[] D;
        final List<int[]> comps;
        Decomposition(int[] D, List<int[]> comps) {
            this.D = D;
            this.comps = comps;
        }
    }

    private static Decomposition dimred(double[][] Lv, int qd, int p, Options opt) {
        List<int[]> comps = new ArrayList<int[]>();
        if (!opt.dimred || p <= 2) {
            comps.add(seq(p));
            return new Decomposition(new int[0], comps);
        }
        boolean[][] adj = new boolean[p][p];
        for (int i = 0; i < qd; i++) {
            for (int a = 0; a < p; a++) {
                if (Lv[i][a] == 0) continue;
                for (int b = 0; b < p; b++) {
                    if (b != a && Lv[i][b] != 0) {
                        adj[a][b] = true;                // each factor is a clique (eq. 3.1)
                    }
                }
            }
        }
        boolean[] none = new boolean[p];
        List<int[]> bestC = components(adj, none, p);
        int best = maxSize(bestC);                       // |D| = 0
        int[] bestD = new int[0];
        int maxd = Math.min(opt.dimredMaxD, p - 1);
        for (int dd = 1; dd <= maxd; dd++) {
            if (dd >= best) {
                break;                                   // dimension is at least |D|
            }
            if (binom(p, dd) > 2e5) {
                break;                                   // (3.3) is solved by enumeration only
            }
            int[] sub = new int[dd];
            for (int t = 0; t < dd; t++) {
                sub[t] = t;
            }
            while (true) {
                boolean[] mask = new boolean[p];
                for (int t = 0; t < dd; t++) {
                    mask[sub[t]] = true;
                }
                List<int[]> cc = components(adj, mask, p);
                int mx = maxSize(cc);
                if (dd + mx < best) {
                    best = dd + mx;
                    bestD = Arrays.copyOf(sub, dd);
                    bestC = cc;
                }
                int t = dd - 1;
                while (t >= 0 && sub[t] == p - dd + t) {
                    t--;
                }
                if (t < 0) {
                    break;
                }
                sub[t]++;
                for (int u = t + 1; u < dd; u++) {
                    sub[u] = sub[u - 1] + 1;
                }
            }
        }
        if (best < p) {
            comps.addAll(bestC);
            return new Decomposition(bestD, comps);
        }
        comps.add(seq(p));
        return new Decomposition(new int[0], comps);
    }

    private static double binom(int n, int k) {
        double v = 1.0;
        for (int i = 1; i <= k; i++) {
            v = v * (n - k + i) / i;
        }
        return v;
    }

    private static int[] seq(int p) {
        int[] v = new int[p];
        for (int j = 0; j < p; j++) {
            v[j] = j;
        }
        return v;
    }

    private static int maxSize(List<int[]> comps) {
        int mx = 0;
        for (int c = 0; c < comps.size(); c++) {
            mx = Math.max(mx, comps.get(c).length);
        }
        return mx;
    }

    // ---- connected components of the graph with the nodes in mask removed ----
    private static List<int[]> components(boolean[][] adj, boolean[] mask, int p) {
        int[] lab = new int[p];
        Arrays.fill(lab, -1);
        List<int[]> out = new ArrayList<int[]>();
        int nc = 0;
        for (int s = 0; s < p; s++) {
            if (mask[s] || lab[s] >= 0) {
                continue;
            }
            lab[s] = nc;
            int[] stack = new int[p];
            int top = 0;
            stack[top++] = s;
            while (top > 0) {
                int v = stack[--top];
                for (int u = 0; u < p; u++) {
                    if (adj[v][u] && !mask[u] && lab[u] < 0) {
                        lab[u] = nc;
                        stack[top++] = u;
                    }
                }
            }
            nc++;
        }
        for (int c = 0; c < nc; c++) {
            int n = 0;
            for (int j = 0; j < p; j++) {
                if (lab[j] == c) n++;
            }
            int[] v = new int[n];
            int t = 0;
            for (int j = 0; j < p; j++) {
                if (lab[j] == c) v[t++] = j;
            }
            out.add(v);
        }
        return out;
    }

    // ---- scaled generating function Gbar restricted to a group of factors ----
    // Gbar_S(w) = exp( sum_{j in chains} alpha_j rho_{j0} (w_j - 1) )
    //             / prod_{i in queues} (1 - sum_j alpha_j rho_{ji} w_j)^{m_i}
    private static double[] gbarSub(Ctx ctx, int[] queues, int[] chains, double off,
                                    double[] wRe, double[] wIm) {
        double expoRe = 0.0, expoIm = 0.0;
        for (int t = 0; t < chains.length; t++) {
            int j = chains[t];
            expoRe += ctx.arho0[j] * (wRe[j] - 1.0);
            expoIm += ctx.arho0[j] * wIm[j];
        }
        double logdenRe = 0.0, logdenIm = 0.0;
        for (int t = 0; t < queues.length; t++) {
            int i = queues[t];
            double aRe = 0.0, aIm = 0.0;
            for (int j = 0; j < ctx.p; j++) {
                aRe += ctx.rhoS[i][j] * wRe[j];
                aIm += ctx.rhoS[i][j] * wIm[j];
            }
            double oRe = 1.0 - aRe;
            double oIm = -aIm;
            double lr = Math.log(Math.hypot(oRe, oIm));
            double li = Math.atan2(oIm, oRe);
            logdenRe += ctx.m[i] * lr;
            logdenIm += ctx.m[i] * li;
        }
        double eRe = expoRe - logdenRe - off;
        double eIm = expoIm - logdenIm;
        double ex = Math.exp(eRe);
        return new double[]{ex * Math.cos(eIm), ex * Math.sin(eIm)};
    }
}
