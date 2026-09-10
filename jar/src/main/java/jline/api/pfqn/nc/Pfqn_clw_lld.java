/**
 * Choudhury-Leung-Whitt normalization constant by numerical inversion of the
 * generating function (JACM 42(5):935-970, 1995), extended to limited
 * load-dependent (LLD) stations via the per-center transforms of Bertozzi and
 * McKenna (SIAM Review 35(2):239-268, 1993).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.matrix.Matrix;

import java.util.Arrays;

/**
 * Computes g(K) of a multichain closed product-form network with limited
 * load-dependent (LLD) stations and (optionally) infinite-server delay by
 * numerically inverting its p-dimensional generating function
 * (Bertozzi-McKenna eqs. 2.17/2.23)
 *
 *   G(z) = exp(sum_j rho_{j0} z_j) prod_i F_i(sum_j rho_{ji} z_j)
 *
 * where F_i is the transform of the station factor of queue i (eq. 2.16) with
 * load-dependent rate scalings S_i(k) = mu(i,k). For an LLD queue, S_i(k) =
 * c_i constant for k &gt;= l_i, and F_i is the rational function (eq. 2.19)
 *
 *   F_i(x) = [c_i + sum_{n=1}^{l_i-1} (c_i - S_i(n)) / prod_{k=1}^n S_i(k) x^n]
 *            / (c_i - x),
 *
 * analytic except for a simple pole at x = c_i. Multiserver and
 * load-independent queues are special cases. Since g(K) depends on S_i(k)
 * only for k &lt;= sum(K), general load-dependent input is truncated to LLD at
 * sum(K) without loss of exactness.
 *
 * g(K) is recovered by p nested one-dimensional lattice-Poisson inversions
 * (CLW eq. 2.3) with restrictive static scaling adapted from CLW eqs.
 * 5.41-5.46 (each queue normalized by its pole c_i, simple pole) and
 * log-domain recovery (eq. 7.1). Cost is prod_j 2 l_j K_j contour points,
 * each of cost O(sum_i l_i); practical for moderate populations and few
 * chains.
 */
public final class Pfqn_clw_lld {
    private Pfqn_clw_lld() {}

    public static Ret.pfqnNc pfqn_clw_lld(Matrix L, Matrix N) {
        return pfqn_clw_lld(L, N, null, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw_lld(Matrix L, Matrix N, Matrix Z) {
        return pfqn_clw_lld(L, N, Z, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw_lld(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        return pfqn_clw_lld(L, N, Z, mu, null, null);
    }

    /**
     * @param L      (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
     * @param N      (1 x p or p x 1) closed-chain population vector K.
     * @param Z      (1 x p) aggregate infinite-server relative intensities rho_{j0}; null = 0.
     * @param mu     (q' x n) load-dependent rate scalings mu(i,k) = S_i(k); if fewer
     *               than sum(K) columns are given the last column is extended (LLD
     *               assumption); null = ones (all queues load-independent).
     * @param lpar   (1 x p) inner lattice parameters l_j; null = 1,2,2,3,3,...
     * @param gampar (1 x p) aliasing parameters gamma_j; null = 11,13,13,15,15,...
     * @return normalization constant g(K) (Inf on overflow) and its natural log.
     */
    public static Ret.pfqnNc pfqn_clw_lld(Matrix L, Matrix N, Matrix Z, Matrix mu,
                                          Matrix lpar, Matrix gampar) {
        int qd = L.getNumRows();
        int pFull = L.getNumCols();

        int[] NvFull = new int[pFull];
        for (int j = 0; j < pFull; j++) {
            NvFull[j] = (int) Math.round(N.get(j));
        }
        double[] ZvFull = new double[pFull];
        if (Z != null && !Z.isEmpty()) {
            for (int j = 0; j < pFull; j++) {
                ZvFull[j] = Z.get(j);
            }
        }
        int[] lvFull = new int[pFull];
        if (lpar != null && !lpar.isEmpty()) {
            for (int j = 0; j < pFull; j++) {
                lvFull[j] = (int) Math.round(lpar.get(j));
            }
        } else {
            for (int j = 0; j < pFull; j++) {
                lvFull[j] = 3;
            }
            lvFull[0] = 1;
            if (pFull >= 2) lvFull[1] = 2;
            if (pFull >= 3) lvFull[2] = 2;
        }
        double[] gamvFull = new double[pFull];
        if (gampar != null && !gampar.isEmpty()) {
            for (int j = 0; j < pFull; j++) {
                gamvFull[j] = gampar.get(j);
            }
        } else {
            for (int j = 0; j < pFull; j++) {
                gamvFull[j] = 15.0;
            }
            gamvFull[0] = 11.0;
            if (pFull >= 2) gamvFull[1] = 13.0;
            if (pFull >= 3) gamvFull[2] = 13.0;
        }

        boolean anyNeg = false;
        int sumN = 0;
        for (int j = 0; j < pFull; j++) {
            if (NvFull[j] < 0) anyNeg = true;
            sumN += NvFull[j];
        }
        if (anyNeg) {
            return new Ret.pfqnNc(0.0, Double.NEGATIVE_INFINITY);
        }
        if (sumN == 0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        // load-dependent scalings extended/truncated to sum(N) columns
        double[][] muv = new double[qd][sumN];
        if (mu != null && !mu.isEmpty()) {
            int muCols = mu.getNumCols();
            for (int i = 0; i < qd; i++) {
                for (int k = 0; k < sumN; k++) {
                    muv[i][k] = mu.get(i, Math.min(k, muCols - 1));
                    if (muv[i][k] <= 0) {
                        throw new RuntimeException("pfqn_clw_lld: load-dependent rates mu(i,k) must be positive.");
                    }
                }
            }
        } else {
            for (int i = 0; i < qd; i++) {
                Arrays.fill(muv[i], 1.0);
            }
        }

        // drop zero-population chains: the coefficient of z_j^0 equals the pgf
        // restricted to z_j = 0, so chain j is removed exactly
        int p = 0;
        for (int j = 0; j < pFull; j++) {
            if (NvFull[j] > 0) p++;
        }
        int[] Nv = new int[p];
        double[] Zv = new double[p];
        int[] lv = new int[p];
        double[] gamv = new double[p];
        double[][] Lv = new double[qd][p];
        int c = 0;
        for (int j = 0; j < pFull; j++) {
            if (NvFull[j] > 0) {
                Nv[c] = NvFull[j];
                Zv[c] = ZvFull[j];
                lv[c] = lvFull[j];
                gamv[c] = gamvFull[j];
                for (int i = 0; i < qd; i++) {
                    Lv[i][c] = L.get(i, j);
                }
                c++;
            }
        }

        // pole c_i and LLD cutoff l_i of each queue: S_i(k) = c_i for k >= l_i
        double[] cpole = new double[qd];
        double[][] numc = new double[qd][];   // numerator coefficients a_0..a_{l_i-1}
        for (int i = 0; i < qd; i++) {
            cpole[i] = muv[i][sumN - 1];
            int last = 0;
            for (int k = sumN - 1; k >= 1; k--) {
                if (muv[i][k - 1] != cpole[i]) {
                    last = k;
                    break;
                }
            }
            int li = last + 1;
            double[] a = new double[li];
            a[0] = cpole[i];
            double cp = 1.0;
            for (int n = 1; n < li; n++) {
                cp *= muv[i][n - 1];               // prod_{k=1}^n S_i(k)
                a[n] = (cpole[i] - muv[i][n - 1]) / cp;
            }
            numc[i] = a;
        }

        // contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (CLW eq. 2.7)
        double[] r = new double[p];
        for (int j = 0; j < p; j++) {
            r[j] = Math.pow(10.0, -gamv[j] / (2.0 * lv[j] * Nv[j]));
        }

        // see _kb/03-api-layer.md for rationale
        double[][] Lt = new double[qd][p];
        for (int i = 0; i < qd; i++) {
            for (int j = 0; j < p; j++) {
                Lt[i][j] = Lv[i][j] / cpole[i];
            }
        }
        double[] alpha = new double[p];
        double[] used = new double[qd];
        for (int j = 0; j < p; j++) {
            int Kj = Nv[j];
            int lj = lv[j];
            int nc = 0;
            for (int i = 0; i < qd; i++) {
                if (Lt[i][j] > 0) nc++;
            }
            double aj = Double.POSITIVE_INFINITY;
            if (nc > 0) {
                Integer[] idx = new Integer[nc];
                double[] eAll = new double[qd];
                int cc = 0;
                for (int i = 0; i < qd; i++) {
                    double denom = 1.0 - used[i];
                    if (denom <= 0) denom = Double.MIN_VALUE;
                    eAll[i] = Lt[i][j] / denom;
                    if (Lt[i][j] > 0) {
                        idx[cc++] = i;
                    }
                }
                final double[] eRef = eAll;
                Arrays.sort(idx, (x, y) -> Double.compare(eRef[y], eRef[x]));  // descending
                double sumE = 0.0;
                for (int n = 0; n < nc; n++) {
                    int qi = idx[n];
                    sumE += eAll[qi];
                    double cumrho = sumE / (n + 1);
                    // N_{ij} = n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43, m_i = 1)
                    double Nnd = n;
                    for (int k = j + 1; k < p; k++) {
                        if (Lv[qi][k] != 0) Nnd += Nv[k];
                    }
                    int Nn = (int) Math.round(Nnd);
                    double an;
                    if (Nn <= 0) {
                        an = 1.0;
                    } else {
                        // in the log domain: the product runs over N_{ij}
                        // factors below one and underflows to zero at a few
                        // hundred of them, which would silently set alpha_j = 0
                        // and lG = NaN
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
            alpha[j] = aj;
            for (int i = 0; i < qd; i++) {
                used[i] += aj * Lt[i][j] * r[j];
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

        Ctx ctx = new Ctx(qd, p, Nv, lv, r, arho0, rhoS, cpole, numc);
        double[] wRe = new double[p];
        double[] wIm = new double[p];
        double[] gbar = invert(0, ctx, wRe, wIm);   // real part carries gbar(K)

        double lsum = 0.0;
        for (int j = 0; j < p; j++) {
            lsum += arho0[j] - Nv[j] * Math.log(alpha[j]);
        }
        double lG = Math.log(gbar[0]) + lsum;
        double G = (lG > 709) ? Double.POSITIVE_INFINITY : Math.exp(lG);
        return new Ret.pfqnNc(G, lG);
    }

    // immutable context shared (read-only) across the recursion
    private static final class Ctx {
        final int qd, p;
        final int[] N, l;
        final double[] r, arho0, cpole;
        final double[][] rhoS, numc;
        Ctx(int qd, int p, int[] N, int[] l, double[] r, double[] arho0,
            double[][] rhoS, double[] cpole, double[][] numc) {
            this.qd = qd; this.p = p; this.N = N; this.l = l;
            this.r = r; this.arho0 = arho0; this.rhoS = rhoS;
            this.cpole = cpole; this.numc = numc;
        }
    }

    // one-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled; returns {re,im}.
    // wRe/wIm hold the outer contour coordinates fixed at indices < j.
    private static double[] invert(int j, Ctx ctx, double[] wRe, double[] wIm) {
        int Kj = ctx.N[j];
        int lj = ctx.l[j];
        double rj = ctx.r[j];
        double accRe = 0.0, accIm = 0.0;
        for (int k1 = 0; k1 < lj; k1++) {
            double phRe = Math.cos(-Math.PI * k1 / lj);
            double phIm = Math.sin(-Math.PI * k1 / lj);
            double innerRe = 0.0, innerIm = 0.0;
            for (int k = -Kj; k < Kj; k++) {
                double sign = (k % 2 == 0) ? 1.0 : -1.0;
                double theta = Math.PI * (k1 + (double) lj * k) / ((double) lj * Kj);
                wRe[j] = rj * Math.cos(theta);
                wIm[j] = rj * Math.sin(theta);
                double vRe, vIm;
                if (j == ctx.p - 1) {
                    double[] g = gbarEval(ctx, wRe, wIm);
                    vRe = g[0];
                    vIm = g[1];
                } else {
                    double[] sub = invert(j + 1, ctx, wRe, wIm);
                    vRe = sub[0];
                    vIm = sub[1];
                }
                innerRe += sign * vRe;
                innerIm += sign * vIm;
            }
            accRe += phRe * innerRe - phIm * innerIm;
            accIm += phRe * innerIm + phIm * innerRe;
        }
        double denom = 2.0 * lj * Kj * Math.pow(rj, Kj);
        double valRe = accRe / denom;
        double valIm = accIm / denom;
        if (j == 0) {
            valIm = 0.0;
        }
        return new double[]{valRe, valIm};
    }

    // scaled generating function Gbar(w), returns {re,im}:
    // Gbar(w) = exp(sum_j arho0_j (w_j-1)) prod_i F_i(sum_j rhoS_ij w_j) with
    // F_i(x) = N_i(x)/(c_i - x) (Bertozzi-McKenna eq. 2.19); F_i(0) = 1.
    // exp(log a + log b) = a*b for the principal complex log, so branch
    // choices in the per-queue logs are immaterial.
    private static double[] gbarEval(Ctx ctx, double[] wRe, double[] wIm) {
        double eRe = 0.0, eIm = 0.0;
        for (int j = 0; j < ctx.p; j++) {
            eRe += ctx.arho0[j] * (wRe[j] - 1.0);
            eIm += ctx.arho0[j] * wIm[j];
        }
        for (int i = 0; i < ctx.qd; i++) {
            double xRe = 0.0, xIm = 0.0;
            for (int j = 0; j < ctx.p; j++) {
                xRe += ctx.rhoS[i][j] * wRe[j];
                xIm += ctx.rhoS[i][j] * wIm[j];
            }
            double[] a = ctx.numc[i];
            double numRe = a[a.length - 1], numIm = 0.0;   // Horner on N_i(x)
            for (int k = a.length - 2; k >= 0; k--) {
                double tRe = numRe * xRe - numIm * xIm + a[k];
                double tIm = numRe * xIm + numIm * xRe;
                numRe = tRe;
                numIm = tIm;
            }
            double dRe = ctx.cpole[i] - xRe;
            double dIm = -xIm;
            eRe += Math.log(Math.hypot(numRe, numIm)) - Math.log(Math.hypot(dRe, dIm));
            eIm += Math.atan2(numIm, numRe) - Math.atan2(dIm, dRe);
        }
        double ex = Math.exp(eRe);
        return new double[]{ex * Math.cos(eIm), ex * Math.sin(eIm)};
    }
}
