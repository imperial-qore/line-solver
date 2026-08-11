/**
 * Choudhury-Leung-Whitt normalization constant by numerical inversion of the
 * generating function (JACM 42(5):935-970, 1995).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.Arrays;

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
 * multiplicity m_i. g(K) is recovered by p nested one-dimensional lattice-Poisson
 * inversions (eq. 2.3) with restrictive static scaling (eqs. 5.41-5.46) and
 * log-domain recovery (eq. 7.1).
 *
 * Exact nested inversion of cost prod_j 2 l_j K_j; practical for moderate
 * populations and few chains. The paper's Euler summation and dimension
 * reduction speed-ups are not applied here.
 */
public final class Pfqn_clw {
    private Pfqn_clw() {}

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N) {
        return pfqn_clw(L, N, null, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z) {
        return pfqn_clw(L, N, Z, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z, Matrix m) {
        return pfqn_clw(L, N, Z, m, null, null);
    }

    /**
     * @param L     (q' x p) single-server relative traffic intensities, L(i,j)=rho_{ji}.
     * @param N     (1 x p or p x 1) closed-chain population vector K.
     * @param Z     (1 x p) aggregate infinite-server relative intensities rho_{j0}; null = 0.
     * @param m     (q' x 1) queue multiplicities m_i; null = ones.
     * @param lpar  (1 x p) inner lattice parameters l_j; null = 1,2,2,3,3,...
     * @param gampar(1 x p) aliasing parameters gamma_j; null = 11,13,13,15,15,...
     * @return normalization constant g(K) (Inf on overflow) and its natural log.
     */
    public static Ret.pfqnNc pfqn_clw(Matrix L, Matrix N, Matrix Z, Matrix m, Matrix lpar, Matrix gampar) {
        int qd = L.getNumRows();
        int p = L.getNumCols();

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
        int[] lv = new int[p];
        if (lpar != null && !lpar.isEmpty()) {
            for (int j = 0; j < p; j++) {
                lv[j] = (int) Math.round(lpar.get(j));
            }
        } else {
            for (int j = 0; j < p; j++) {
                lv[j] = 3;
            }
            lv[0] = 1;
            if (p >= 2) lv[1] = 2;
            if (p >= 3) lv[2] = 2;
        }
        double[] gamv = new double[p];
        if (gampar != null && !gampar.isEmpty()) {
            for (int j = 0; j < p; j++) {
                gamv[j] = gampar.get(j);
            }
        } else {
            for (int j = 0; j < p; j++) {
                gamv[j] = 15.0;
            }
            gamv[0] = 11.0;
            if (p >= 2) gamv[1] = 13.0;
            if (p >= 3) gamv[2] = 13.0;
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

        // contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (eq. 2.7)
        double[] r = new double[p];
        for (int j = 0; j < p; j++) {
            r[j] = (Nv[j] == 0) ? 1.0
                    : Math.pow(10.0, -gamv[j] / (2.0 * lv[j] * Nv[j]));
        }

        // restrictive static scaling (eqs. 5.41-5.46), outer vars at |z_k| = r_k
        double[] alpha = new double[p];
        double[] used = new double[qd];
        for (int j = 0; j < p; j++) {
            int Kj = Nv[j];
            int lj = lv[j];
            // effective intensities e_i = L_ij / (1 - used_i), positive queues only
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
                    // N_{ij} = mbar_n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43)
                    double Nnd = cummb - 1.0;
                    for (int k = j + 1; k < p; k++) {
                        if (Lv[qi][k] != 0) Nnd += Nv[k];
                    }
                    int Nn = (int) Math.round(Nnd);
                    double an;
                    if (Nn <= 0) {
                        an = 1.0;
                    } else {
                        double prod = 1.0;
                        double twolK = 2.0 * lj * Kj;
                        for (int ll = 1; ll <= Nn; ll++) {
                            prod *= (Kj + ll) / (Kj + twolK + ll);
                        }
                        an = Math.pow(prod, 1.0 / twolK);
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
                used[i] += aj * Lv[i][j] * r[j];
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

        Ctx ctx = new Ctx(qd, p, Nv, lv, r, arho0, rhoS, mv);
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
        final double[] r, arho0;
        final double[][] rhoS;
        final double[] m;
        Ctx(int qd, int p, int[] N, int[] l, double[] r, double[] arho0, double[][] rhoS, double[] m) {
            this.qd = qd; this.p = p; this.N = N; this.l = l;
            this.r = r; this.arho0 = arho0; this.rhoS = rhoS; this.m = m;
        }
    }

    // one-dimensional lattice-Poisson inversion (eq. 2.3), scaled; returns {re,im}.
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

    // scaled generating function Gbar(w), returns {re,im}
    private static double[] gbarEval(Ctx ctx, double[] wRe, double[] wIm) {
        double expoRe = 0.0, expoIm = 0.0;
        for (int j = 0; j < ctx.p; j++) {
            expoRe += ctx.arho0[j] * (wRe[j] - 1.0);
            expoIm += ctx.arho0[j] * wIm[j];
        }
        double logdenRe = 0.0, logdenIm = 0.0;
        for (int i = 0; i < ctx.qd; i++) {
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
        double eRe = expoRe - logdenRe;
        double eIm = expoIm - logdenIm;
        double ex = Math.exp(eRe);
        return new double[]{ex * Math.cos(eIm), ex * Math.sin(eIm)};
    }
}
