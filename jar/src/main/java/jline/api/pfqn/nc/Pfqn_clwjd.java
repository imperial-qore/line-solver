/**
 * Limited joint-dependent (LJD) normalization constant by numerical inversion of
 * the multichain generating function (Choudhury-Leung-Whitt, J. ACM 42(5):935-970,
 * 1995), the joint-dependent generalization of {@link Pfqn_clwoi}.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.ToDoubleFunction;

/**
 * Computes G(N) of a closed product-form network made of an aggregated
 * infinite-server (delay) node and any number of LIMITED JOINT-DEPENDENT
 * stations. This stands to {@link Pfqn_clwoi} as {@link Pfqn_clw_lld} stands to
 * {@link Pfqn_clw}: a per-station cutoff beyond which the rate stops changing
 * turns an infinite series into a rational function of the same denominators.
 *
 * <p>Station i has a rate mu_i(n) that reads the whole per-class occupancy but
 * saturates coordinatewise: with a cutoff vector l_i,
 * mu_i(n) = c_{i,t}, t = (min(n_1,l_{i,1}), ..., min(n_R,l_{i,R})), so past
 * l_{i,r} further class-r jobs no longer change the rate. Order independence is
 * l_i = 1 (t is the support indicator); a multiserver station with c servers is
 * l_i = c, min(sum n, c) being a function of the clipped vector once every
 * l_{i,r} &gt;= c.
 *
 * <p>Splitting the count lattice by clipped region, on which mu_i is constant,
 * and writing F_{i,t} for the part of F_i carried by the states with t_i(n) = t,
 *
 *   (mu_{i,t} - sum_{r: t_r = l_{i,r}} v_{i,r} z_r) F_{i,t}(z)
 *        = sum_{r: t_r &gt;= 1} v_{i,r} z_r F_{i,t-e_r}(z),   F_{i,0} = 1,
 *   F_i(z) = sum_t F_{i,t}(z),
 *
 * the two sides differing because removing a class-r job leaves the region only
 * on an UNsaturated coordinate: for t_r &lt; l_{i,r} the region pins n_r = t_r so
 * n - e_r lands in t - e_r, while for t_r = l_{i,r} the region is n_r &gt;= l and
 * n - e_r lands in t or in t - e_r. The singularities are therefore the
 * hyperplanes sum_{r in S} v_{i,r} z_r = mu_{i,t} over the SATURATED sets
 * S = {r : t_r = l_{i,r}}: at most 2^R per station, however large the cutoffs.
 *
 * <p>The restrictive static scaling of CLW eqs. 5.41-5.46 runs on one row per
 * (station, saturated set), carrying v_{i,r}/min{mu_{i,t} : saturated set of t
 * is S}, the smallest rate over regions sharing a saturated set being the
 * binding one; rows dominated by a superset of no larger rate are dropped.
 *
 * <p>The rate must be constant on each clipped region, which is checked on probe
 * states. Any rate is admissible with lcut = N (the default), the clipping being
 * vacuous on the reachable lattice. Cost is prod_r 2 l_r N_r contour points, each
 * costing O(M R prod_r (l_{i,r}+1)), against O(M prod_r (N_r+1)(N_r+2)/2) for
 * {@link Pfqn_ncjd}: the inversion pays off when the joint dependence saturates
 * early, and loses outright at lcut = N.
 */
public final class Pfqn_clwjd {
    private Pfqn_clwjd() {}

    public static Ret.pfqnNc pfqn_clwjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        return pfqn_clwjd(Z, N, mu, null, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clwjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                        double[][] visits) {
        return pfqn_clwjd(Z, N, mu, visits, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clwjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                        double[][] visits, int[][] lcut) {
        return pfqn_clwjd(Z, N, mu, visits, lcut, null, null);
    }

    /**
     * @param Z      (R) think-time demand of the aggregated delay node.
     * @param N      (R) closed population, finite.
     * @param mu     one rate handle per joint-dependent station, mapping a per-class
     *               count vector to the total service rate. May be null/empty for a
     *               pure delay network.
     * @param visits (M x R) per-station class visit ratios weighting the balance
     *               recursion; null = unit visits.
     * @param lcut   (M x R) per-station per-class saturation cutoffs l_{i,r} &gt;= 1;
     *               entries are clipped to N_r, which is exact because a rate
     *               difference at n_r &gt; N_r can only move coefficients with
     *               n_r &gt; N_r. null = N (no truncation).
     * @param lpar   (R) inner lattice parameters l_j; null = 1,2,2,3,3,...
     * @param gampar (R) aliasing parameters gamma_j; null = 11,13,13,15,15,...
     * @return G(N) (Inf on overflow) and its natural logarithm.
     */
    public static Ret.pfqnNc pfqn_clwjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                        double[][] visits, int[][] lcut, int[] lpar,
                                        double[] gampar) {
        int R = N.length;
        List<ToDoubleFunction<int[]>> muList = (mu == null) ? new ArrayList<ToDoubleFunction<int[]>>() : mu;
        int M = muList.size();

        double[] Zfull = (Z == null || Z.length == 0) ? new double[R] : Z;
        if (Zfull.length != R) {
            throw new RuntimeException("pfqn_clwjd: Z and N must have the same number of classes.");
        }

        boolean anyNeg = false;
        int sumN = 0;
        for (int r = 0; r < R; r++) {
            if (N[r] < 0) anyNeg = true;
            sumN += Math.max(N[r], 0);
        }
        if (anyNeg) {
            return new Ret.pfqnNc(0.0, Double.NEGATIVE_INFINITY);
        }
        if (sumN == 0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        int[] lvFull = new int[R];
        if (lpar != null && lpar.length == R) {
            System.arraycopy(lpar, 0, lvFull, 0, R);
        } else {
            for (int r = 0; r < R; r++) lvFull[r] = 3;
            lvFull[0] = 1;
            if (R >= 2) lvFull[1] = 2;
            if (R >= 3) lvFull[2] = 2;
        }
        double[] gamvFull = new double[R];
        if (gampar != null && gampar.length == R) {
            System.arraycopy(gampar, 0, gamvFull, 0, R);
        } else {
            for (int r = 0; r < R; r++) gamvFull[r] = 15.0;
            gamvFull[0] = 11.0;
            if (R >= 2) gamvFull[1] = 13.0;
            if (R >= 3) gamvFull[2] = 13.0;
        }

        // drop zero-population chains: the coefficient of z_r^0 is the generating
        // function at z_r = 0, which kills every F_{i,t} with t_r >= 1
        int Rk = 0;
        for (int r = 0; r < R; r++) {
            if (N[r] > 0) Rk++;
        }
        int[] keep = new int[Rk];
        int[] Nk = new int[Rk];
        double[] Zk = new double[Rk];
        int[] lv = new int[Rk];
        double[] gamv = new double[Rk];
        int c = 0;
        for (int r = 0; r < R; r++) {
            if (N[r] > 0) {
                keep[c] = r;
                Nk[c] = N[r];
                Zk[c] = Zfull[r];
                lv[c] = lvFull[r];
                gamv[c] = gamvFull[r];
                c++;
            }
        }

        // saturation cutoffs, broadcast and clipped to the reachable lattice
        int[][] Lk = new int[Math.max(M, 1)][Rk];
        for (int i = 0; i < M; i++) {
            for (int b = 0; b < Rk; b++) {
                int li = (lcut == null) ? Nk[b] : lcut[i][keep[b]];
                if (li < 1) li = 1;
                if (li > Nk[b]) li = Nk[b];
                Lk[i][b] = li;
            }
        }

        // per-station region tables over the clipped box prod_r {0,...,l_{i,r}},
        // in mixed radix so that t - e_r always precedes t
        int[] ntreg = new int[Math.max(M, 1)];
        double[][] muT = new double[Math.max(M, 1)][];
        int[][][] regDec = new int[Math.max(M, 1)][][];   // [i][tl] = flattened (chain, subcol) pairs
        int[][][] regSat = new int[Math.max(M, 1)][][];   // [i][tl] = saturated chains
        int[][] satMask = new int[Math.max(M, 1)][];
        for (int i = 0; i < M; i++) {
            int[] st = new int[Rk];
            int nt = 1;
            for (int b = 0; b < Rk; b++) {
                st[b] = nt;
                nt *= Lk[i][b] + 1;
            }
            ntreg[i] = nt;
            muT[i] = new double[nt];
            regDec[i] = new int[nt][];
            regSat[i] = new int[nt][];
            satMask[i] = new int[nt];
            int[] t = new int[Rk];
            for (int tl = 0; tl < nt; tl++) {
                for (int b = 0; b < Rk; b++) {
                    t[b] = (tl / st[b]) % (Lk[i][b] + 1);
                }
                int ndec = 0, nsat = 0, smask = 0;
                for (int b = 0; b < Rk; b++) {
                    if (t[b] >= 1) ndec++;
                    if (t[b] == Lk[i][b]) {
                        nsat++;
                        smask |= (1 << b);
                    }
                }
                int[] dec = new int[2 * ndec];
                int[] sat = new int[nsat];
                int kd = 0, ks = 0;
                for (int b = 0; b < Rk; b++) {
                    if (t[b] >= 1) {
                        dec[kd++] = b;
                        dec[kd++] = tl - st[b];
                    }
                    if (t[b] == Lk[i][b]) sat[ks++] = b;
                }
                regDec[i][tl] = dec;
                regSat[i][tl] = sat;
                satMask[i][tl] = smask;
                muT[i][tl] = (tl == 0) ? 1.0 : regionRate(muList.get(i), t, keep, R, N, Lk[i], i);
            }
        }

        // per-station visit vectors restricted to the retained chains
        double[][] V = new double[Math.max(M, 1)][Rk];
        for (int i = 0; i < M; i++) {
            for (int b = 0; b < Rk; b++) {
                V[i][b] = (visits == null) ? 1.0 : visits[i][keep[b]];
            }
        }

        // contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (CLW eq. 2.7)
        double[] rad = new double[Rk];
        for (int j = 0; j < Rk; j++) {
            rad[j] = Math.pow(10.0, -gamv[j] / (2.0 * lv[j] * Nk[j]));
        }

        // binding rate of each saturated set: regions sharing a saturated set share
        // the hyperplane sum_{r in S} v_r z_r = mu, so the smallest rate constrains
        int nmask = 1 << Rk;
        double[][] muS = new double[Math.max(M, 1)][nmask];
        for (int i = 0; i < M; i++) {
            Arrays.fill(muS[i], Double.POSITIVE_INFINITY);
            for (int tl = 0; tl < ntreg[i]; tl++) {
                int sm = satMask[i][tl];
                if (sm > 0 && muT[i][tl] < muS[i][sm]) muS[i][sm] = muT[i][tl];
            }
        }

        // restrictive static scaling (CLW eqs. 5.41-5.46), one row per (station,
        // saturated set), with rows dominated by a superset of no larger rate dropped
        List<double[]> rows = new ArrayList<double[]>();
        for (int i = 0; i < M; i++) {
            for (int mask = 1; mask < nmask; mask++) {
                if (!Double.isFinite(muS[i][mask])) continue;
                boolean dominated = false;
                for (int mask2 = 1; mask2 < nmask; mask2++) {
                    if (mask2 != mask && (mask & mask2) == mask && Double.isFinite(muS[i][mask2])
                            && muS[i][mask2] <= muS[i][mask] * (1 + 1e-12)) {
                        dominated = true;
                        break;
                    }
                }
                if (dominated) continue;
                double[] row = new double[Rk];
                for (int b = 0; b < Rk; b++) {
                    if ((mask & (1 << b)) != 0) row[b] = V[i][b] / muS[i][mask];
                }
                rows.add(row);
            }
        }
        int nrow = Math.max(rows.size(), 1);
        double[][] Lt = new double[nrow][Rk];
        for (int i = 0; i < rows.size(); i++) {
            Lt[i] = rows.get(i);
        }

        double[] alpha = new double[Rk];
        double[] used = new double[nrow];
        for (int j = 0; j < Rk; j++) {
            int Kj = Nk[j];
            int lj = lv[j];
            int nc = 0;
            for (int i = 0; i < nrow; i++) {
                if (Lt[i][j] > 0) nc++;
            }
            double aj = Double.POSITIVE_INFINITY;
            if (nc > 0) {
                Integer[] idx = new Integer[nc];
                double[] eAll = new double[nrow];
                int cc = 0;
                for (int i = 0; i < nrow; i++) {
                    double denom = 1.0 - used[i];
                    if (denom <= 0) denom = Double.MIN_VALUE;
                    eAll[i] = Lt[i][j] / denom;
                    if (Lt[i][j] > 0) idx[cc++] = i;
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
                    for (int k = j + 1; k < Rk; k++) {
                        if (Lt[qi][k] != 0) Nnd += Nk[k];
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
            if (Zk[j] > 0) {
                aj = Math.min(aj, Kj / Zk[j]);
            }
            if (!Double.isFinite(aj)) {
                aj = 1.0;
            }
            alpha[j] = aj;
            for (int i = 0; i < nrow; i++) {
                used[i] += aj * Lt[i][j] * rad[j];
            }
        }

        double[] arho0 = new double[Rk];
        double[][] vs = new double[Math.max(M, 1)][Rk];
        for (int j = 0; j < Rk; j++) {
            arho0[j] = alpha[j] * Zk[j];
            for (int i = 0; i < M; i++) {
                vs[i][j] = V[i][j] * alpha[j];
            }
        }

        Ctx ctx = new Ctx(M, Rk, Nk, lv, rad, arho0, vs, muT, regDec, regSat, ntreg);
        double[] wRe = new double[Rk];
        double[] wIm = new double[Rk];
        double[] gbar = invert(0, ctx, wRe, wIm);

        double lsum = 0.0;
        for (int j = 0; j < Rk; j++) {
            lsum += arho0[j] - Nk[j] * Math.log(alpha[j]);
        }
        double lG = Math.log(gbar[0]) + lsum;
        double G = (lG > 709) ? Double.POSITIVE_INFINITY : Math.exp(lG);
        return new Ret.pfqnNc(G, lG);
    }

    // rate of one clipped region, with the constancy check: the region pins every
    // unsaturated coordinate and leaves the saturated ones free above the cutoff,
    // so the rate is probed at the representative and above the cutoff
    private static double regionRate(ToDoubleFunction<int[]> murate, int[] t, int[] keep,
                                     int R, int[] N, int[] lrow, int i) {
        int[] nrep = new int[R];
        for (int b = 0; b < t.length; b++) nrep[keep[b]] = t[b];
        double rate = murate.applyAsDouble(nrep);
        if (!(rate > 0)) {
            throw new RuntimeException("pfqn_clwjd: station " + (i + 1)
                    + " has a non-positive rate on a reachable region.");
        }
        for (int pass = 0; pass < 2; pass++) {
            int[] probe = nrep.clone();
            boolean differs = false;
            for (int b = 0; b < t.length; b++) {
                if (t[b] != lrow[b]) continue;
                int nr = N[keep[b]];
                probe[keep[b]] = (pass == 0) ? nr : Math.max(t[b], (t[b] + nr) / 2);
                if (probe[keep[b]] != nrep[keep[b]]) differs = true;
            }
            if (!differs) continue;
            double rt = murate.applyAsDouble(probe);
            if (Math.abs(rt - rate) > 1e-9 * Math.max(1.0, Math.abs(rate))) {
                throw new RuntimeException("pfqn_clwjd: station " + (i + 1)
                        + " has a rate that varies within a clipped region (mu=" + rate
                        + " at the region representative, " + rt + " above the cutoff);"
                        + " raise lcut or use pfqn_ncjd.");
            }
        }
        return rate;
    }

    // immutable context shared (read-only) across the recursion
    private static final class Ctx {
        final int M, p;
        final int[] N, l, ntreg;
        final double[] r, arho0;
        final double[][] vs, muT;
        final int[][][] regDec, regSat;
        Ctx(int M, int p, int[] N, int[] l, double[] r, double[] arho0, double[][] vs,
            double[][] muT, int[][][] regDec, int[][][] regSat, int[] ntreg) {
            this.M = M; this.p = p; this.N = N; this.l = l; this.r = r;
            this.arho0 = arho0; this.vs = vs; this.muT = muT;
            this.regDec = regDec; this.regSat = regSat; this.ntreg = ntreg;
        }
    }

    // one-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled; returns {re,im}
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
    // Gbar(w) = exp(sum_r arho0_r (w_r-1)) prod_i F_i(alpha_r v_{i,r} w_r), with F_i
    // assembled by the clipped-region recursion. Summing logs is legitimate for the
    // principal complex log because exp(log a + log b) = a b.
    private static double[] gbarEval(Ctx ctx, double[] wRe, double[] wIm) {
        double eRe = 0.0, eIm = 0.0;
        for (int j = 0; j < ctx.p; j++) {
            eRe += ctx.arho0[j] * (wRe[j] - 1.0);
            eIm += ctx.arho0[j] * wIm[j];
        }
        for (int i = 0; i < ctx.M; i++) {
            int nt = ctx.ntreg[i];
            double[] fRe = new double[nt];
            double[] fIm = new double[nt];
            fRe[0] = 1.0;
            double sumRe = 1.0, sumIm = 0.0;
            for (int tl = 1; tl < nt; tl++) {
                double numRe = 0.0, numIm = 0.0;
                int[] dec = ctx.regDec[i][tl];
                for (int k = 0; k < dec.length; k += 2) {
                    int b = dec[k];
                    int sc = dec[k + 1];
                    double xRe = ctx.vs[i][b] * wRe[b];
                    double xIm = ctx.vs[i][b] * wIm[b];
                    numRe += xRe * fRe[sc] - xIm * fIm[sc];
                    numIm += xRe * fIm[sc] + xIm * fRe[sc];
                }
                double denRe = ctx.muT[i][tl], denIm = 0.0;
                int[] sat = ctx.regSat[i][tl];
                for (int k = 0; k < sat.length; k++) {
                    int b = sat[k];
                    denRe -= ctx.vs[i][b] * wRe[b];
                    denIm -= ctx.vs[i][b] * wIm[b];
                }
                double d2 = denRe * denRe + denIm * denIm;
                fRe[tl] = (numRe * denRe + numIm * denIm) / d2;
                fIm[tl] = (numIm * denRe - numRe * denIm) / d2;
                sumRe += fRe[tl];
                sumIm += fIm[tl];
            }
            eRe += Math.log(Math.hypot(sumRe, sumIm));
            eIm += Math.atan2(sumIm, sumRe);
        }
        double ex = Math.exp(eRe);
        return new double[]{ex * Math.cos(eIm), ex * Math.sin(eIm)};
    }
}
