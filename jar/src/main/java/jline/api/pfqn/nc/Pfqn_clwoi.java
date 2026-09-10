/**
 * Order-independent (OI) normalization constant by numerical inversion of the
 * multichain generating function (Choudhury-Leung-Whitt, J. ACM 42(5):935-970,
 * 1995), the OI counterpart of {@link Pfqn_clw_lld}.
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
 * infinite-server (delay) node and any number of order-independent (OI) /
 * pass-and-swap stations with empty swap graph, by inverting
 *
 *   G(z) = exp(sum_r Z_r z_r) prod_i F_i(z),   F_i(z) = sum_n Phi_i(n) z^n,
 *
 * with Phi_i the v-weighted balanced-fairness balance function of station i,
 * mu_i(n) Phi_i(n) = sum_{r: n_r&gt;0} v_{i,r} Phi_i(n - e_r).
 *
 * <p>Unlike a load-dependent station, whose factor collapses to a function of
 * the single argument sum_r rho_{ri} z_r, an OI station factor depends on the
 * whole vector z, because mu_i(n) depends on the occupancy through its SUPPORT
 * supp(n) = {r : n_r &gt; 0}. It is nevertheless rational and available in closed
 * form: splitting the count lattice by support, on which mu_i(n) = mu_{i,S} is
 * constant, and writing F_{i,S} for the part of F_i carried by the states of
 * support S, the balance recursion gives
 *
 *   (mu_{i,S} - sum_{r in S} v_{i,r} z_r) F_{i,S}(z)
 *        = sum_{r in S} v_{i,r} z_r F_{i,S\{r}}(z),   F_{i,{}} = 1,
 *   F_i(z) = sum_S F_{i,S}(z),
 *
 * since removing a class-r job from a state of support S lands on support S
 * when n_r &gt;= 2 and on S\{r} when n_r = 1. The singularities are the
 * hyperplanes sum_{r in S} v_{i,r} z_r = mu_{i,S}, one per support, in place of
 * the single pole x = c_i of the load-dependent case. A load-independent
 * single-server queue (mu_{i,S} = 1) gives back 1/(1 - sum_r v_{i,r} z_r).
 *
 * <p>G(N) is the coefficient of prod_r z_r^{N_r}, recovered by R nested
 * one-dimensional lattice-Poisson inversions (CLW eq. 2.3). The restrictive
 * static scaling of eqs. 5.41-5.46 is applied to the expanded constraint matrix
 * that lists one row per (station, nonempty support) pair with unit-pole
 * intensities v_{i,r}/mu_{i,S}, after dropping the rows dominated by a superset
 * of no larger rate; each surviving row is a binding singular hyperplane.
 *
 * <p>Rates must be support-only, mu_i(n) = mu_i(supp(n)); this is the defining
 * property of an OI station and what makes the transform a finite rational
 * function. Every rate handle is verified EXHAUSTIVELY on the count lattice
 * 0 &lt; n &lt;= N before the inversion, at prod_r (N_r+1) evaluations per station
 * (below the contour points spent afterwards), and a state whose rate differs
 * from that of its support is rejected with that state named: a rate that
 * varies inside a support is a general balanced-fairness station and belongs to
 * {@link Pfqn_ncoi}. It is refused rather than warned-and-inverted because the
 * inversion would otherwise return a plausible but wrong G(N). Cost is prod_r 2 l_r N_r contour points, each costing
 * O(M R 2^R), against O(M prod_r (N_r+1)(N_r+2)/2) for the convolution of
 * {@link Pfqn_ncoi}: linear rather than quadratic in each population, and so
 * preferable on large populations with few chains.
 */
public final class Pfqn_clwoi {
    private Pfqn_clwoi() {}

    public static Ret.pfqnNc pfqn_clwoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        return pfqn_clwoi(Z, N, mu, null, null, null);
    }

    public static Ret.pfqnNc pfqn_clwoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                         double[][] visits) {
        return pfqn_clwoi(Z, N, mu, visits, null, null);
    }

    /**
     * @param Z      (R) think-time demand of the aggregated delay node.
     * @param N      (R) closed population, finite.
     * @param mu     list of OI station rate handles; mu.get(i) maps a per-class count
     *               vector n (length R) to the total service rate of station i, and
     *               must depend on n only through its support. May be null/empty for a
     *               pure delay network.
     * @param visits (M x R) per-station class visit ratios weighting the balance
     *               recursion; null = unit visits.
     * @param lpar   (R) inner lattice parameters l_j; null = 1,2,2,3,3,...
     * @param gampar (R) aliasing parameters gamma_j; null = 11,13,13,15,15,...
     * @return G(N) (Inf on overflow) and its natural logarithm.
     */
    public static Ret.pfqnNc pfqn_clwoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                         double[][] visits, int[] lpar, double[] gampar) {
        int R = N.length;
        List<ToDoubleFunction<int[]>> muList = (mu == null) ? new ArrayList<ToDoubleFunction<int[]>>() : mu;
        int M = muList.size();

        double[] Zfull = (Z == null || Z.length == 0) ? new double[R] : Z;
        if (Zfull.length != R) {
            throw new RuntimeException("pfqn_clwoi: Z and N must have the same number of classes.");
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
        // function at z_r = 0, which kills every F_{i,S} with r in S
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
        int nmask = 1 << Rk;

        // support rate table mu_{i,S}, S a bitmask over the retained chains
        double[][] muS = new double[Math.max(M, 1)][nmask];
        for (int i = 0; i < M; i++) {
            for (int mask = 1; mask < nmask; mask++) {
                int[] chi = new int[R];
                for (int b = 0; b < Rk; b++) {
                    if ((mask & (1 << b)) != 0) chi[keep[b]] = 1;
                }
                muS[i][mask] = supportRate(muList.get(i), chi, i);
            }
        }

        // exhaustive support-only verification: the transform is exact only if
        // mu_i is constant on each support, and a rate that violates that returns
        // a wrong G with no other symptom, so it is refused rather than inverted
        checkSupportOnly(muList, muS, keep, Nk, Rk, R);

        // per-mask chain lists and the S\{r} column indices of the recursion
        int[][] bits = new int[nmask][];
        int[][] subcol = new int[nmask][];
        for (int mask = 1; mask < nmask; mask++) {
            int nb = Integer.bitCount(mask);
            bits[mask] = new int[nb];
            subcol[mask] = new int[nb];
            int t = 0;
            for (int b = 0; b < Rk; b++) {
                if ((mask & (1 << b)) != 0) {
                    bits[mask][t] = b;
                    subcol[mask][t] = mask & ~(1 << b);
                    t++;
                }
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

        // Restrictive static scaling (CLW eqs. 5.41-5.46) on one row per (station,
        // nonempty support), holding the unit-pole intensities v_{i,r}/mu_{i,S} of
        // that singular hyperplane. Support S is dominated by a superset S' with
        // mu_{i,S'} <= mu_{i,S}, since then v/mu_{S'} >= v/mu_S on all of S; keeping
        // such slack rows would perturb the group averages of eq. 5.44.
        List<double[]> rows = new ArrayList<double[]>();
        for (int i = 0; i < M; i++) {
            for (int mask = 1; mask < nmask; mask++) {
                boolean dominated = false;
                for (int mask2 = 1; mask2 < nmask; mask2++) {
                    if (mask2 != mask && (mask & mask2) == mask
                            && muS[i][mask2] <= muS[i][mask] * (1 + 1e-12)) {
                        dominated = true;
                        break;
                    }
                }
                if (dominated) continue;
                double[] row = new double[Rk];
                for (int t = 0; t < bits[mask].length; t++) {
                    int b = bits[mask][t];
                    row[b] = V[i][b] / muS[i][mask];
                }
                rows.add(row);
            }
        }
        int nrow = rows.size();
        double[][] Lt = new double[nrow][];
        for (int i = 0; i < nrow; i++) {
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

        Ctx ctx = new Ctx(M, Rk, nmask, Nk, lv, rad, arho0, vs, muS, bits, subcol);
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

    // support rate of an OI station, read at its 0/1 support indicator, itself a
    // lattice point of that support since every retained chain has N_r >= 1
    private static double supportRate(ToDoubleFunction<int[]> murate, int[] chi, int i) {
        double rate = murate.applyAsDouble(chi);
        if (!(rate > 0)) {
            throw new RuntimeException("pfqn_clwoi: station " + (i + 1)
                    + " has a non-positive rate on a reachable support.");
        }
        return rate;
    }

    // exhaustive support-only check: every state 0 < n <= N is compared against the
    // rate of its own support. The scan costs prod_r (N_r+1) rate evaluations per
    // station, below the prod_r 2 l_r N_r contour points the inversion itself
    // spends, since l_r >= 1 gives 2 l_r N_r >= N_r + 1
    private static void checkSupportOnly(List<ToDoubleFunction<int[]>> muList, double[][] muS,
                                         int[] keep, int[] Nk, int Rk, int R) {
        int M = muList.size();
        if (M == 0) return;
        long L = 1;
        for (int j = 0; j < Rk; j++) L *= (long) Nk[j] + 1;
        for (int i = 0; i < M; i++) {
            int[] n = new int[R];
            for (long idx = 1; idx < L; idx++) {   // idx 0 is the empty support, unused
                long rem = idx;
                int mask = 0;
                for (int j = 0; j < Rk; j++) {
                    long base = (long) Nk[j] + 1;
                    int nj = (int) (rem % base);
                    rem /= base;
                    n[keep[j]] = nj;
                    if (nj > 0) mask |= 1 << j;
                }
                double rate = muList.get(i).applyAsDouble(n);
                double ref = muS[i][mask];
                if (Math.abs(rate - ref) > 1e-9 * Math.max(1.0, Math.abs(ref))) {
                    int[] chi = new int[R];
                    for (int j = 0; j < Rk; j++) {
                        if ((mask & (1 << j)) != 0) chi[keep[j]] = 1;
                    }
                    throw new RuntimeException("pfqn_clwoi: station " + (i + 1)
                            + " has a rate that varies within a support (mu=" + rate
                            + " at n=" + Arrays.toString(n) + " but mu=" + ref
                            + " at the indicator " + Arrays.toString(chi) + " of the same"
                            + " support). pfqn_clwoi requires order-independent (support-only)"
                            + " rates, mu(n)=mu(supp(n)); a rate that varies inside a support is"
                            + " a general balanced-fairness station and must use pfqn_ncoi.");
                }
            }
        }
    }

    // immutable context shared (read-only) across the recursion
    private static final class Ctx {
        final int M, p, nmask;
        final int[] N, l;
        final double[] r, arho0;
        final double[][] vs, muS;
        final int[][] bits, subcol;
        Ctx(int M, int p, int nmask, int[] N, int[] l, double[] r, double[] arho0,
            double[][] vs, double[][] muS, int[][] bits, int[][] subcol) {
            this.M = M; this.p = p; this.nmask = nmask; this.N = N; this.l = l;
            this.r = r; this.arho0 = arho0; this.vs = vs; this.muS = muS;
            this.bits = bits; this.subcol = subcol;
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
    // Gbar(w) = exp(sum_r arho0_r (w_r-1)) prod_i F_i(alpha_r v_{i,r} w_r), with F_i
    // assembled by the support recursion. Summing logs is legitimate for the
    // principal complex log because exp(log a + log b) = a b.
    private static double[] gbarEval(Ctx ctx, double[] wRe, double[] wIm) {
        double eRe = 0.0, eIm = 0.0;
        for (int j = 0; j < ctx.p; j++) {
            eRe += ctx.arho0[j] * (wRe[j] - 1.0);
            eIm += ctx.arho0[j] * wIm[j];
        }
        double[] fRe = new double[ctx.nmask];
        double[] fIm = new double[ctx.nmask];
        for (int i = 0; i < ctx.M; i++) {
            fRe[0] = 1.0;
            fIm[0] = 0.0;
            double sumRe = 1.0, sumIm = 0.0;
            for (int mask = 1; mask < ctx.nmask; mask++) {
                double numRe = 0.0, numIm = 0.0;
                double denRe = ctx.muS[i][mask], denIm = 0.0;
                for (int t = 0; t < ctx.bits[mask].length; t++) {
                    int b = ctx.bits[mask][t];
                    double xRe = ctx.vs[i][b] * wRe[b];
                    double xIm = ctx.vs[i][b] * wIm[b];
                    int sc = ctx.subcol[mask][t];
                    numRe += xRe * fRe[sc] - xIm * fIm[sc];
                    numIm += xRe * fIm[sc] + xIm * fRe[sc];
                    denRe -= xRe;
                    denIm -= xIm;
                }
                double d2 = denRe * denRe + denIm * denIm;
                fRe[mask] = (numRe * denRe + numIm * denIm) / d2;
                fIm[mask] = (numIm * denRe - numRe * denIm) / d2;
                sumRe += fRe[mask];
                sumIm += fIm[mask];
            }
            eRe += Math.log(Math.hypot(sumRe, sumIm));
            eIm += Math.atan2(sumIm, sumRe);
        }
        double ex = Math.exp(eRe);
        return new double[]{ex * Math.cos(eIm), ex * Math.sin(eIm)};
    }
}
