/**
 * MVAC (Mean Value Analysis by Chain) for load-dependent closed product-form networks
 *
 * Implements the Section V extension of MVAC to queue-length dependent service centers,
 * which propagates the marginal queue-length distributions instead of the mean
 * queue-lengths of {@link jline.api.pfqn.mva.Pfqn_mvac}.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mvacld {
    private Pfqn_mvacld() {}

    /**
     * MVAC for a load-dependent closed network with no infinite-server center.
     */
    public static Ret.pfqnMVACLD pfqn_mvacld(Matrix L, Matrix N, Matrix mu) {
        return pfqn_mvacld(L, N, new Matrix(1, L.getNumCols()), mu);
    }

    /**
     * MVAC (Mean Value Analysis by Chain) for a closed multichain product-form queueing
     * network that may contain queue-length dependent (QLD) service centers. This is the
     * Section V extension of Conway, de Souza e Silva and Lavenberg, IEEE Trans.
     * Computers 38(3):432-442, 1989; {@link jline.api.pfqn.mva.Pfqn_mvac} implements
     * Sections II-IV, which cover single-server fixed-rate (SSFR) and infinite-server
     * (IS) centers only.
     *
     * <p>Where {@code Pfqn_mvac} propagates the MEAN queue-lengths L^k_j(v) through
     * eq. (7) and closes the recursion with the arrival-theorem identity (10), the QLD
     * extension propagates the MARGINAL queue-length DISTRIBUTIONS P^k_j(n,v) instead.
     * That is forced by load dependence -- the rate seen by a job depends on the whole
     * occupancy, so a mean no longer suffices -- but it also SIMPLIFIES the recursion:
     * eq. (21)-(25) read level k-1 only at the shifted vectors v + 1_i, so the basic step
     * sweeps v in I_k alone, where {@code Pfqn_mvac} must sweep the larger
     * I_k u ... u I_K. The marginals come almost for free and are returned as a
     * first-class output.</p>
     *
     * <p>Notation follows {@code Pfqn_mvac}: j, i index centers (j = 1,...,J1 the QLD
     * centers of L, j = J1+1,...,J the IS centers of Z), k, l index the single-customer
     * chains, a_jk = theta_jk T_jk is the demand of chain k at center j, K = sum(N) and
     * I_k = {v : sum_j v_j = K - k} with v_j the number of self-looping single-customer
     * (SCSL) chains pinned at center j. P^k_j(n,v) is the probability of n customers at
     * center j -- EXCLUDING the v_j SCSL customers there -- in the network with
     * normalizing constant G_k(v). Writing tau_k(v,i) for the throughput of an SCSL chain
     * that replaces chain k at center i:</p>
     *
     * <pre>
     *   tau_k(v,i)  = T_ik^-1 sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i)
     *                              * mu_i(n+v_i+1)/(n+v_i+1)                     (21)
     *   tau_k(v,i)  = T_ik^-1,                                     i IS          (22)
     *   L^k_{jk}(v) = theta_jk tau_k(v,j)^-1 / sum_m theta_mk tau_k(v,m)^-1      (23)
     *   lambda^k_k(v) = tau_k(v,j(k)) L^k_{j(k)k}(v)                             (24)
     *   P^k_j(n,v)  = L^k_{jk}(v) P^{k-1}_j(n-1,v+1_j)
     *                 + sum_{m != j} L^k_{mk}(v) P^{k-1}_j(n,v+1_m)              (25)
     * </pre>
     *
     * <p>with P^0_j(0,v) = 1 and P^{k-1}_j(n,.) = 0 for n &lt; 0 or n &gt; k-1. Eq. (21)
     * is just "the mean rate at which an SCSL chain is served": given n other customers
     * the processor-sharing rate share is mu_i(n+v_i+1)/(n+v_i+1), averaged over the
     * distribution of those others. The queueing discipline may be assumed PS with no
     * loss of generality, since product-form measures do not depend on it.</p>
     *
     * <p>This implementation writes (23)-(24) in the reference-station-free form</p>
     *
     * <pre>
     *   c_i(k,v)      = sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i) mu_i(n+v_i+1)/(n+v_i+1)
     *   L^k_{jk}(v)   = (a_jk/c_j) / sum_m (a_mk/c_m)
     *   lambda^k_k(v) = 1 / sum_m (a_mk/c_m)
     * </pre>
     *
     * <p>which follows from theta_jk tau_k(v,j)^-1 = a_jk/c_j and theta_{j(k)k} = 1, so
     * only the demands a_jk are needed and the visit ratios never appear separately. For
     * an IS center c_i = 1 identically, which is exactly (22). Eq. (25) is
     * self-normalizing, sum_n P^k_j(n,v) = sum_m L^k_{mk}(v) = 1, so no normalizing
     * constant is formed and the recursion involves only positive quantities: unlike the
     * classic load-dependent MVA of {@link Pfqn_mvald} it cannot produce negative
     * probabilities and needs no stabilization.</p>
     *
     * <p>Parts 2 and 3 are unchanged from {@code Pfqn_mvac}, since eq. (6) holds verbatim
     * in the presence of QLD centers: part 2 resolves the chains that visit at least one
     * IS center, and the chains that visit no IS center are resolved by re-executing
     * part 1 with their label interchanged with K.</p>
     *
     * @param L  service demand matrix of the QLD centers (M x R)
     * @param N  population vector (1 x R), finite and nonnegative
     * @param Z  demand matrix of the IS centers, (1 x R) or (Mz x R), one row per center
     * @param mu load-dependent rates (M x Nt) with Nt &gt;= sum(N); mu(j,n-1) is the total
     *           service rate of center j with n jobs present. mu(j,:) = 1 is a
     *           single-server fixed-rate queue, mu(j,n-1) = min(n,c) a c-server queue,
     *           mu(j,n-1) = n an infinite server. A null or empty matrix means all
     *           centers are SSFR, in which case results agree with {@code Pfqn_mvac}
     * @return per-class throughput, queue-length, per-station utilization, per-class
     *         cycle time and the marginal queue-length distributions
     */
    public static Ret.pfqnMVACLD pfqn_mvacld(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        Matrix Zl = (Z == null || Z.isEmpty()) ? new Matrix(1, R) : Z;
        if (Zl.getNumCols() != R) {
            throw new IllegalArgumentException(
                    "the think time matrix and the demand matrix have a different number of classes");
        }
        double[] Nd = new double[Math.max(N.getNumRows(), N.getNumCols())];
        if (Nd.length != R) {
            throw new IllegalArgumentException(
                    "the demand matrix and the population vector have a different number of classes");
        }
        for (int r = 0; r < R; r++) {
            Nd[r] = (N.getNumRows() == 1) ? N.get(0, r) : N.get(r, 0);
        }
        int[] Nv = new int[R];
        int K = 0; // number of single-customer chains in the original network
        for (int r = 0; r < R; r++) {
            if (Nd[r] < 0 || Double.isInfinite(Nd[r]) || Double.isNaN(Nd[r])) {
                throw new IllegalArgumentException("the population vector must be finite and nonnegative");
            }
            Nv[r] = (int) Math.round(Nd[r]);
            K += Nv[r];
        }

        Matrix mul = mu;
        if (mul == null || mul.isEmpty()) {
            mul = Matrix.ones(M, Math.max(K, 1));
        }
        if (mul.getNumRows() != M) {
            throw new IllegalArgumentException(
                    "the rate matrix and the demand matrix have a different number of centers");
        }
        if (K > 0 && mul.getNumCols() < K) {
            throw new IllegalArgumentException(
                    "the rate matrix must supply a rate for every population up to sum(N)");
        }

        Matrix XN = new Matrix(1, R);
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, 1);
        Matrix CN = new Matrix(1, R);
        Matrix pij = new Matrix(M, K + 1);
        for (int i = 0; i < M; i++) {
            pij.set(i, 0, 1.0); // an unvisited center holds no jobs
        }
        if (K == 0) {
            return new Ret.pfqnMVACLD(XN, QN, UN, CN, pij);
        }

        // Discard the centers that no chain visits: they carry no customers and would only
        // inflate the multiplicity vector v.
        List<Integer> ldIdx = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) > 0) {
                    ldIdx.add(i);
                    break;
                }
            }
        }
        List<Integer> isIdx = new ArrayList<Integer>();
        for (int i = 0; i < Zl.getNumRows(); i++) {
            for (int r = 0; r < R; r++) {
                if (Zl.get(i, r) > 0) {
                    isIdx.add(i);
                    break;
                }
            }
        }
        int J1 = ldIdx.size();
        int J = J1 + isIdx.size();
        if (J == 0) {
            throw new IllegalArgumentException("all service demands are zero, the throughput is unbounded");
        }
        // A[j][r] = a_jr, centers ordered QLD first, then IS, as required by (21)-(22)
        double[][] A = new double[J][R];
        for (int j = 0; j < J1; j++) {
            for (int r = 0; r < R; r++) {
                A[j][r] = L.get(ldIdx.get(j), r);
            }
        }
        for (int j = J1; j < J; j++) {
            for (int r = 0; r < R; r++) {
                A[j][r] = Zl.get(isIdx.get(j - J1), r);
            }
        }
        double[][] MU = new double[J1][K]; // rates of the QLD centers
        for (int j = 0; j < J1; j++) {
            for (int n = 0; n < K; n++) {
                MU[j][n] = mul.get(ldIdx.get(j), n);
                if (MU[j][n] <= 0) {
                    throw new IllegalArgumentException("the service rates must be strictly "
                            + "positive for every population up to sum(N)");
                }
            }
        }

        // Partition the chains into subsets of identical single-customer chains.
        List<Integer> posr = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (Nv[r] > 0) {
                posr.add(r);
            }
        }
        Map<String, Integer> seen = new HashMap<String, Integer>();
        List<double[]> distinct = new ArrayList<double[]>();
        int[] grpOfClass = new int[posr.size()];
        for (int c = 0; c < posr.size(); c++) {
            int r = posr.get(c);
            double[] col = new double[J];
            StringBuilder key = new StringBuilder();
            for (int j = 0; j < J; j++) {
                col[j] = A[j][r];
                key.append(Double.toString(col[j])).append(',');
            }
            Integer g = seen.get(key.toString());
            if (g == null) {
                g = distinct.size();
                seen.put(key.toString(), g);
                distinct.add(col);
            }
            grpOfClass[c] = g;
        }
        int D = distinct.size();
        // see _kb/03-api-layer.md for rationale
        boolean[] visitsIS = new boolean[D];
        for (int g = 0; g < D; g++) {
            for (int j = J1; j < J; j++) {
                if (distinct.get(g)[j] > 0) {
                    visitsIS[g] = true;
                    break;
                }
            }
        }
        int[] gorder = new int[D];
        int pos = 0;
        for (int g = 0; g < D; g++) {
            if (visitsIS[g]) {
                gorder[pos++] = g;
            }
        }
        int S = 0;
        for (int g = 0; g < D; g++) {
            if (!visitsIS[g]) {
                gorder[pos++] = g;
                S++;
            }
        }
        int[] mult = new int[D];
        for (int g = 0; g < D; g++) {
            for (int c = 0; c < posr.size(); c++) {
                if (grpOfClass[c] == gorder[g]) {
                    mult[g] += Nv[posr.get(c)];
                }
            }
        }

        // Chain labels (0-based here, 1-based in the paper): the D representatives take
        // the labels K-D,...,K-1, the remaining K-D identical chains fill 0,...,K-D-1.
        int[] chainGroup = new int[K];
        for (int g = 0; g < D; g++) {
            chainGroup[K - D + g] = g;
        }
        int p = 0;
        for (int g = 0; g < D; g++) {
            for (int c = 0; c < mult[g] - 1; c++) {
                chainGroup[p++] = g;
            }
        }
        double[][] a = new double[J][K]; // relative utilizations a_jk
        for (int k = 0; k < K; k++) {
            double[] col = distinct.get(gorder[chainGroup[k]]);
            for (int j = 0; j < J; j++) {
                a[j][k] = col[j];
            }
        }

        // see _kb/03-api-layer.md for rationale
        List<int[]> Vlist = new ArrayList<int[]>();
        int[] cnt = new int[K + 1];
        int[] off = new int[K + 1];
        for (int t = 0; t <= K; t++) {
            Matrix mc = Maths.multichoose((double) J, (double) t);
            off[t] = Vlist.size();
            cnt[t] = mc.getNumRows();
            for (int i = 0; i < mc.getNumRows(); i++) {
                int[] v = new int[J];
                for (int j = 0; j < J; j++) {
                    v[j] = (int) Math.round(mc.get(i, j));
                }
                Vlist.add(v);
            }
        }
        int nv = Vlist.size();
        Map<String, Integer> keyOf = new HashMap<String, Integer>();
        for (int vi = 0; vi < nv; vi++) {
            keyOf.put(vecKey(Vlist.get(vi)), vi);
        }
        // succ[vi][i] is the global index of v + 1_i, or -1 when out of range
        int[][] succ = new int[nv][J];
        for (int vi = 0; vi < nv; vi++) {
            int[] v = Vlist.get(vi);
            for (int j = 0; j < J; j++) {
                v[j]++;
                Integer idx = keyOf.get(vecKey(v));
                succ[vi][j] = (idx == null) ? -1 : idx.intValue();
                v[j]--;
            }
        }

        // MVAC recursion
        double[][][][] Pall = new double[K + 1][][][]; // Pall[k][j][vloc][n] = P^k_j(n,v)
        double[][][] Ljkall = new double[K + 1][][]; // Ljkall[k][j][vloc] = L^k_{jk}(v)
        double[][] lamall = new double[K + 1][]; // lamall[k][vloc] = lambda^k_k(v)
        Pall[0] = new double[J1][cnt[K]][1]; // P^0_j(0,v) = 1, I_0 is the block of sum K
        for (int j = 0; j < J1; j++) {
            for (int vloc = 0; vloc < cnt[K]; vloc++) {
                Pall[0][j][vloc][0] = 1.0;
            }
        }

        double[] lamChain = new double[K]; // per-chain throughput, by ORIGINAL chain label
        double[][] Lchain = new double[J][K]; // per-chain mean queue-length, same indexing

        // First execution of the basic step, k0 = 1
        part1(1, K, J, J1, a, MU, Vlist, off, cnt, succ, Pall, Ljkall, lamall);
        lamChain[K - 1] = lamall[K][0];
        for (int j = 0; j < J; j++) {
            Lchain[j][K - 1] = Ljkall[K][j][0];
        }
        // The marginals of the ORIGINAL network are read at k = K and v = 0; the label
        // interchanges below overwrite this level, so capture them now.
        for (int j = 0; j < J1; j++) {
            for (int n = 0; n <= K; n++) {
                pij.set(ldIdx.get(j), n, Pall[K][j][0][n]);
            }
        }

        // Part 2: measures of the chains that visit at least one IS center. Eq. (6) is
        // unchanged by load dependence.
        int lmaxK = Math.min(K - 1, K - S);
        if (D >= 2 && lmaxK >= K - D + 1) {
            double[][][] L2prev = new double[K + 1][][];
            for (int k = K - D + 2; k <= K; k++) {
                int t = K - k;
                double[][][] L2cur = new double[K + 1][][];
                for (int l = K - D + 1; l <= Math.min(k - 1, K - S); l++) {
                    double[][] acc = new double[J][cnt[t]];
                    for (int vloc = 0; vloc < cnt[t]; vloc++) {
                        int vi = off[t] + vloc;
                        for (int j = 0; j < J; j++) {
                            int sloc = succ[vi][j] - off[t + 1];
                            // base case of (6): L^{k-1}_{i,k-1} comes from (23)
                            double[][] prev = (l == k - 1) ? Ljkall[k - 1] : L2prev[l];
                            for (int i = 0; i < J; i++) {
                                acc[i][vloc] += Ljkall[k][j][vloc] * prev[i][sloc];
                            }
                        }
                    }
                    L2cur[l] = acc;
                }
                if (k == K) {
                    for (int l = K - D + 1; l <= lmaxK; l++) {
                        for (int i = 0; i < J; i++) {
                            Lchain[i][l - 1] = L2cur[l][i][0];
                        }
                        // Little's law at an IS center visited by chain l
                        int jIS = -1;
                        for (int j = J1; j < J; j++) {
                            if (a[j][l - 1] > 0) {
                                jIS = j;
                                break;
                            }
                        }
                        lamChain[l - 1] = Lchain[jIS][l - 1] / a[jIS][l - 1];
                    }
                }
                L2prev = L2cur;
            }
        }

        // see _kb/03-api-layer.md for rationale
        int[] perm = new int[K];
        for (int k = 0; k < K; k++) {
            perm[k] = k;
        }
        for (int l = 1; l <= S - 1; l++) {
            int lo = K - l - 1;
            int hi = K - 1;
            int tmpi = perm[lo];
            perm[lo] = perm[hi];
            perm[hi] = tmpi;
            for (int j = 0; j < J; j++) {
                double tmp = a[j][lo];
                a[j][lo] = a[j][hi];
                a[j][hi] = tmp;
            }
            part1(K - l, K, J, J1, a, MU, Vlist, off, cnt, succ, Pall, Ljkall, lamall);
            lamChain[perm[hi]] = lamall[K][0];
            for (int j = 0; j < J; j++) {
                Lchain[j][perm[hi]] = Ljkall[K][j][0];
            }
        }

        // see _kb/03-api-layer.md for rationale
        for (int g = 0; g < D; g++) {
            int kg = K - D + g; // original label of the representative of subset g
            for (int c = 0; c < posr.size(); c++) {
                if (grpOfClass[c] != gorder[g]) {
                    continue;
                }
                int r = posr.get(c);
                XN.set(0, r, Nv[r] * lamChain[kg]);
                for (int j = 0; j < J1; j++) {
                    int i = ldIdx.get(j);
                    QN.set(i, r, Nv[r] * Lchain[j][kg]);
                }
            }
        }
        // Utilization of a load-dependent center is 1-P_j(0), NOT XN(r)*L(j,r).
        for (int i = 0; i < M; i++) {
            UN.set(i, 0, 1.0 - pij.get(i, 0));
        }
        // Cycle time exclusive of think time, as in Pfqn_mvald and Pfqn_dac. An empty
        // class has zero throughput, hence no cycle time.
        for (int r = 0; r < R; r++) {
            if (Nv[r] > 0) {
                double zr = 0;
                for (int i = 0; i < Zl.getNumRows(); i++) {
                    zr += Zl.get(i, r);
                }
                CN.set(0, r, Nv[r] / XN.get(0, r) - zr);
            }
        }

        return new Ret.pfqnMVACLD(XN, QN, UN, CN, pij);
    }

    /**
     * Part 1 of the basic step for networks with QLD centers: evaluate (21)-(25) for
     * k = k0,...,K over v in I_k, the block of multiplicity vectors of sum K-k.
     */
    private static void part1(int k0, int K, int J, int J1, double[][] a, double[][] MU,
                              List<int[]> Vlist, int[] off, int[] cnt, int[][] succ,
                              double[][][][] Pall, double[][][] Ljkall, double[][] lamall) {
        for (int k = k0; k <= K; k++) {
            int t = K - k;
            int nvk = cnt[t];
            double[][][] Pp = Pall[k - 1]; // P^{k-1}, over I_{k-1}, last index n = 0..k-1
            double[][] Lk = new double[J][nvk];
            double[] lam = new double[nvk];
            double[][][] Pk = new double[J1][nvk][k + 1];
            int[] sloc = new int[J];
            for (int vloc = 0; vloc < nvk; vloc++) {
                int vi = off[t] + vloc;
                int[] v = Vlist.get(vi);
                for (int i = 0; i < J; i++) {
                    sloc[i] = succ[vi][i] - off[t + 1]; // local index of v + 1_i in I_{k-1}
                }
                // see _kb/03-api-layer.md for rationale
                double[] c = new double[J];
                for (int i = 0; i < J; i++) {
                    c[i] = 1.0;
                }
                for (int i = 0; i < J1; i++) {
                    double ci = 0;
                    for (int n = 0; n <= k - 1; n++) {
                        ci += Pp[i][sloc[i]][n] * MU[i][n + v[i]] / (n + v[i] + 1);
                    }
                    c[i] = ci;
                }
                // (23)-(24) in reference-station-free form: theta_jk/tau_k(v,j) = a_jk/c_j
                double[] w = new double[J];
                double sw = 0;
                for (int j = 0; j < J; j++) {
                    if (a[j][k - 1] > 0) {
                        w[j] = a[j][k - 1] / c[j];
                        sw += w[j];
                    }
                }
                lam[vloc] = 1.0 / sw;
                for (int j = 0; j < J; j++) {
                    Lk[j][vloc] = w[j] / sw;
                }
                // (25): condition on the center holding the single chain-k customer
                for (int j = 0; j < J1; j++) {
                    for (int n = 0; n <= k; n++) {
                        double s = 0;
                        if (n >= 1) {
                            // chain k is at j, so n-1 of the k-1 others are there too
                            s = Lk[j][vloc] * Pp[j][sloc[j]][n - 1];
                        }
                        if (n <= k - 1) {
                            // chain k is elsewhere, so all n are from the k-1 others
                            for (int m = 0; m < J; m++) {
                                if (m != j) {
                                    s += Lk[m][vloc] * Pp[j][sloc[m]][n];
                                }
                            }
                        }
                        Pk[j][vloc][n] = s;
                    }
                }
            }
            Pall[k] = Pk;
            Ljkall[k] = Lk;
            lamall[k] = lam;
        }
    }

    private static String vecKey(int[] v) {
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < v.length; j++) {
            sb.append(v[j]).append(',');
        }
        return sb.toString();
    }
}
