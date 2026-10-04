/**
 * MVAC (Mean Value Analysis by Chain) for closed product-form queueing networks
 *
 * Implements the exact mean value analysis by chain of Conway, de Souza e Silva and
 * Lavenberg, which recurs on the chains rather than on the population vector and is
 * therefore attractive for networks with few service centers and many chains.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mvac {
    private Pfqn_mvac() {}

    /**
     * MVAC algorithm for a closed product-form network with no infinite-server center.
     */
    public static Ret.pfqnMVAC pfqn_mvac(Matrix L, Matrix N) {
        return pfqn_mvac(L, N, new Matrix(1, L.getNumCols()));
    }

    /**
     * MVAC (Mean Value Analysis by Chain) algorithm for closed multichain product-form
     * queueing networks composed of single-server fixed-rate (SSFR) queues and
     * infinite-server (IS) centers.
     *
     * <p>Unlike the classic MVA recursion of {@link Pfqn_mva}, which recurs on the
     * population vector and costs O(prod(N+1)), MVAC recurs on the chains: each class is
     * reduced to single-customer chains and the removed chains are replaced by
     * self-looping single-customer (SCSL) chains pinned at a service center, so the
     * multiplicity vector v = (v_1,...,v_J), with v_j the number of SCSL chains at center
     * j, indexes the recursion in place of the population vector. MVAC is thus attractive
     * for networks with few centers and many chains, and is the mean-value counterpart of
     * the RECAL normalizing-constant recursion of
     * {@link jline.api.pfqn.nc.Pfqn_recal}. Since no normalizing constant is formed, MVAC
     * does not suffer the floating-point underflow/overflow that complicates RECAL and
     * convolution.</p>
     *
     * <p>With j, i indexing centers (j = 1,...,J1 SSFR and j = J1+1,...,J IS), k, l
     * indexing single-customer chains, a_jk the relative utilization (service demand) of
     * chain k at center j, a_k = sum_j a_jk, K the total number of single-customer chains
     * and I_k = {v : sum_j v_j = K - k}, the recursion of Section II of the paper reads</p>
     *
     * <pre>
     *   lambda^k_k(v) = 1 / (a_k + sum_{j=1}^{J1} a_jk (L^{k-1}_j(v) + v_j))     (10)
     *   L^k_{jk}(v)   = lambda^k_k(v) a_jk (1 + L^{k-1}_j(v) + v_j),  j SSFR     (9a)
     *   L^k_{jk}(v)   = lambda^k_k(v) a_jk,                           j IS       (9b)
     *   L^k_i(v)      = sum_j L^k_{jk}(v) L^{k-1}_i(v + 1_j) + L^k_{ik}(v)       (7)
     *   L^k_{il}(v)   = sum_j L^k_{jk}(v) L^{k-1}_{il}(v + 1_j),  l = 1,...,k-1  (6)
     * </pre>
     *
     * <p>with L^0_j(v) = 0, where L^k_j(v) is the mean number of customers at center j
     * (SCSL customers excluded), L^k_{jl}(v) the mean number of chain-l customers at
     * center j, and lambda^k_k(v) the throughput of chain k, all for the network with
     * normalizing constant G_k(v). Equation (10) is the arrival-theorem closure obtained
     * by summing (9a)-(9b) over all centers, since chain k holds a single customer. The
     * measures of the original network are read off at k = K and v = 0.</p>
     *
     * <p>Part 1 of the basic step evaluates (10), (9) and (7) and yields the measures of
     * chain K; part 2 evaluates (6) and yields those of the chains that visit at least one
     * IS center, whose throughput then follows from Little's law at that center. Chains
     * that visit only SSFR centers require a re-execution of part 1 with their label
     * interchanged with K, which is cheap because the levels below the interchanged label
     * are unaffected and are reused from the first execution. Classes with N_r &gt; 1, and
     * classes with identical demand columns, collapse into a single subset of identical
     * single-customer chains: only one representative per subset is analyzed and its
     * per-chain measures are scaled by the class population, so the cost depends on the
     * number D of distinct chains rather than on K.</p>
     *
     * <p>Reference: A. E. Conway, E. de Souza e Silva and S. S. Lavenberg, "Mean Value
     * Analysis by Chain of Product Form Queueing Networks", IEEE Trans. Computers,
     * 38(3):432-442, 1989.</p>
     *
     * @param L service demand matrix of the SSFR queues (M x R)
     * @param N population vector (1 x R), finite and nonnegative
     * @param Z service demand matrix of the IS centers, (1 x R) or (Mz x R), one row per
     *          IS center
     * @return the per-class throughput, queue-length, utilization and residence time
     */
    public static Ret.pfqnMVAC pfqn_mvac(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        Matrix Zl = (Z == null || Z.isEmpty()) ? new Matrix(1, R) : Z;
        if (Zl.getNumCols() != R) {
            throw new IllegalArgumentException(
                    "the think time matrix and the demand matrix have a different number of classes");
        }
        // Accept the population vector in either orientation
        double[] Nd = new double[Math.max(N.getNumRows(), N.getNumCols())];
        if (Nd.length != R) {
            throw new IllegalArgumentException(
                    "the demand matrix and the population vector have a different number of classes");
        }
        for (int r = 0; r < R; r++) {
            Nd[r] = (N.getNumRows() == 1) ? N.get(0, r) : N.get(r, 0);
        }
        int[] Nv = new int[R];
        for (int r = 0; r < R; r++) {
            if (Nd[r] < 0 || Double.isInfinite(Nd[r]) || Double.isNaN(Nd[r])) {
                throw new IllegalArgumentException("the population vector must be finite and nonnegative");
            }
            Nv[r] = (int) Math.round(Nd[r]);
        }

        Matrix XN = new Matrix(1, R);
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        Matrix CN = new Matrix(M, R);

        int K = 0; // number of single-customer chains in the original network
        for (int r = 0; r < R; r++) {
            K += Nv[r];
        }
        if (K == 0) {
            return new Ret.pfqnMVAC(XN, QN, UN, CN);
        }

        // Discard the centers that no chain visits: they carry no customers and would only
        // inflate the multiplicity vector v.
        List<Integer> ssfrIdx = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) > 0) {
                    ssfrIdx.add(i);
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
        int J1 = ssfrIdx.size();
        int J = J1 + isIdx.size();
        if (J == 0) {
            throw new IllegalArgumentException("all service demands are zero, the throughput is unbounded");
        }
        // A[j][r] = a_jr, centers ordered SSFR first as required by (9a)-(9b) and (10)
        double[][] A = new double[J][R];
        for (int j = 0; j < J1; j++) {
            for (int r = 0; r < R; r++) {
                A[j][r] = L.get(ssfrIdx.get(j), r);
            }
        }
        for (int j = J1; j < J; j++) {
            for (int r = 0; r < R; r++) {
                A[j][r] = Zl.get(isIdx.get(j - J1), r);
            }
        }

        // see _kb/03-api-layer.md for rationale
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

        // see _kb/03-api-layer.md for rationale
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
        double[] ak = new double[K]; // a_k = sum_j a_jk, over ALL centers
        for (int k = 0; k < K; k++) {
            for (int j = 0; j < J; j++) {
                ak[k] += a[j][k];
            }
        }

        // see _kb/03-api-layer.md for rationale
        List<int[]> Vlist = new ArrayList<int[]>();
        for (int t = 0; t <= K; t++) {
            Matrix mc = Maths.multichoose((double) J, (double) t);
            for (int i = 0; i < mc.getNumRows(); i++) {
                int[] v = new int[J];
                for (int j = 0; j < J; j++) {
                    v[j] = (int) Math.round(mc.get(i, j));
                }
                Vlist.add(v);
            }
        }
        int nv = Vlist.size();
        int[] vsum = new int[nv];
        Map<String, Integer> keyOf = new HashMap<String, Integer>();
        for (int vi = 0; vi < nv; vi++) {
            int[] v = Vlist.get(vi);
            int s = 0;
            for (int j = 0; j < J; j++) {
                s += v[j];
            }
            vsum[vi] = s;
            keyOf.put(vecKey(v), vi);
        }
        // succ[vi][j] is the index of v + 1_j, or -1 when out of range (sum(v) == K), in
        // which case it is never read
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
        int z0 = -1; // index of v = 0
        for (int vi = 0; vi < nv; vi++) {
            if (vsum[vi] == 0) {
                z0 = vi;
                break;
            }
        }

        // MVAC recursion
        double[][][] Lall = new double[K + 1][][]; // Lall[k][i][vi] = L^k_i(v)
        double[][][] Ljkall = new double[K + 1][][]; // Ljkall[k][j][vi] = L^k_{jk}(v)
        double[][] lamall = new double[K + 1][]; // lamall[k][vi] = lambda^k_k(v)
        Lall[0] = new double[J][nv]; // L^0_i(v) = 0

        double[] lamChain = new double[K]; // per-chain throughput, by ORIGINAL chain label
        double[][] Lchain = new double[J][K]; // per-chain mean queue-length, same indexing

        // First execution of the basic step, k0 = 1
        part1(1, K, J, J1, a, ak, Vlist, vsum, succ, Lall, Ljkall, lamall);
        lamChain[K - 1] = lamall[K][z0];
        for (int j = 0; j < J; j++) {
            Lchain[j][K - 1] = Ljkall[K][j][z0];
        }

        // Part 2: measures of the chains that visit at least one IS center
        int lmaxK = Math.min(K - 1, K - S);
        if (D >= 2 && lmaxK >= K - D + 1) {
            double[][][] L2prev = new double[K + 1][][];
            for (int k = K - D + 2; k <= K; k++) {
                List<Integer> idxI = new ArrayList<Integer>(); // v in I_k
                for (int vi = 0; vi < nv; vi++) {
                    if (vsum[vi] == K - k) {
                        idxI.add(vi);
                    }
                }
                double[][][] L2cur = new double[K + 1][][];
                for (int l = K - D + 1; l <= Math.min(k - 1, K - S); l++) {
                    double[][] acc = new double[J][nv];
                    for (int ii = 0; ii < idxI.size(); ii++) {
                        int vi = idxI.get(ii);
                        for (int j = 0; j < J; j++) {
                            int sj = succ[vi][j];
                            // base case of (6): L^{k-1}_{i,k-1} comes from (9)
                            double[][] prev = (l == k - 1) ? Ljkall[k - 1] : L2prev[l];
                            for (int i = 0; i < J; i++) {
                                acc[i][vi] += Ljkall[k][j][vi] * prev[i][sj];
                            }
                        }
                    }
                    L2cur[l] = acc;
                }
                if (k == K) {
                    for (int l = K - D + 1; l <= lmaxK; l++) {
                        for (int i = 0; i < J; i++) {
                            Lchain[i][l - 1] = L2cur[l][i][z0];
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
        int[] perm = new int[K]; // perm[k] is the original label of the chain now at k
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
            double tmpd = ak[lo];
            ak[lo] = ak[hi];
            ak[hi] = tmpd;
            part1(K - l, K, J, J1, a, ak, Vlist, vsum, succ, Lall, Ljkall, lamall);
            lamChain[perm[hi]] = lamall[K][z0];
            for (int j = 0; j < J; j++) {
                Lchain[j][perm[hi]] = Ljkall[K][j][z0];
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
                    int i = ssfrIdx.get(j);
                    QN.set(i, r, Nv[r] * Lchain[j][kg]);
                    UN.set(i, r, XN.get(0, r) * L.get(i, r));
                    CN.set(i, r, QN.get(i, r) / XN.get(0, r));
                }
            }
        }
        // An empty class has zero throughput and queue-length, so its residence time is
        // reported as the bare service demand, as in Pfqn_mva.
        for (int r = 0; r < R; r++) {
            if (Nv[r] == 0) {
                for (int i = 0; i < M; i++) {
                    CN.set(i, r, L.get(i, r));
                }
            }
        }

        return new Ret.pfqnMVAC(XN, QN, UN, CN);
    }

    /**
     * Part 1 of the MVAC basic step: evaluate (10), (9a)-(9b) and (7) for k = k0,...,K
     * over v in V_k = I_k u ... u I_K.
     */
    private static void part1(int k0, int K, int J, int J1, double[][] a, double[] ak,
                              List<int[]> Vlist, int[] vsum, int[][] succ,
                              double[][][] Lall, double[][][] Ljkall, double[][] lamall) {
        int nv = Vlist.size();
        for (int k = k0; k <= K; k++) {
            double[][] Lp = Lall[k - 1]; // L^{k-1}
            double[][] Lk = new double[J][nv];
            double[][] Ljk = new double[J][nv];
            double[] lamv = new double[nv];
            for (int vi = 0; vi < nv; vi++) {
                if (vsum[vi] > K - k) {
                    continue;
                }
                int[] v = Vlist.get(vi);
                // (10): the chain-k customer is at some center, so sum_j L^k_{jk}(v) = 1
                double den = ak[k - 1];
                for (int j = 0; j < J1; j++) {
                    den += a[j][k - 1] * (Lp[j][vi] + v[j]);
                }
                double lam = 1.0 / den;
                lamv[vi] = lam;
                for (int j = 0; j < J1; j++) {
                    Ljk[j][vi] = lam * a[j][k - 1] * (1.0 + Lp[j][vi] + v[j]); // (9a)
                }
                for (int j = J1; j < J; j++) {
                    Ljk[j][vi] = lam * a[j][k - 1]; // (9b)
                }
                // (7)
                for (int i = 0; i < J; i++) {
                    double s = Ljk[i][vi];
                    for (int j = 0; j < J; j++) {
                        s += Ljk[j][vi] * Lp[i][succ[vi][j]];
                    }
                    Lk[i][vi] = s;
                }
            }
            Lall[k] = Lk;
            Ljkall[k] = Ljk;
            lamall[k] = lamv;
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
