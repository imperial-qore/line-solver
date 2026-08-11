/**
 * Distribution Analysis by Chain (DAC)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * DAC (Distribution Analysis by Chain) method for closed product-form queueing
 * networks with single-server fixed-rate, infinite-server and queue-dependent
 * service centers.
 *
 * <p>Unlike MVA, RECAL or MVAC, the recursion returns the whole set of joint
 * queue-length probabilities, which are required e.g. in availability modeling.
 * It proceeds chain by chain over a related network in which every chain holds a
 * single customer, a transformation that leaves the aggregate queue-length
 * distribution unchanged. Given the distribution of a network with k-1 such
 * chains, adding one customer of a chain with demands r gives
 *
 * <pre>
 *   c_j      = sum_{n=1..k} (n/mu_j(n)) * P_j^{k-1}(n-1)
 *   lambda_k = 1 / sum_j r_j c_j
 *   P^k(n)   = lambda_k * sum_j r_j (n_j/mu_j(n_j)) * P^{k-1}(n-e_j)
 * </pre>
 *
 * where lambda_k is the throughput of the customer being added and 1/c_j is the
 * throughput of a chain visiting center j only. The recursion conserves
 * probability mass by construction, hence it is numerically stable.
 *
 * <p>References: E. de Souza e Silva, "Distribution Analysis of Product Form
 * Queueing Networks", UCLA Computer Science Department, CSD-870023, April 1987.
 */
public final class Pfqn_dac {
    private Pfqn_dac() {}

    /**
     * DAC with default think times (zero) and default rates (single-server fixed rate).
     *
     * @param L service demand matrix (M x R)
     * @param N population vector (1 x R)
     * @return the joint queue-length distribution and the mean measures
     */
    public static Ret.pfqnDAC pfqn_dac(Matrix L, Matrix N) {
        return pfqn_dac(L, N, null, null);
    }

    /**
     * DAC (Distribution Analysis by Chain) for joint queue-length distributions.
     *
     * @param L  service demand matrix (M x R)
     * @param N  population vector (1 x R)
     * @param Z  think time vector (1 x R), or null for zero think times. If the
     *           think times are non-zero an extra infinite-server station is
     *           appended, so that the states matrix has M+1 columns and its last
     *           column holds the think-station population.
     * @param mu load-dependent rate matrix (M x Nt), Nt=sum(N), or null for
     *           single-server fixed rate. Use mu(j,n)=n+1 for infinite server and
     *           mu(j,n)=min(n+1,c) for a c-server station.
     * @return the joint queue-length distribution and the mean measures
     */
    public static Ret.pfqnDAC pfqn_dac(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        double[] Nv = new double[R];
        for (int r = 0; r < R; r++) {
            Nv[r] = N.get(r);
            if (Nv[r] < 0) {
                throw new IllegalArgumentException("pfqn_dac: population vector must be non-negative.");
            }
        }
        double[] Zv = new double[R];
        if (Z != null) {
            for (int r = 0; r < R; r++) {
                Zv[r] = Z.get(r);
            }
        }
        int Nt = 0;
        for (int r = 0; r < R; r++) {
            Nt += (int) Math.round(Nv[r]);
        }
        int W = Math.max(Nt, 1);

        if (mu != null) {
            if (mu.getNumRows() != M) {
                throw new IllegalArgumentException("pfqn_dac: the mu matrix must have one row per station.");
            }
            if (Nt > 0 && mu.getNumCols() < Nt) {
                throw new IllegalArgumentException("pfqn_dac: the mu matrix must have at least sum(N) columns.");
            }
        }

        // A non-zero think time is modelled as an appended infinite-server center.
        double zsum = 0;
        for (int r = 0; r < R; r++) {
            zsum += Zv[r];
        }
        boolean hasZ = zsum > 0;
        int J = hasZ ? M + 1 : M;

        double[][] Lx = new double[J][R];
        for (int j = 0; j < M; j++) {
            for (int r = 0; r < R; r++) {
                Lx[j][r] = L.get(j, r);
            }
        }
        double[][] mux = new double[J][W];
        for (int j = 0; j < M; j++) {
            for (int n = 0; n < W; n++) {
                mux[j][n] = (mu == null) ? 1.0 : mu.get(j, n);
            }
        }
        if (hasZ) {
            for (int r = 0; r < R; r++) {
                Lx[M][r] = Zv[r];
            }
            for (int n = 0; n < W; n++) {
                mux[M][n] = n + 1;
            }
        }

        if (Nt == 0) {
            Matrix pi0 = Matrix.zeros(M, 1);
            for (int j = 0; j < M; j++) {
                pi0.set(j, 0, 1.0);
            }
            Matrix pj = new Matrix(1, 1);
            pj.set(0, 0, 1.0);
            return new Ret.pfqnDAC(pj, Matrix.zeros(1, J), Matrix.zeros(1, R), Matrix.zeros(M, R),
                    Matrix.zeros(M, 1), Matrix.zeros(1, R), pi0);
        }

        List<Integer> active = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (Nv[r] > 0) {
                active.add(Integer.valueOf(r));
            }
        }
        int D = active.size();
        for (int idx = 0; idx < D; idx++) {
            int r = active.get(idx).intValue();
            boolean any = false;
            for (int j = 0; j < J; j++) {
                if (Lx[j][r] > 0) {
                    any = true;
                    break;
                }
            }
            if (!any) {
                throw new IllegalArgumentException("pfqn_dac: chain " + (r + 1) + " has null demand at every center.");
            }
        }

        Lattice lat = new Lattice(J, Nt);

        // Chains are ordered so that the last D added are distinct, which lets all
        // per-chain measures reuse the common prefix of the recursion.
        List<Integer> prefix = new ArrayList<Integer>();
        for (int idx = 0; idx < D; idx++) {
            int r = active.get(idx).intValue();
            for (int c = 0; c < (int) Math.round(Nv[r]) - 1; c++) {
                prefix.add(Integer.valueOf(r));
            }
        }

        // Prefix of the recursion, shared by every per-chain run.
        double[] p = new double[]{1.0};
        int k = 0;
        for (int idx = 0; idx < prefix.size(); idx++) {
            Step st = step(p, Lx, prefix.get(idx).intValue(), mux, k, lat);
            p = st.p;
            k++;
        }

        // Base run over the tail, saving the intermediate distributions S_0..S_{D-1}.
        double[][] Sp = new double[D][];
        Sp[0] = p;
        double[] pb = p;
        int kb = k;
        Step last = null;
        for (int idx = 0; idx < D; idx++) {
            last = step(pb, Lx, active.get(idx).intValue(), mux, kb, lat);
            pb = last.p;
            kb++;
            if (idx < D - 1) {
                Sp[idx + 1] = pb;
            }
        }

        Matrix XN = Matrix.zeros(1, R);
        Matrix QN = Matrix.zeros(M, R);
        // The base run already places the last active chain last.
        int rlast = active.get(D - 1).intValue();
        XN.set(0, rlast, Nv[rlast] * last.lam);
        for (int j = 0; j < M; j++) {
            QN.set(j, rlast, Nv[rlast] * last.Lq[j]);
        }

        // Re-run the tail with chain active(idx) moved last, restarting from S_idx.
        for (int idx = 0; idx < D - 1; idx++) {
            double[] pc = Sp[idx];
            int kc = k + idx;
            List<Integer> order = new ArrayList<Integer>();
            for (int t = idx + 1; t < D; t++) {
                order.add(active.get(t));
            }
            order.add(active.get(idx));
            Step cur = null;
            for (int t = 0; t < order.size(); t++) {
                cur = step(pc, Lx, order.get(t).intValue(), mux, kc, lat);
                pc = cur.p;
                kc++;
            }
            int r = active.get(idx).intValue();
            XN.set(0, r, Nv[r] * cur.lam);
            for (int j = 0; j < M; j++) {
                QN.set(j, r, Nv[r] * cur.Lq[j]);
            }
        }

        int[][] states = lat.states[Nt];
        Matrix Pjoint = new Matrix(states.length, 1);
        Matrix stM = new Matrix(states.length, J);
        for (int i = 0; i < states.length; i++) {
            Pjoint.set(i, 0, pb[i]);
            for (int j = 0; j < J; j++) {
                stM.set(i, j, states[i][j]);
            }
        }

        // Marginal queue-length probabilities at the full population.
        Matrix pi = Matrix.zeros(M, Nt + 1);
        for (int i = 0; i < states.length; i++) {
            for (int j = 0; j < M; j++) {
                pi.set(j, states[i][j], pi.get(j, states[i][j]) + pb[i]);
            }
        }
        Matrix UN = Matrix.zeros(M, 1);
        for (int j = 0; j < M; j++) {
            UN.set(j, 0, 1.0 - pi.get(j, 0));
        }
        Matrix CN = Matrix.zeros(1, R);
        for (int idx = 0; idx < D; idx++) {
            int r = active.get(idx).intValue();
            CN.set(0, r, Nv[r] / XN.get(0, r) - Zv[r]);
        }

        return new Ret.pfqnDAC(Pjoint, stM, XN, QN, UN, CN, pi);
    }

    /**
     * Result of a single step of the recursion.
     */
    private static class Step {
        double[] p;
        double lam;
        double[] Lq;

        Step(double[] p, double lam, double[] Lq) {
            this.p = p;
            this.lam = lam;
            this.Lq = Lq;
        }
    }

    /**
     * One step of the recursion: add a single customer of chain c to a network
     * holding k customers, returning the distribution at level k+1, the throughput
     * of the added customer and its per-center presence probabilities.
     */
    private static Step step(double[] p, double[][] Lx, int c, double[][] mux, int k, Lattice lat) {
        int J = Lx.length;
        int[][] S = lat.states[k];

        // Marginal queue lengths of the network with k customers.
        double[][] marg = new double[J][k + 1];
        for (int i = 0; i < S.length; i++) {
            for (int j = 0; j < J; j++) {
                marg[j][S[i][j]] += p[i];
            }
        }
        double[] cj = new double[J];
        for (int j = 0; j < J; j++) {
            for (int n = 1; n <= k + 1; n++) {
                cj[j] += (n / mux[j][n - 1]) * marg[j][n - 1];
            }
        }
        double den = 0;
        for (int j = 0; j < J; j++) {
            den += Lx[j][c] * cj[j];
        }
        double lam = 1.0 / den;
        double[] Lq = new double[J];
        for (int j = 0; j < J; j++) {
            Lq[j] = lam * Lx[j][c] * cj[j];
        }

        int[][] ix = lat.succ[k];
        double[] pn = new double[lat.states[k + 1].length];
        for (int i = 0; i < S.length; i++) {
            for (int j = 0; j < J; j++) {
                if (Lx[j][c] <= 0) {
                    continue;
                }
                int nj = S[i][j] + 1;
                pn[ix[i][j]] += lam * Lx[j][c] * (nj / mux[j][nj - 1]) * p[i];
            }
        }
        return new Step(pn, lam, Lq);
    }

    /**
     * Aggregate state space, enumerated level by level, together with the index of
     * each state with one more job at a given center.
     */
    private static class Lattice {
        final int[][][] states;
        final int[][][] succ;
        private final int J;
        private final long[][] C;

        Lattice(int J, int Nt) {
            this.J = J;
            this.states = new int[Nt + 1][][];
            for (int k = 0; k <= Nt; k++) {
                this.states[k] = compositions(J, k);
            }
            // Binomial table, C[a][b] = binom(a,b), for ranking compositions.
            this.C = new long[Nt + J + 1][J + 2];
            for (int a = 0; a <= Nt + J; a++) {
                for (int b = 0; b <= Math.min(a, J); b++) {
                    if (b == 0) {
                        C[a][b] = 1;
                    } else {
                        C[a][b] = C[a - 1][b - 1] + C[a - 1][b];
                    }
                }
            }
            this.succ = new int[Nt + 1][][];
            for (int k = 0; k < Nt; k++) {
                int[][] S = states[k];
                int[][] ix = new int[S.length][J];
                for (int i = 0; i < S.length; i++) {
                    for (int j = 0; j < J; j++) {
                        int[] t = S[i].clone();
                        t[j]++;
                        ix[i][j] = rank(t, k + 1);
                    }
                }
                succ[k] = ix;
            }
        }

        /**
         * Lexicographic rank of a composition of k into J parts.
         */
        private int rank(int[] n, int k) {
            int idx = 0;
            int rem = k;
            for (int j = 0; j < J - 1; j++) {
                int parts = J - 1 - j;
                for (int v = 0; v < n[j]; v++) {
                    idx += (int) C[(rem - v) + parts - 1][parts - 1];
                }
                rem -= n[j];
            }
            return idx;
        }

        /**
         * All J-part compositions of k, in lexicographic order.
         */
        private static int[][] compositions(int J, int k) {
            if (J == 1) {
                return new int[][]{{k}};
            }
            List<int[]> out = new ArrayList<int[]>();
            for (int v = 0; v <= k; v++) {
                int[][] B = compositions(J - 1, k - v);
                for (int i = 0; i < B.length; i++) {
                    int[] row = new int[J];
                    row[0] = v;
                    System.arraycopy(B[i], 0, row, 1, J - 1);
                    out.add(row);
                }
            }
            return out.toArray(new int[out.size()][]);
        }
    }
}
