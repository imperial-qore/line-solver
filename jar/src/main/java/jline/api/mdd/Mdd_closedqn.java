package jline.api.mdd;

import java.util.ArrayList;
import java.util.List;

import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_solve;
import jline.util.matrix.Matrix;

/**
 * Exact solve of a single-class closed exponential queueing network whose CTMC
 * state space (reachable occupancy vectors) is stored in a Multi-valued
 * Decision Diagram instead of an explicit state list.
 *
 * <p>Port of matlab/src/api/mdd/mdd_closedqn.m. The reachable set is generated
 * with {@link Mdd_reachset} and the generator matrix is assembled using the
 * MDD's O(K) state indexing ({@link MDD#index}), so no explicit (|S| x width)
 * state matrix is ever materialised during assembly -- the diagram is the
 * store. For single-class exponential stations the aggregated (occupancy)
 * chain is exact: the rate from n to n-e_i+e_j is mu_i * min(n_i, c_i) *
 * P(i,j) for n_i &gt; 0, matching SolverCTMC on the same model, which makes
 * this the live exact oracle the {@link Mdd_mcd} aggregation is validated
 * against.</p>
 *
 * <p>The reference's 'verbose' knob is not carried (the api layer is silent);
 * the same storage numbers are returned in {@link MddClosedQnResult#stats}.</p>
 */
public class Mdd_closedqn {

    private Mdd_closedqn() {}

    /**
     * @param mu per-station exponential service rates, length M
     * @param P M x M Markovian routing matrix (row-stochastic, irreducible)
     * @param servers servers per station; Double.POSITIVE_INFINITY for a delay
     * @param N closed population
     */
    public static MddClosedQnResult mdd_closedqn(double[] mu, Matrix P, double[] servers, int N) {
        return mdd_closedqn(mu, P, servers, N, null);
    }

    /**
     * @param reuse an already-built reachable set (the {@code mdd} of a previous
     *              result on the same mu/P/servers/N) to skip regeneration
     */
    public static MddClosedQnResult mdd_closedqn(double[] mu, Matrix P, double[] servers, int N,
                                                 MDD reuse) {
        final int M = mu.length;
        if (P.getNumRows() != M || P.getNumCols() != M)
            throw new IllegalArgumentException("mdd_closedqn: routing matrix must be M x M");
        if (servers.length != M)
            throw new IllegalArgumentException("mdd_closedqn: servers must have one entry per station");
        if (N < 1)
            throw new IllegalArgumentException("mdd_closedqn: the closed population must be positive");

        // events: a completion at station i (rate mu_i * min(n_i, c_i)) routes
        // to station j with probability P(i,j); self-routing i == j leaves the
        // occupancy vector unchanged and is skipped.
        List<int[]> ev = new ArrayList<int[]>();
        List<Double> pr = new ArrayList<Double>();
        for (int i = 0; i < M; i++)
            for (int j = 0; j < M; j++)
                if (i != j && P.get(i, j) != 0.0) {
                    ev.add(new int[]{i, j});
                    pr.add(P.get(i, j));
                }
        final int E = ev.size();
        final int[] ii = new int[E];
        final int[] jj = new int[E];
        final double[] pij = new double[E];
        for (int a = 0; a < E; a++) {
            ii[a] = ev.get(a)[0];
            jj[a] = ev.get(a)[1];
            pij[a] = pr.get(a);
        }

        // all jobs start at station 1; an irreducible routing chain makes every
        // composition of N over M stations reachable.
        int[] domain = new int[M];
        int[] init = new int[M];
        for (int i = 0; i < M; i++) domain[i] = N + 1;
        init[0] = N;
        MddNextState nextfun = new MddNextState() {
            @Override
            public int[][] next(int[] s) {
                List<int[]> succ = new ArrayList<int[]>();
                for (int a = 0; a < E; a++) {
                    if (s[ii[a]] > 0) {
                        int[] t = s.clone();
                        t[ii[a]]--;
                        t[jj[a]]++;
                        succ.add(t);
                    }
                }
                return succ.toArray(new int[succ.size()][]);
            }
        };

        long t0 = System.nanoTime();
        MDD mdd;
        double timeReach;
        if (reuse != null) {
            mdd = reuse;
            timeReach = 0.0;
        } else {
            mdd = Mdd_reachset.mdd_reachset(domain, init, nextfun);
            timeReach = (System.nanoTime() - t0) / 1e9;
        }

        t0 = System.nanoTime();
        final long n = mdd.cardinality();
        final int ni = (int) n;
        int[][] states = mdd.enumerate();

        // assemble the generator directly from the MDD state indexing
        Matrix Q = new Matrix(ni, ni);
        for (int s = 0; s < ni; s++) {
            int[] st = states[s];
            int row = (int) mdd.index(st);
            for (int a = 0; a < E; a++) {
                int i = ii[a];
                if (st[i] > 0) {
                    int busy = Double.isInfinite(servers[i]) ? st[i]
                            : Math.min(st[i], (int) servers[i]);
                    double rate = mu[i] * busy * pij[a];
                    int[] t = st.clone();
                    t[i]--;
                    t[jj[a]]++;
                    int col = (int) mdd.index(t);
                    Q.set(row, col, Q.get(row, col) + rate);
                }
            }
        }
        Q = Ctmc_makeinfgen.ctmc_makeinfgen(Q);
        double timeGen = (System.nanoTime() - t0) / 1e9;

        t0 = System.nanoTime();
        Matrix pi = Ctmc_solve.ctmc_solve(Q);
        double timeSolve = (System.nanoTime() - t0) / 1e9;

        // performance metrics
        t0 = System.nanoTime();
        double[] QLen = new double[M];
        double[] U = new double[M];
        double[] X = new double[M];
        for (int i = 0; i < M; i++) {
            double qi = 0.0, busyMean = 0.0;
            for (int s = 0; s < ni; s++) {
                int cnt = states[s][i];
                int busy = Double.isInfinite(servers[i]) ? cnt : Math.min(cnt, (int) servers[i]);
                double p = pi.get(s);
                qi += p * cnt;
                busyMean += p * busy;
            }
            QLen[i] = qi;
            X[i] = mu[i] * busyMean;  // throughput = mean completion rate
            if (Double.isInfinite(servers[i]))
                U[i] = qi;  // mean number busy (IS station)
            else
                U[i] = busyMean / servers[i];
        }
        double timeMetrics = (System.nanoTime() - t0) / 1e9;

        return new MddClosedQnResult(mdd, Q, pi, states, QLen, U, X, mdd.stats(), timeReach,
                timeGen, timeSolve, timeMetrics);
    }
}
