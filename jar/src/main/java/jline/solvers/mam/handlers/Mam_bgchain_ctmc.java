/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_solve;
import jline.io.InputOutput;
import jline.lang.constant.SchedStrategy;
import jline.lang.state.State;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Background modulating chain of a mixed network: the continuous-time Markov
 * chain of the closed-class population vector, with the open classes present
 * only through the server capacity they leave free.
 *
 * <p>The chain carries B = 1 or 2 background classes. B = 1 is the single closed
 * chain of the model; B = 2 is the tagged/aggregate pair built by
 * {@link Solver_mam_bgchain} when the model has several closed chains (class 1
 * is the tagged chain, class 2 the flow-equivalent aggregate of the rest).</p>
 *
 * <p>A station holding n[i][1] + n[i][2] = e closed jobs serves background class
 * b at rate n[i][b]/ST[i][b] when it is an infinite server, and
 * cshare[i][e] * (n[i][b]/e) / ST[i][b] otherwise, splitting the capacity the
 * closed jobs hold over the background classes in proportion to their counts.
 * That is exact under PS and is the random-order surrogate under FCFS.</p>
 *
 * <p>The open classes enter ONLY through cshare, which is what makes this a
 * MODULATING chain rather than a joint model. The exchanged quantity is the
 * SHARE, already averaged over the open occupancy, and not the mean open
 * occupancy itself: e/(e+k) is convex in k, so rebuilding the share from a mean
 * k would bias the closed service rate downwards by Jensen's inequality, and the
 * closed throughput with it.</p>
 *
 * <p>The block of one background class is {@link State#spaceClosedSinglePublic},
 * the lattice primitive the CTMC solver enumerates a closed population over, and
 * the joint space is their cartesian product in class-major order. Only the
 * SUPPORT differs: the columns are the class's own stations rather than every
 * station. Keeping the reference primitive is what stops the enumeration, the
 * size counter in mam_bgchain_states and the CTMC solver from drifting apart.</p>
 *
 * @see Solver_mam_bgchain
 * @see Mam_bgchain_env
 */
public final class Mam_bgchain_ctmc {
    private Mam_bgchain_ctmc() {}

    private static final String MFILENAME = "mam_bgchain_ctmc";

    /** Solved background chain. */
    public static final class Result {
        /** space[s][i][b]: class-b jobs held by station i in state s. */
        public int[][][] space;
        /** totocc[s][i]: closed jobs of any background class held by station i. */
        public int[][] totocc;
        /** Stationary distribution, one row. */
        public Matrix pi;
        /** Generator. */
        public Matrix Q;
        /** Mean queue length per station per background class. */
        public double[][] QLen;
        /** Mean completion rate per station per background class. */
        public double[][] Tput;
        /** Mean fraction of the servers held, per station per background class. */
        public double[][] Ubusy;
        public int nstates;
    }

    /**
     * Builds and solves the background chain.
     *
     * @param Nb       population of each background class, length B in {1,2}
     * @param STb      STb[i][b], mean service time per station per background class
     * @param Pb       row-stochastic (Mc x Mc) routing matrix per background class
     * @param sched    discipline of each station of the chain's support
     * @param nservers servers of each station of the chain's support
     * @param cshare   cshare[i][e]: mean number of servers of station i that its
     *                 e closed jobs hold once the open work has taken its share
     * @param supp     supp[i][b]: station i is on the route of background class b.
     *                 NOT an optimization -- a chain that never visits a station
     *                 cannot hold jobs there, and enumerating the union of every
     *                 chain's stations puts probability on unreachable
     *                 configurations that also ABSORB, because the chain's routing
     *                 matrix has a zero row at an unvisited station which
     *                 row-normalizes to a self-loop. The generator turns reducible
     *                 and population conservation silently fails
     * @param options  solver options; reads config.bgstates_max
     * @return the solved chain
     */
    public static Result mam_bgchain_ctmc(int[] Nb, double[][] STb, List<Matrix> Pb,
                                          SchedStrategy[] sched, double[] nservers,
                                          double[][] cshare, boolean[][] supp,
                                          SolverOptions options) {
        int Mc = STb.length;
        int B = Nb.length;

        int bgstatesMax = 20000;
        Object cfg = options.config.get("bgstates_max");
        if (cfg instanceof Number) {
            bgstatesMax = ((Number) cfg).intValue();
        }

        // Each class is enumerated over ITS OWN stations only; see the supp javadoc.
        List<int[][]> sp = new ArrayList<int[][]>();      // expanded to Mc columns
        List<int[][]> compb = new ArrayList<int[][]>();   // reduced to the class's stations
        List<int[]> idxb = new ArrayList<int[]>();
        int[] nst = new int[B];
        long nstatesLong = 1;
        for (int b = 0; b < B; b++) {
            List<Integer> idxList = new ArrayList<Integer>();
            for (int i = 0; i < Mc; i++) {
                if (supp == null || supp[i][b]) idxList.add(Integer.valueOf(i));
            }
            if (idxList.isEmpty()) {
                if (Nb[b] > 0) {
                    InputOutput.line_error(MFILENAME, "Background class " + (b + 1) + " holds "
                            + Nb[b] + " jobs but visits no station.");
                }
                idxList.add(Integer.valueOf(0));  // an empty class still needs one slot
            }
            int[] idx = new int[idxList.size()];
            for (int t = 0; t < idx.length; t++) idx[t] = idxList.get(t).intValue();
            int[][] comp = closedBlock(idx.length, Nb[b]);
            int[][] exp = new int[comp.length][Mc];
            for (int r = 0; r < comp.length; r++) {
                for (int t = 0; t < idx.length; t++) exp[r][idx[t]] = comp[r][t];
            }
            idxb.add(idx);
            compb.add(comp);
            sp.add(exp);
            nst[b] = comp.length;
            nstatesLong *= nst[b];
        }
        if (nstatesLong > bgstatesMax) {
            InputOutput.line_error(MFILENAME, "The background chain of this model has " + nstatesLong
                    + " states, above the limit of " + bgstatesMax + ". The chain enumerates the closed-class "
                    + "population vector over the " + Mc + " stations the closed classes visit, so its size grows "
                    + "as nchoosek(N+Mc-1,Mc-1) per class, and it carries " + B + " classes. Lower "
                    + "options.config.bgaggr to aggregate more of the closed chains into fewer classes, raise "
                    + "options.config.bgstates_max to solve it anyway, or reduce the closed populations.");
        }
        int nstates = (int) nstatesLong;

        // Per-class transition targets: tgt[b][s][i*Mc+j] is the class-b state
        // reached from s when one job moves from station i to station j, or -1.
        List<int[][]> tgt = new ArrayList<int[][]>();
        for (int b = 0; b < B; b++) {
            int[][] comp = compb.get(b);
            int[] idx = idxb.get(b);
            int mb = idx.length;
            long[] keys = new long[comp.length];
            for (int s = 0; s < comp.length; s++) {
                keys[s] = encode(comp[s], Nb[b] + 1);
            }
            long[] sorted = keys.clone();
            Arrays.sort(sorted);
            int[] ord = argsortByKey(keys);
            int[][] t = new int[comp.length][Mc * Mc];
            for (int s = 0; s < comp.length; s++) {
                Arrays.fill(t[s], -1);
            }
            int[] work = new int[mb];
            for (int s = 0; s < comp.length; s++) {
                for (int ii = 0; ii < mb; ii++) {
                    if (comp[s][ii] == 0) continue;
                    for (int jj = 0; jj < mb; jj++) {
                        if (jj == ii) continue;
                        System.arraycopy(comp[s], 0, work, 0, mb);
                        work[ii]--;
                        work[jj]++;
                        int pos = binarySearch(sorted, encode(work, Nb[b] + 1));
                        if (pos >= 0) {
                            t[s][idx[ii] * Mc + idx[jj]] = ord[pos];
                        }
                    }
                }
            }
            tgt.add(t);
        }

        // Joint space, class 1 outermost
        int[] strideb = new int[B];
        for (int b = 0; b < B; b++) {
            int s = 1;
            for (int b2 = b + 1; b2 < B; b2++) s *= nst[b2];
            strideb[b] = s;
        }
        int[][][] space = new int[nstates][Mc][B];
        int[][] subidx = new int[nstates][B];
        int[][] totocc = new int[nstates][Mc];
        for (int s = 0; s < nstates; s++) {
            int rem = s;
            for (int b = B - 1; b >= 0; b--) {
                subidx[s][b] = rem % nst[b];
                rem = rem / nst[b];
            }
            for (int b = 0; b < B; b++) {
                int[] vec = sp.get(b)[subidx[s][b]];
                for (int i = 0; i < Mc; i++) {
                    space[s][i][b] = vec[i];
                    totocc[s][i] += vec[i];
                }
            }
        }

        double[][] mu = new double[Mc][B];
        for (int b = 0; b < B; b++) {
            for (int i = 0; i < Mc; i++) {
                if (STb[i][b] > 0 && Double.isFinite(STb[i][b])) {
                    mu[i][b] = 1.0 / STb[i][b];
                }
            }
        }

        Matrix Q = new Matrix(nstates, nstates);
        double[][][] rateFull = new double[nstates][Mc][B];
        double[][] capBusy = new double[nstates][Mc];
        for (int s = 0; s < nstates; s++) {
            for (int i = 0; i < Mc; i++) {
                int eclosed = totocc[s][i];
                if (eclosed == 0) continue;
                double held;
                if (sched[i] == SchedStrategy.INF) {
                    held = eclosed;
                } else {
                    int ecap = Math.min(eclosed, cshare[i].length - 1);
                    held = cshare[i][ecap];
                }
                if (held <= 0) continue;
                capBusy[s][i] = held;
                for (int b = 0; b < B; b++) {
                    if (space[s][i][b] == 0 || mu[i][b] == 0) continue;
                    double r = held * ((double) space[s][i][b] / eclosed) * mu[i][b];
                    rateFull[s][i][b] = r;
                    Matrix P = Pb.get(b);
                    for (int j = 0; j < Mc; j++) {
                        if (j == i || P.get(i, j) <= 0) continue;
                        int tsub = tgt.get(b)[subidx[s][b]][i * Mc + j];
                        if (tsub < 0) continue;
                        int sdest = s + (tsub - subidx[s][b]) * strideb[b];
                        Q.set(s, sdest, Q.get(s, sdest) + r * P.get(i, j));
                    }
                }
            }
        }

        Q = Ctmc_makeinfgen.ctmc_makeinfgen(Q);
        Matrix pi;
        if (nstates == 1) {
            pi = new Matrix(1, 1);
            pi.set(0, 0, 1.0);
        } else {
            pi = Ctmc_solve.ctmc_solve(Q);
        }
        double tot = 0;
        for (int s = 0; s < nstates; s++) {
            double v = Math.max(pi.get(0, s), 0.0);
            pi.set(0, s, v);
            tot += v;
        }
        if (tot > 0) pi.scaleEq(1.0 / tot);

        Result res = new Result();
        res.space = space;
        res.totocc = totocc;
        res.pi = pi;
        res.Q = Q;
        res.nstates = nstates;
        res.QLen = new double[Mc][B];
        res.Tput = new double[Mc][B];
        res.Ubusy = new double[Mc][B];
        for (int s = 0; s < nstates; s++) {
            double p = pi.get(0, s);
            if (p == 0) continue;
            for (int i = 0; i < Mc; i++) {
                for (int b = 0; b < B; b++) {
                    res.QLen[i][b] += p * space[s][i][b];
                    res.Tput[i][b] += p * rateFull[s][i][b];
                }
            }
        }
        for (int i = 0; i < Mc; i++) {
            if (sched[i] == SchedStrategy.INF) {
                for (int b = 0; b < B; b++) res.Ubusy[i][b] = res.QLen[i][b];
            } else {
                for (int s = 0; s < nstates; s++) {
                    double p = pi.get(0, s);
                    int occ = totocc[s][i];
                    if (p == 0 || occ == 0) continue;
                    for (int b = 0; b < B; b++) {
                        res.Ubusy[i][b] += p * capBusy[s][i] * ((double) space[s][i][b] / occ) / nservers[i];
                    }
                }
            }
        }
        return res;
    }

    /**
     * Block of one background class: the ways to place N jobs over the m stations
     * that class visits. This is {@link State#spaceClosedSinglePublic}, the lattice
     * primitive the CTMC solver enumerates a closed population over, so the row
     * order and the row count are the reference's rather than this class's; only
     * the support differs, being the class's own stations rather than all of them.
     */
    static int[][] closedBlock(int m, int N) {
        Matrix ss = State.spaceClosedSinglePublic(m, N);
        int[][] comp = new int[ss.getNumRows()][m];
        for (int r = 0; r < ss.getNumRows(); r++) {
            for (int t = 0; t < m; t++) {
                comp[r][t] = (int) Math.round(ss.get(r, t));
            }
        }
        return comp;
    }

    private static long encode(int[] vec, int base) {
        long key = 0;
        long mult = 1;
        for (int i = 0; i < vec.length; i++) {
            key += vec[i] * mult;
            mult *= base;
        }
        return key;
    }

    private static int[] argsortByKey(long[] keys) {
        int n = keys.length;
        Integer[] idx = new Integer[n];
        for (int i = 0; i < n; i++) idx[i] = Integer.valueOf(i);
        final long[] k = keys;
        Arrays.sort(idx, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Long.compare(k[a.intValue()], k[b.intValue()]);
            }
        });
        int[] ord = new int[n];
        for (int i = 0; i < n; i++) ord[i] = idx[i].intValue();
        return ord;
    }

    private static int binarySearch(long[] sorted, long key) {
        int lo = 0;
        int hi = sorted.length - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid] == key) return mid;
            if (sorted[mid] < key) lo = mid + 1; else hi = mid - 1;
        }
        return -1;
    }
}
