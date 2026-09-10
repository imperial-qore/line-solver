/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba.analyzers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.api.sn.SnRtStations;
import jline.api.snc.Snc_env_map;
import jline.api.snc.Snc_env_poisson;
import jline.api.snc.Snc_leftover;
import jline.api.snc.Snc_mean_delay;
import jline.api.snc.Snc_output;
import jline.api.snc.Snc_srv_exp;
import jline.api.snc.SncEnvelope;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Stochastic network calculus UPPER bound on the mean response times and queue
 * lengths of a feed-forward open network, valid for EVERY work-conserving
 * scheduling policy at every station.
 *
 * <p>Port of matlab/src/solvers/BA/solver_ba_snc_analyzer.m. The api domain
 * {@link jline.api.snc} supplies the envelope algebra; this analyzer maps the
 * LINE model onto it, propagates envelopes hop by hop, and reads the bound back
 * per station and class.</p>
 *
 * <p>UNITS ARE JOBS, NOT WORK. The arrival envelope counts jobs and the service
 * element is {@link Snc_srv_exp}, the counting process of an Exp(mu) server.
 * That is what lets a departure envelope from one station be the arrival
 * envelope of the next: a service-time work unit differs from station to
 * station, a job does not. On a single M/M/1 the resulting backlog bound decays
 * as (lambda/mu)^n and the delay bound as exp(-(mu-lambda)*d), both exact
 * rates.</p>
 *
 * <p>BOUND CONVENTION. R(i,r) is {@link Snc_mean_delay} of the (arrival,
 * service) envelope pair at that station, i.e. the integral of the delay tail
 * bound, so each entry is a valid upper bound on its own. Q follows by Little's
 * law from the bounded R and the EXACT throughput T (an open network's
 * per-class rates are fixed by the traffic equations, not by the policy), and
 * so does C. U is exact for the same reason.</p>
 *
 * <p>WHAT IS ASSUMED, and refused when it does not hold: a fully open network
 * with a Source, no delay station, one server per station, exponential service
 * (the Source may be any Markovian (D0,D1) process), a FEED-FORWARD station
 * graph, DETERMINISTIC routing downstream of the Source (splitting AT a Poisson
 * Source is exact and is allowed), and EQUAL service rates among the classes
 * sharing a station.</p>
 *
 * <p>TIGHTNESS. This is a policy-robust bound, so it is loose on the mean: 2.4x
 * the exact M/M/1 mean response time at rho = 0.1 and 10.4x at rho = 0.95. Its
 * sharp object is the TAIL, whose decay rate it reproduces exactly; reach it
 * through the SolverBA quantile accessors rather than through the mean
 * columns.</p>
 *
 * <p>Reference: M. Fidler, A. Rizk (2015). A Guide to the Stochastic Network
 * Calculus. IEEE Communications Surveys and Tutorials 17(1), 92-105.</p>
 */
public final class Solver_ba_snc_analyzer {

    private Solver_ba_snc_analyzer() {}

    private static final double TOL = 1e-10;

    /**
     * The per-pair envelopes the analyzer built, kept so that the quantile
     * accessors can read the tail without a second pass over the network.
     */
    public static final class Envelopes {
        /** Arrival envelope of each (station, class) pair that carries traffic. */
        public final Map<Long, SncEnvelope> arv = new HashMap<Long, SncEnvelope>();
        /** Service envelope of the same pairs, cross traffic already removed. */
        public final Map<Long, SncEnvelope> srv = new HashMap<Long, SncEnvelope>();
        /** Exact arrival rate of each pair. */
        public final Map<Long, Double> lam = new HashMap<Long, Double>();
        /** Service rate of each pair. */
        public final Map<Long, Double> mu = new HashMap<Long, Double>();
        /** Number of stations of the model. */
        public int M;
        /** Number of classes of the model. */
        public int K;

        /** @param i station index @param r class index @return the map key */
        public long key(int i, int r) {
            return (long) i * 1000000L + r;
        }
    }

    /**
     * @param sn      the model
     * @param options solver options (unused: the family has no tunable)
     * @return the (Q,U,R,T,C,X) block of the 'snc.upper' family
     */
    public static MVAResult solver_ba_snc_analyzer(NetworkStruct sn, SolverOptions options) {
        final int M = sn.nstations;
        final int K = sn.nclasses;
        Envelopes env = envelopes(sn);

        MVAResult ret = new MVAResult();
        ret.QN = new Matrix(M, K);
        ret.UN = new Matrix(M, K);
        ret.RN = new Matrix(M, K);
        ret.TN = new Matrix(M, K);
        ret.CN = new Matrix(1, K);
        ret.XN = new Matrix(1, K);
        ret.logNormConstAggr = Double.NaN;
        ret.iter = 1;

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                long k = env.key(i, r);
                if (!env.arv.containsKey(k)) {
                    continue;
                }
                ret.RN.set(i, r, Snc_mean_delay.snc_mean_delay(env.arv.get(k), env.srv.get(k)).value);
                ret.TN.set(i, r, env.lam.get(k));
                ret.UN.set(i, r, env.lam.get(k) / env.mu.get(k));
            }
        }

        // ---- exact open-network quantities ----
        for (int i = 0; i < M; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) != NodeType.Source) {
                continue;
            }
            for (int r = 0; r < K; r++) {
                double arr = sn.rates.get(i, r);
                if (Double.isFinite(arr) && arr > 0) {
                    ret.TN.set(i, r, ret.TN.get(i, r) + arr);
                    ret.XN.set(0, r, ret.XN.get(0, r) + arr);
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                ret.QN.set(i, r, ret.TN.get(i, r) * ret.RN.get(i, r));
            }
        }
        for (int r = 0; r < K; r++) {
            if (ret.XN.get(0, r) > 0) {
                double sum = 0;
                for (int i = 0; i < M; i++) {
                    sum += ret.QN.get(i, r);
                }
                ret.CN.set(0, r, sum / ret.XN.get(0, r));
            }
        }
        return ret;
    }

    /**
     * Builds the per-pair (arrival, service) envelopes of a feed-forward model.
     *
     * <p>Every gate of the family is checked here rather than in the caller, so
     * that the quantile accessors refuse an unsupported model with the same
     * reason as the mean columns.</p>
     *
     * @param sn the model
     * @return the envelopes of the pairs that carry traffic
     */
    public static Envelopes envelopes(NetworkStruct sn) {
        final int M = sn.nstations;
        final int K = sn.nclasses;
        Envelopes env = new Envelopes();
        env.M = M;
        env.K = K;

        // ---- model gates ----
        for (int r = 0; r < sn.njobs.length(); r++) {
            if (Double.isFinite(sn.njobs.get(r))) {
                throw new RuntimeException(
                        "Method 'snc.upper' supports fully open networks only (no closed classes).");
            }
        }
        boolean[] isSource = new boolean[M];
        List<Integer> srcList = new ArrayList<Integer>();
        List<Integer> qstat = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            isSource[i] = sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source;
            if (isSource[i]) {
                srcList.add(i);
            } else {
                qstat.add(i);
            }
        }
        if (srcList.isEmpty()) {
            throw new RuntimeException(
                    "Method 'snc.upper' requires an open network with a Source station.");
        }
        for (int a = 0; a < qstat.size(); a++) {
            int i = qstat.get(a);
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                throw new RuntimeException("Method 'snc.upper' does not support delay "
                        + "(infinite-server) stations: the service envelope is that of a single "
                        + "busy server.");
            }
            if (sn.nservers.get(i, 0) > 1) {
                throw new RuntimeException(
                        "Method 'snc.upper' does not support multi-server stations.");
            }
        }

        // ---- station-space routing, with the Source absorbed into the injections ----
        Pair<Matrix, Matrix> rtstPair = SnRtStations.snRtStations(sn);
        Matrix rtst = rtstPair.getLeft();

        int nqK = qstat.size() * K;
        int[] pairStation = new int[nqK];
        int[] pairClass = new int[nqK];
        int[] pairFlat = new int[nqK];
        int np = 0;
        for (int a = 0; a < qstat.size(); a++) {
            int i = qstat.get(a);
            for (int r = 0; r < K; r++) {
                pairStation[np] = i;
                pairClass[np] = r;
                pairFlat[np] = i * K + r;
                np++;
            }
        }

        // Injections are kept per (source, class) so that the exogenous process
        // of each stream is still identifiable once the pairs are known.
        int ncols = srcList.size() * K;
        int[] srcOfCol = new int[ncols];
        int[] clsOfCol = new int[ncols];
        double[][] inject = new double[np][ncols];
        double[] lambda0 = new double[np];
        int col = 0;
        for (int si = 0; si < srcList.size(); si++) {
            int s = srcList.get(si);
            for (int r0 = 0; r0 < K; r0++) {
                srcOfCol[col] = s;
                clsOfCol[col] = r0;
                double arr = sn.rates.get(s, r0);
                if (Double.isFinite(arr) && arr > 0) {
                    int srow = s * K + r0;
                    for (int p = 0; p < np; p++) {
                        inject[p][col] = arr * rtst.get(srow, pairFlat[p]);
                        lambda0[p] += inject[p][col];
                    }
                }
                col++;
            }
        }

        Matrix Pfull = new Matrix(np, np);
        for (int p = 0; p < np; p++) {
            for (int q = 0; q < np; q++) {
                double v = rtst.get(pairFlat[p], pairFlat[q]);
                if (v != 0) {
                    Pfull.set(p, q, v);
                }
            }
        }

        // ---- restrict to the pairs that actually carry traffic ----
        Matrix ImPt = new Matrix(np, np);
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < np; j++) {
                ImPt.set(i, j, (i == j ? 1.0 : 0.0) - Pfull.get(j, i));
            }
        }
        Matrix rhs0 = new Matrix(np, 1);
        for (int p = 0; p < np; p++) {
            rhs0.set(p, 0, lambda0[p]);
        }
        Matrix lamAll = ImPt.inv().mult(rhs0);
        double lamMax = 0;
        for (int p = 0; p < np; p++) {
            lamMax = Math.max(lamMax, lamAll.get(p, 0));
        }
        double lamTol = 1e-12 * Math.max(1.0, lamMax);
        List<Integer> keep = new ArrayList<Integer>();
        for (int p = 0; p < np; p++) {
            if (lamAll.get(p, 0) > lamTol) {
                keep.add(p);
            }
        }
        if (keep.isEmpty()) {
            throw new RuntimeException("The model carries no open traffic.");
        }

        final int nk = keep.size();
        double[] lam = new double[nk];
        double[] mu = new double[nk];
        int[] statk = new int[nk];
        int[] clsk = new int[nk];
        double[][] injk = new double[nk][ncols];
        Matrix Pk = new Matrix(nk, nk);
        for (int a = 0; a < nk; a++) {
            int p = keep.get(a);
            lam[a] = lamAll.get(p, 0);
            statk[a] = pairStation[p];
            clsk[a] = pairClass[p];
            for (int c = 0; c < ncols; c++) {
                injk[a][c] = inject[p][c];
            }
            for (int b = 0; b < nk; b++) {
                double v = Pfull.get(p, keep.get(b));
                if (v != 0) {
                    Pk.set(a, b, v);
                }
            }
            mu[a] = sn.rates.get(statk[a], clsk[a]);
            if (!Double.isFinite(mu[a]) || mu[a] <= 0) {
                throw new RuntimeException("Station " + (statk[a] + 1) + " has no service rate for "
                        + "class " + (clsk[a] + 1) + " but carries its traffic.");
            }
            ProcessType pt = sn.procid.get(sn.stations.get(statk[a])).get(sn.jobclasses.get(clsk[a]));
            if (pt != ProcessType.EXP) {
                throw new RuntimeException("Method 'snc.upper' requires exponential service: "
                        + "station " + (statk[a] + 1) + " class " + (clsk[a] + 1) + " is " + pt + ".");
            }
        }

        // ---- routing restrictions: no split downstream of the Source ----
        for (int a = 0; a < nk; a++) {
            int nsucc = 0;
            int succ = -1;
            for (int b = 0; b < nk; b++) {
                if (Pk.get(a, b) > TOL) {
                    nsucc++;
                    succ = b;
                }
            }
            if (nsucc > 1) {
                throw new RuntimeException("Method 'snc.upper' requires deterministic routing "
                        + "downstream of the Source: station " + (statk[a] + 1) + " class "
                        + (clsk[a] + 1) + " splits its flow over " + nsucc + " destinations.");
            }
            if (nsucc == 1 && Math.abs(Pk.get(a, succ) - 1.0) > 1e-8) {
                throw new RuntimeException("Method 'snc.upper' requires deterministic routing "
                        + "downstream of the Source: station " + (statk[a] + 1) + " class "
                        + (clsk[a] + 1) + " routes onward with probability " + Pk.get(a, succ) + ".");
            }
        }
        // A Source that splits is exact only when its process is Poisson, since
        // a Bernoulli thinning of a Poisson stream is again Poisson.
        for (int c = 0; c < ncols; c++) {
            int ndest = 0;
            for (int a = 0; a < nk; a++) {
                if (injk[a][c] > TOL) {
                    ndest++;
                }
            }
            if (ndest <= 1) {
                continue;
            }
            int s = srcOfCol[c];
            int r0 = clsOfCol[c];
            ProcessType pt = sn.procid.get(sn.stations.get(s)).get(sn.jobclasses.get(r0));
            if (pt != ProcessType.EXP) {
                throw new RuntimeException("Method 'snc.upper' can split only a Poisson Source: "
                        + "source " + (s + 1) + " class " + (r0 + 1) + " is " + pt + " and feeds "
                        + ndest + " stations.");
            }
        }

        // ---- one service rate per station, and a feed-forward station graph ----
        List<Integer> stationsUsed = new ArrayList<Integer>();
        for (int a = 0; a < nk; a++) {
            if (!stationsUsed.contains(statk[a])) {
                stationsUsed.add(statk[a]);
            }
        }
        for (int u = 0; u < stationsUsed.size(); u++) {
            int i = stationsUsed.get(u);
            double lo = Double.POSITIVE_INFINITY;
            double hi = 0;
            for (int a = 0; a < nk; a++) {
                if (statk[a] == i) {
                    lo = Math.min(lo, mu[a]);
                    hi = Math.max(hi, mu[a]);
                }
            }
            if (hi - lo > 1e-8 * Math.max(1.0, hi)) {
                throw new RuntimeException("Method 'snc.upper' requires the classes sharing a "
                        + "station to have equal service rates: station " + (i + 1)
                        + " carries rates in [" + lo + ", " + hi + "].");
            }
        }

        int ns = stationsUsed.size();
        boolean[][] adj = new boolean[ns][ns];
        for (int a = 0; a < nk; a++) {
            for (int b = 0; b < nk; b++) {
                if (Pk.get(a, b) > TOL) {
                    adj[stationsUsed.indexOf(statk[a])][stationsUsed.indexOf(statk[b])] = true;
                }
            }
        }
        int[] order = topoOrder(adj);
        if (order == null) {
            throw new RuntimeException("Method 'snc.upper' requires a feed-forward network: the "
                    + "station graph has a cycle, so a station's cross traffic is not determined "
                    + "upstream of it.");
        }

        // ---- envelope propagation, station by station in feed-forward order ----
        SncEnvelope[] arvH = new SncEnvelope[nk];
        SncEnvelope[] srvH = new SncEnvelope[nk];
        final SncEnvelope[] outH = new SncEnvelope[nk];
        for (int oi = 0; oi < ns; oi++) {
            int i = stationsUsed.get(order[oi]);
            List<Integer> here = new ArrayList<Integer>();
            for (int a = 0; a < nk; a++) {
                if (statk[a] == i) {
                    here.add(a);
                }
            }
            for (int h = 0; h < here.size(); h++) {
                final int a = here.get(h);
                List<SncEnvelope> parts = new ArrayList<SncEnvelope>();
                for (int c = 0; c < ncols; c++) {
                    if (injk[a][c] > TOL) {
                        int s = srcOfCol[c];
                        int r0 = clsOfCol[c];
                        ProcessType pt =
                                sn.procid.get(sn.stations.get(s)).get(sn.jobclasses.get(r0));
                        if (pt == ProcessType.EXP) {
                            parts.add(Snc_env_poisson.of(injk[a][c]));
                        } else {
                            MatrixCell proc = sn.proc.get(sn.stations.get(s))
                                    .get(sn.jobclasses.get(r0));
                            parts.add(Snc_env_map.of(proc.get(0), proc.get(1)));
                        }
                    }
                }
                for (int b = 0; b < nk; b++) {
                    if (Pk.get(b, a) > TOL) {
                        parts.add(indirect(outH, b));
                    }
                }
                if (parts.isEmpty()) {
                    throw new RuntimeException("Station " + (i + 1) + " class " + (clsk[a] + 1)
                            + " carries traffic with no identifiable source.");
                }
                arvH[a] = sum(parts);
            }
            for (int h = 0; h < here.size(); h++) {
                int a = here.get(h);
                List<SncEnvelope> cross = new ArrayList<SncEnvelope>();
                for (int g = 0; g < here.size(); g++) {
                    if (here.get(g) != a) {
                        cross.add(arvH[here.get(g)]);
                    }
                }
                srvH[a] = leftover(mu[a], cross);
                outH[a] = Snc_output.of(arvH[a], srvH[a]);
            }
        }

        for (int a = 0; a < nk; a++) {
            long k = env.key(statk[a], clsk[a]);
            env.arv.put(k, arvH[a]);
            env.srv.put(k, srvH[a]);
            env.lam.put(k, lam[a]);
            env.mu.put(k, mu[a]);
        }
        return env;
    }

    /**
     * Superposition of independent flows: the exponential forms multiply, so
     * the envelope terms add.
     */
    private static SncEnvelope sum(final List<SncEnvelope> parts) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                double sigma = 0;
                double rho = 0;
                for (int k = 0; k < parts.size(); k++) {
                    double[] e = parts.get(k).eval(theta);
                    sigma += e[0];
                    rho += e[1];
                }
                return new double[] {sigma, rho};
            }
        };
    }

    /** The busy-server counting process minus what blind multiplexing gives away. */
    private static SncEnvelope leftover(final double mu, final List<SncEnvelope> cross) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                double[] s = Snc_srv_exp.snc_srv_exp(mu, theta);
                if (cross.isEmpty()) {
                    return s;
                }
                double[] x = sum(cross).eval(theta);
                return Snc_leftover.snc_leftover(s[0], s[1], x[0], x[1]);
            }
        };
    }

    /**
     * A departure envelope read through its slot rather than by value: the
     * predecessor's entry is filled later in the same station pass, so
     * capturing it eagerly would capture a null.
     */
    private static SncEnvelope indirect(final SncEnvelope[] slots, final int idx) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return slots[idx].eval(theta);
            }
        };
    }

    /**
     * Kahn's algorithm. Returns null when the graph has a cycle, which is the
     * feed-forward refusal.
     */
    private static int[] topoOrder(boolean[][] adj) {
        int n = adj.length;
        int[] indeg = new int[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (adj[i][j]) {
                    indeg[j]++;
                }
            }
        }
        boolean[] done = new boolean[n];
        int[] order = new int[n];
        int placed = 0;
        while (true) {
            int cand = -1;
            for (int i = 0; i < n && cand < 0; i++) {
                if (!done[i] && indeg[i] == 0) {
                    cand = i;
                }
            }
            if (cand < 0) {
                break;
            }
            order[placed++] = cand;
            done[cand] = true;
            for (int j = 0; j < n; j++) {
                if (adj[cand][j]) {
                    indeg[j]--;
                }
            }
            indeg[cand] = 1; // keep it out of the candidate set
        }
        return placed < n ? null : order;
    }
}
