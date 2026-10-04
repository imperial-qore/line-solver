/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import jline.GlobalConstants;
import jline.api.dpfqn.Dpfqn_nc;
import jline.api.dpfqn.Dpfqn_ncld;
import jline.api.dpfqn.DpfqnNcLdResult;
import jline.api.dpfqn.DpfqnNcResult;
import jline.api.dqsys.Bernoulli1Result;
import jline.api.dqsys.Dqsys_bernoulli1;
import jline.api.sn.SnRtStations;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Exact normalizing-constant analysis of a discrete-time (slotted) queueing
 * model, selected by {@code options.config.slotted}.
 *
 * <p>Two families are covered, both from H. Daduna, Queueing Networks with
 * Discrete Time Scale, LNCS 2046, Springer, 2001:</p>
 *
 * <ul>
 *   <li>chapter 2, a Bernoulli server fed by a Bernoulli arrival stream, with
 *       an unbounded buffer (theorem 2.3, corollary 2.7), a finite buffer
 *       (corollary 2.8) or a load-dependent service probability (example
 *       2.10), evaluated by {@link Dqsys_bernoulli1};</li>
 *   <li>chapter 3, a closed cycle of Bernoulli servers (theorem 3.2, corollary
 *       3.4), evaluated by {@link Dpfqn_nc} when the service probabilities are
 *       state independent and by {@link Dpfqn_ncld} otherwise.</li>
 * </ul>
 *
 * <p>Every metric is expressed on the slot lattice: a rate is a per-slot
 * probability and a time is a number of slots. {@code options.config.slotlength}
 * rescales both to model time units.</p>
 *
 * <p>On the cycle route the per-class split is proportional to the per-class
 * population. Service in the cycle is type independent and FCFS forbids
 * overtaking, so the cyclic order of the jobs is frozen; the marginal law of
 * the queue lengths carries no class information, and the long-run share of
 * station j held by chain g is its population share N_g/N. That is the sense in
 * which section 3.2 of the reference calls the multichain case a direct
 * adaptation of the unichain one.</p>
 */
public final class Solver_nc_dt_analyzer {

    private Solver_nc_dt_analyzer() {}

    /** Classification of a model against the discrete-time product form. */
    public static final class DtModel {
        /** One of "bernoulli1", "cycle" or "none". */
        public String kind = "none";
        /** Why the model is not admissible, when kind is "none". */
        public String reason = "";
        /** Queueing station index (bernoulli1). */
        public int station = -1;
        /** Source station index (bernoulli1). */
        public int source = -1;
        /** Offered arrival probability (bernoulli1). */
        public double arrivalProb;
        /** Service probabilities p(n), n = 1..L (bernoulli1). */
        public double[] serviceProbSingle;
        /** Buffer capacity, {@link Integer#MAX_VALUE} when unbounded. */
        public int capacity = Integer.MAX_VALUE;
        /** Station indices in cycle order (cycle). */
        public int[] order;
        /** Service probabilities p_j(n), indexed [cycle position][n-1]. */
        public double[][] serviceProb;
        /** Total closed population (cycle). */
        public int population;
    }

    /** True when the caller asked for the discrete-time route. */
    public static boolean isSlotted(SolverOptions options) {
        return options != null && options.config != null && options.config.slotted;
    }

    private static double slotLength(SolverOptions options) {
        double d = 1.0;
        if (options != null && options.config != null && options.config.slotlength != null) {
            d = options.config.slotlength;
        }
        if (Double.isNaN(d) || Double.isInfinite(d) || d <= 0) {
            throw new RuntimeException("options.config.slotlength must be a positive finite scalar");
        }
        return d;
    }

    /**
     * Classify sn against the two discrete-time product-form families.
     *
     * <p>The admissible feature set is narrow because the discrete-time product
     * form is narrow. Beyond the geometric service requirement a cycle is the
     * only topology (section 4.1 records that general discrete-time topologies
     * of FCFS Bernoulli servers have no product form), every station must be a
     * single server (Pestien and Ramakrishnan, quoted before example 2.10,
     * proved that a multiserver node inside a cycle of geometrical queues
     * destroys the product form for any finite server count), and class
     * switching is rejected.</p>
     */
    public static DtModel ncIsDtModel(NetworkStruct sn) {
        DtModel dt = new DtModel();
        if (sn == null) {
            dt.reason = "no network structure available";
            return dt;
        }
        int M = sn.nstations;
        int R = sn.nclasses;

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(i, r);
                if (!Double.isNaN(rate) && !Double.isInfinite(rate) && rate > 0) {
                    ProcessType pt = procType(sn, i, r);
                    if (pt != ProcessType.GEOMETRIC) {
                        dt.reason = String.format("station %d serves class %d with a %s process; "
                                + "a discrete-time model needs Geometric service and interarrival times",
                                i, r, pt);
                        return dt;
                    }
                }
            }
        }

        if (sn.csmask != null && sn.csmask.getNumRows() == R && sn.csmask.getNumCols() == R) {
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    if (r != s && sn.csmask.get(r, s) > 0) {
                        dt.reason = "class switching is not covered by the discrete-time product form";
                        return dt;
                    }
                }
            }
        }
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            dt.reason = "class-dependent scaling is not covered by the discrete-time product form";
            return dt;
        }
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            dt.reason = "joint-dependent scaling is not covered by the discrete-time product form";
            return dt;
        }

        boolean isOpen = false;
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                isOpen = true;
                break;
            }
        }
        return isOpen ? dtOpenSingle(sn, dt) : dtClosedCycle(sn, dt);
    }

    private static ProcessType procType(NetworkStruct sn, int i, int r) {
        if (sn.procid == null || sn.stations == null || i >= sn.stations.size()) {
            return null;
        }
        Station st = sn.stations.get(i);
        if (!sn.procid.containsKey(st)) {
            return null;
        }
        java.util.Map<jline.lang.JobClass, ProcessType> row = sn.procid.get(st);
        for (java.util.Map.Entry<jline.lang.JobClass, ProcessType> e : row.entrySet()) {
            if (e.getKey().getIndex() - 1 == r) {
                return e.getValue();
            }
        }
        return null;
    }

    private static DtModel dtOpenSingle(NetworkStruct sn, DtModel dt) {
        int M = sn.nstations;
        int R = sn.nclasses;
        if (R != 1) {
            dt.reason = "the discrete-time single-node route handles one open class";
            return dt;
        }
        int src = -1;
        int nsrc = 0;
        for (int i = 0; i < M; i++) {
            if (sched(sn, i) == SchedStrategy.EXT) {
                src = i;
                nsrc++;
            }
        }
        if (nsrc != 1) {
            dt.reason = "an open discrete-time model needs exactly one Source";
            return dt;
        }
        int ist = -1;
        int nq = 0;
        for (int i = 0; i < M; i++) {
            if (i != src) {
                ist = i;
                nq++;
            }
        }
        if (nq != 1) {
            dt.reason = String.format(
                    "the discrete-time single-node route handles one queueing station, found %d", nq);
            return dt;
        }
        if (sched(sn, ist) != SchedStrategy.FCFS) {
            dt.reason = "a Bernoulli server is a FCFS station";
            return dt;
        }
        double c = sn.nservers.get(ist, 0);
        if (!Double.isInfinite(c) && c != 1) {
            dt.reason = "a Bernoulli server is a single-server station; use load dependence "
                    + "for the multiserver approximation of example 2.10";
            return dt;
        }

        double b = sn.rates.get(src, 0);
        double p = sn.rates.get(ist, 0);
        if (Double.isNaN(b) || Double.isInfinite(b) || b <= 0 || b > 1) {
            dt.reason = "the source arrival probability must lie in (0,1]";
            return dt;
        }
        if (Double.isNaN(p) || Double.isInfinite(p) || p <= 0 || p > 1) {
            dt.reason = "the service probability must lie in (0,1]";
            return dt;
        }

        int cap = capacity(sn, ist, 0);
        double[] lld = lld(sn, ist);
        if (cap == Integer.MAX_VALUE && lld != null) {
            dt.reason = "a load-dependent Bernoulli server needs a finite capacity to bound the state space";
            return dt;
        }
        if (cap == Integer.MAX_VALUE && b >= p) {
            dt.reason = "an unbounded discrete-time queue needs an arrival probability below the service probability";
            return dt;
        }

        if (cap == Integer.MAX_VALUE) {
            dt.serviceProbSingle = new double[]{p};
        } else {
            double[] alpha = lldVector(lld, cap);
            dt.serviceProbSingle = new double[cap];
            for (int n = 0; n < cap; n++) {
                dt.serviceProbSingle[n] = p * alpha[n];
                if (dt.serviceProbSingle[n] <= 0 || dt.serviceProbSingle[n] > 1) {
                    dt.reason = "load dependence must keep the service probability inside (0,1]";
                    return dt;
                }
            }
        }

        dt.kind = "bernoulli1";
        dt.station = ist;
        dt.source = src;
        dt.arrivalProb = b;
        dt.capacity = cap;
        return dt;
    }

    private static DtModel dtClosedCycle(NetworkStruct sn, DtModel dt) {
        int M = sn.nstations;
        int R = sn.nclasses;
        double Ntot = 0;
        for (int r = 0; r < R; r++) {
            double v = sn.njobs.get(0, r);
            if (!Double.isInfinite(v)) {
                Ntot += v;
            }
        }
        if (Ntot <= 0 || Ntot != Math.floor(Ntot)) {
            dt.reason = "the closed population must be a positive integer";
            return dt;
        }
        int N = (int) Ntot;

        for (int i = 0; i < M; i++) {
            if (sched(sn, i) != SchedStrategy.FCFS) {
                dt.reason = String.format(
                        "station %d is not FCFS; a cycle of Bernoulli servers is FCFS throughout", i);
                return dt;
            }
            double c = sn.nservers.get(i, 0);
            if (!Double.isInfinite(c) && c != 1) {
                dt.reason = String.format("station %d has %d servers; a multiserver node inside a "
                        + "cycle of geometrical queues has no product form", i, (int) c);
                return dt;
            }
        }

        double[] p = new double[M];
        for (int i = 0; i < M; i++) {
            double lo = Double.POSITIVE_INFINITY;
            double hi = Double.NEGATIVE_INFINITY;
            for (int r = 0; r < R; r++) {
                double v = sn.rates.get(i, r);
                if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                    lo = Math.min(lo, v);
                    hi = Math.max(hi, v);
                }
            }
            if (Double.isInfinite(lo)) {
                dt.reason = String.format("station %d serves no class", i);
                return dt;
            }
            if (hi - lo > GlobalConstants.FineTol * hi) {
                dt.reason = String.format("station %d has a class-dependent service probability; the "
                        + "discrete-time cycle needs one Bernoulli server per node", i);
                return dt;
            }
            p[i] = lo;
            if (p[i] <= 0 || p[i] >= 1) {
                dt.reason = String.format("station %d has service probability %g; the product form "
                        + "of theorem 3.2 needs p in (0,1)", i, p[i]);
                return dt;
            }
        }

        int[] order = cycleOrder(sn);
        if (order == null) {
            dt.reason = "the stations do not form a single deterministic cycle; discrete-time FCFS "
                    + "networks of other topologies have no product form";
            return dt;
        }

        double[][] P = new double[M][N];
        for (int k = 0; k < M; k++) {
            double[] alpha = lldVector(lld(sn, order[k]), N);
            for (int n = 0; n < N; n++) {
                P[k][n] = p[order[k]] * alpha[n];
                if (P[k][n] <= 0 || P[k][n] > 1) {
                    dt.reason = "load dependence must keep every service probability inside (0,1]";
                    return dt;
                }
            }
        }

        dt.kind = "cycle";
        dt.order = order;
        dt.serviceProb = P;
        dt.population = N;
        return dt;
    }

    /**
     * Station order along the cycle, or null when the routing is not a single
     * deterministic cycle visiting every station exactly once.
     */
    private static int[] cycleOrder(NetworkStruct sn) {
        int M = sn.nstations;
        int R = sn.nclasses;
        Pair<Matrix, Matrix> rt = SnRtStations.snRtStations(sn);
        Matrix rtst = rt.getLeft();

        int[] succ = new int[M];
        for (int i = 0; i < M; i++) {
            int tgt = -1;
            for (int j = 0; j < M; j++) {
                double w = 0;
                for (int r = 0; r < R; r++) {
                    for (int s = 0; s < R; s++) {
                        w += rtst.get(i * R + r, j * R + s);
                    }
                }
                if (w > GlobalConstants.FineTol) {
                    if (Math.abs(w - R) > GlobalConstants.CoarseTol
                            && Math.abs(w - 1) > GlobalConstants.CoarseTol) {
                        return null;   // fractional routing out of station i
                    }
                    if (tgt != -1) {
                        return null;   // more than one successor
                    }
                    tgt = j;
                }
            }
            if (tgt == -1 || tgt == i) {
                return null;
            }
            succ[i] = tgt;
        }

        boolean[] visited = new boolean[M];
        int[] order = new int[M];
        int cur = 0;
        for (int k = 0; k < M; k++) {
            if (visited[cur]) {
                return null;
            }
            visited[cur] = true;
            order[k] = cur;
            cur = succ[cur];
        }
        if (cur != 0) {
            return null;
        }
        for (int i = 0; i < M; i++) {
            if (!visited[i]) {
                return null;
            }
        }
        return order;
    }

    private static SchedStrategy sched(NetworkStruct sn, int i) {
        if (sn.sched == null || sn.stations == null || i >= sn.stations.size()) {
            return null;
        }
        return sn.sched.get(sn.stations.get(i));
    }

    /** Effective buffer capacity of station ist for class r. */
    private static int capacity(NetworkStruct sn, int ist, int r) {
        double cap = Double.POSITIVE_INFINITY;
        if (sn.cap != null && sn.cap.getNumRows() > ist) {
            cap = Math.min(cap, sn.cap.get(ist, 0));
        }
        if (sn.classcap != null && sn.classcap.getNumRows() > ist
                && sn.classcap.getNumCols() > r) {
            cap = Math.min(cap, sn.classcap.get(ist, r));
        }
        if (Double.isInfinite(cap) || cap >= Integer.MAX_VALUE) {
            return Integer.MAX_VALUE;
        }
        return (int) cap;
    }

    /** Load-dependent scaling vector of station ist, null when undeclared. */
    private static double[] lld(NetworkStruct sn, int ist) {
        if (sn.lldscaling == null || sn.lldscaling.getNumRows() <= ist
                || sn.lldscaling.getNumCols() == 0) {
            return null;
        }
        int n = sn.lldscaling.getNumCols();
        double[] row = new double[n];
        boolean varies = false;
        for (int k = 0; k < n; k++) {
            row[k] = sn.lldscaling.get(ist, k);
            if (Math.abs(row[k] - 1.0) > GlobalConstants.FineTol) {
                varies = true;
            }
        }
        return varies ? row : null;
    }

    /**
     * Expand a load-dependent scaling to alpha(1..N), holding the last declared
     * value beyond the tabulated range, as the load-dependent NC solvers do.
     */
    private static double[] lldVector(double[] lld, int N) {
        double[] alpha = new double[N];
        java.util.Arrays.fill(alpha, 1.0);
        if (lld == null || lld.length == 0) {
            return alpha;
        }
        int n = Math.min(N, lld.length);
        System.arraycopy(lld, 0, alpha, 0, n);
        for (int k = lld.length; k < N; k++) {
            alpha[k] = lld[lld.length - 1];
        }
        return alpha;
    }

    /** Exact discrete-time analysis of sn. */
    public static NCResult solver_nc_dt_analyzer(NetworkStruct sn, SolverOptions options) {
        long tstart = System.nanoTime();
        DtModel dt = ncIsDtModel(sn);
        NCResult res = new NCResult();
        String method;
        if ("bernoulli1".equals(dt.kind)) {
            method = "dt.bernoulli1";
            dtSingle(sn, dt, res);
        } else if ("cycle".equals(dt.kind)) {
            boolean stateIndependent = true;
            for (int k = 0; k < dt.serviceProb.length && stateIndependent; k++) {
                for (int n = 1; n < dt.serviceProb[k].length; n++) {
                    if (Math.abs(dt.serviceProb[k][n] - dt.serviceProb[k][0])
                            > GlobalConstants.FineTol) {
                        stateIndependent = false;
                        break;
                    }
                }
            }
            method = stateIndependent ? "dt.cycle" : "dt.cycleld";
            dtCycle(sn, dt, method, res);
        } else {
            throw new RuntimeException("options.config.slotted was requested but the model is not "
                    + "a discrete-time product-form model: " + dt.reason);
        }

        double d = slotLength(options);
        if (d != 1.0) {
            res.TN.scale(1.0 / d);
            res.XN.scale(1.0 / d);
            res.RN.scale(d);
            res.CN.scale(d);
        }
        res.method = method;
        res.it = 1;
        res.runtime = (System.nanoTime() - tstart) / 1e9;
        return res;
    }

    private static void dtSingle(NetworkStruct sn, DtModel dt, NCResult res) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix Q = new Matrix(M, K);
        Matrix U = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        Matrix T = new Matrix(M, K);
        Matrix C = new Matrix(1, K);
        Matrix X = new Matrix(1, K);

        Bernoulli1Result b;
        if (dt.capacity == Integer.MAX_VALUE) {
            b = Dqsys_bernoulli1.dqsys_bernoulli1(dt.arrivalProb, dt.serviceProbSingle[0]);
        } else {
            b = Dqsys_bernoulli1.dqsys_bernoulli1(new double[]{dt.arrivalProb},
                    dt.serviceProbSingle, dt.capacity);
        }
        Q.set(dt.station, 0, b.meanQueueLength);
        U.set(dt.station, 0, b.utilization);
        T.set(dt.station, 0, b.throughput);
        R.set(dt.station, 0, b.meanSojournTime);
        X.set(0, 0, b.throughput);
        C.set(0, 0, b.meanSojournTime);
        // The Source row carries the offered stream, as on the continuous-time route.
        T.set(dt.source, 0, dt.arrivalProb);

        res.QN = Q;
        res.UN = U;
        res.RN = R;
        res.TN = T;
        res.CN = C;
        res.XN = X;
        res.lG = Math.log(b.normConst);
        res.STeff = new Matrix(M, K);
    }

    private static void dtCycle(NetworkStruct sn, DtModel dt, String method, NCResult res) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int N = dt.population;
        double[] Qs = new double[M];
        double[] Us = new double[M];
        double[] Ts = new double[M];
        double lG;

        if ("dt.cycle".equals(method)) {
            // State independent: propositions 3.18 and 3.19 end to end, with the
            // index of corollary 3.20(a) corrected (see Dpfqn_nc).
            double[] p = new double[M];
            for (int k = 0; k < M; k++) {
                p[k] = dt.serviceProb[k][0];
            }
            DpfqnNcResult nc = Dpfqn_nc.dpfqn_nc(p, N);
            lG = nc.lG;
            double xtput = nc.throughput();
            for (int k = 0; k < M; k++) {
                double q = 1.0 - p[k];
                double tail = 0.0;
                for (int n = 1; n <= N; n++) {
                    tail += Math.pow(q / p[k], n) / q * nc.G1[N - n + 1] / nc.G;
                }
                Qs[dt.order[k]] = tail;          // E[X_j] = sum_{n>=1} P(X_j>=n)
                Us[dt.order[k]] = xtput / p[k];  // P(X_j >= 1)
                Ts[dt.order[k]] = xtput;
            }
        } else {
            // State dependent: theorem 3.2 through the convolution of Dpfqn_ncld.
            DpfqnNcLdResult nc = Dpfqn_ncld.dpfqn_ncld(dt.serviceProb, N);
            lG = nc.lG;
            for (int k = 0; k < M; k++) {
                double[] marg = nc.marginal(k);
                double q = 0;
                double t = 0;
                for (int n = 0; n <= N; n++) {
                    q += marg[n] * n;
                    if (n >= 1) {
                        t += marg[n] * dt.serviceProb[k][n - 1];
                    }
                }
                Qs[dt.order[k]] = q;
                Us[dt.order[k]] = 1.0 - marg[0];
                Ts[dt.order[k]] = t;
            }
        }

        // Per-class split by population share; see the class comment.
        double[] share = new double[K];
        for (int r = 0; r < K; r++) {
            double v = sn.njobs.get(0, r);
            share[r] = (Double.isInfinite(v) || N == 0) ? 0.0 : v / N;
        }
        Matrix Q = new Matrix(M, K);
        Matrix U = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        Matrix T = new Matrix(M, K);
        Matrix C = new Matrix(1, K);
        Matrix X = new Matrix(1, K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                Q.set(i, r, Qs[i] * share[r]);
                U.set(i, r, Us[i] * share[r]);
                T.set(i, r, Ts[i] * share[r]);
                if (T.get(i, r) > 0) {
                    R.set(i, r, Q.get(i, r) / T.get(i, r));
                }
            }
        }
        for (int r = 0; r < K; r++) {
            int ref = (int) sn.refstat.get(r, 0);
            if (ref >= 0 && ref < M) {
                X.set(0, r, T.get(ref, r));
            }
            if (X.get(0, r) > 0) {
                C.set(0, r, sn.njobs.get(0, r) / X.get(0, r));
            }
        }

        res.QN = Q;
        res.UN = U;
        res.RN = R;
        res.TN = T;
        res.CN = C;
        res.XN = X;
        res.lG = lG;
        res.STeff = new Matrix(M, K);
    }
}
