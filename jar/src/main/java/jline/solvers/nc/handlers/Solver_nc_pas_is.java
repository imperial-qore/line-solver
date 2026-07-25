/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.handlers;

import jline.api.pfqn.nc.Pas_swap2order;
import jline.api.pfqn.nc.Pfqn_pas_is;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.nodes.Node;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

/**
 * Importance-sampling (IS) normalizing-constant analysis of a closed two-station
 * pass-and-swap (P&S) tandem with a non-empty swap graph (Casale, Comte and
 * Dorsman, 2026). Monte-Carlo counterpart of {@link Solver_nc_oi} for the case
 * that the exact OI convolution does not apply (non-empty swap graph): the
 * ordered-state chain is reducible and the recurrent communicating class carries
 * a per-class product form whose constant G_C is estimated by {@link Pfqn_pas_is}.
 *
 * <p>Station 1 is the upstream P&S queue (prefix of the ordering), station 2 the
 * downstream queue (reversed suffix). The global placement order is derived from
 * the swap graph via {@link Pas_swap2order}. Mean per-class queue lengths come
 * directly from the auto-normalized IS; per-class throughput uses the ratio
 * X_r = G(N - e_r)/G(N) with common random numbers; utilization and response
 * time follow the conventions of {@link Solver_nc_oi}.
 *
 * <p>Port of matlab/src/solvers/NC/solver_nc_pas_is_analyzer.m and nc_is_pas_model.m.
 */
public final class Solver_nc_pas_is {
    private Solver_nc_pas_is() {}

    /**
     * True when the model is a closed two-station P&S tandem: both stations
     * OI/PAS, no other station. The swap graph may be empty, since an
     * order-independent queue is exactly the P&S specialization with an empty
     * swap graph, and Pfqn_pas_is with H=0 samples all microstates (the OI
     * case). Used to bind the importance-sampling selectors ('is', and
     * 'sampling' which maps to 'is' when OI/PAS stations are present). A P&S
     * tandem with a NON-EMPTY swap graph is reducible and lies outside the exact
     * path of {@link Solver_nc_oi#nc_is_oi_model}, so it also takes this
     * analyzer on 'default'; a pure-OI tandem on 'default'/'exact' is caught
     * earlier by the exact OI analyzer.
     */
    public static boolean nc_is_pas_model(NetworkStruct sn) {
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                return false;
            }
        }
        if (sn.nstations != 2) {
            return false;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s != SchedStrategy.PAS && s != SchedStrategy.OI) {
                return false;
            }
            Node node = sn.nodes.get((int) sn.stationToNode.get(ist));
            NodeParam np = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
            if (!(np instanceof QueueNodeParam)) {
                return false;
            }
            if (((QueueNodeParam) np).svcRateFun == null) {
                return false;   // no OI rank-rate function: cannot evaluate the balance function
            }
        }
        return true;
    }

    public static SolverNC.SolverNCReturn solver_nc_pas_is(NetworkStruct sn, SolverOptions options) {
        long tStart = System.nanoTime();
        int M = sn.nstations;
        int K = sn.nclasses;

        // ---- reject class switching (P&S rank rates are per raw class) ------
        for (int c = 0; c < sn.nchains; c++) {
            if (sn.inchain.get(c).getNumElements() > 1) {
                throw new RuntimeException("solver_nc_pas_is requires one class per chain (no class switching).");
            }
        }
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                throw new RuntimeException("solver_nc_pas_is requires a closed queueing network.");
            }
        }
        if (M != 2) {
            throw new RuntimeException("solver_nc_pas_is models a two-station pass-and-swap tandem (got " + M + " stations).");
        }
        int[] N = new int[K];
        for (int r = 0; r < K; r++) {
            N[r] = (int) Math.round(sn.njobs.get(r));
        }

        // ---- classify stations and read swap graph + service rate functions -
        SerializableFunction<Matrix, Double>[] svc = getSvc(sn, M);
        Matrix[] swapG = new Matrix[M];
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s != SchedStrategy.PAS && s != SchedStrategy.OI) {
                throw new RuntimeException("solver_nc_pas_is requires both stations to be OI/PAS (station " + ist + " is not).");
            }
            Node node = sn.nodes.get((int) sn.stationToNode.get(ist));
            NodeParam np = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
            if (!(np instanceof QueueNodeParam)) {
                throw new RuntimeException("station " + ist + " has no OI/PAS node parameters.");
            }
            if (svc[ist] == null) {
                throw new RuntimeException("OI/PAS station " + ist + " has no service rate function; set it via setService(c -> ...).");
            }
            swapG[ist] = ((QueueNodeParam) np).swapGraph;
        }

        // The stored swap graph is the raw (undirected) class-compatibility
        // graph; derive the GLOBAL placement-order DAG H (defining the recurrent
        // communicating class D) from the P&S dynamics on the single-job-per-
        // class instance. Station 2 is the reversed suffix inside Pfqn_pas_is.
        int[] N0 = new int[K];
        for (int r = 0; r < K; r++) {
            N0[r] = 1;
        }
        int[][] H = Pas_swap2order.pas_swap2order(swapG[0], swapG[1], svc[0], svc[1], N0);

        // ---- per-class visits (chain == class); require unit visits ---------
        Matrix V = new Matrix(M, K);
        V.zero();
        for (int r = 0; r < K; r++) {
            int c = -1;
            for (int cc = 0; cc < sn.nchains; cc++) {
                if (sn.chains.get(cc, r) != 0) {
                    c = cc;
                    break;
                }
            }
            Matrix vis = sn.visits.get(c);
            for (int ist = 0; ist < M; ist++) {
                int isf = (int) sn.stationToStateful.get(ist);
                V.set(ist, r, vis.get(isf, r));
            }
            double vref = V.get((int) sn.refstat.get(r), r);
            if (vref > 0) {
                for (int ist = 0; ist < M; ist++) {
                    V.set(ist, r, V.get(ist, r) / vref);
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                if (N[r] > 0 && Math.abs(V.get(ist, r) - 1) > 1e-9) {
                    throw new RuntimeException("solver_nc_pas_is requires unit per-class visits (station " + ist + ", class " + r + ", V=" + V.get(ist, r) + ").");
                }
            }
        }

        // ---- OI rank-rate handles on a per-class count vector ---------------
        List<ToDoubleFunction<int[]>> mu = new ArrayList<ToDoubleFunction<int[]>>();
        mu.add(makeRankRate(svc[0]));
        mu.add(makeRankRate(svc[1]));

        long nsamples = options.samples > 0 ? options.samples : 10000L;
        long seed = options.seed;
        boolean verbose = false;

        // ---- normalizing constant and mean queue lengths at population N ----
        Pfqn_pas_is.Result res = Pfqn_pas_is.pfqn_pas_is(N, mu, H, nsamples, seed, verbose);
        double G = res.G;
        double lG = res.lG;

        Matrix Q = new Matrix(M, K);
        Q.zero();
        for (int r = 0; r < K; r++) {
            Q.set(0, r, res.Q[0][r]);
            Q.set(1, r, res.Q[1][r]);
        }

        // ---- per-class throughput X_r = G(N - e_r)/G(N) (common randoms) ----
        Matrix X = new Matrix(1, K);
        X.zero();
        for (int r = 0; r < K; r++) {
            if (N[r] > 0) {
                int[] Nm = N.clone();
                Nm[r]--;
                double Gr = Pfqn_pas_is.pfqn_pas_is(Nm, mu, H, nsamples, seed, false).G;
                if (G > 0) {
                    X.set(r, Gr / G);
                }
            }
        }

        // ---- throughput, utilization, response time -------------------------
        Matrix T = new Matrix(M, K);
        Matrix U = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        T.zero();
        U.zero();
        R.zero();
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                T.set(ist, r, X.get(r) * V.get(ist, r));
            }
        }
        for (int ist = 0; ist < M; ist++) {
            double S = sn.nservers.get(ist);
            if (!Double.isFinite(S) || S <= 0) {
                S = 1;
            }
            for (int r = 0; r < K; r++) {
                if (N[r] > 0) {
                    int[] er = new int[K];
                    er[r] = 1;
                    double muR = mu.get(ist).applyAsDouble(er);
                    if (muR > 0) {
                        U.set(ist, r, T.get(ist, r) / muR / S);
                    }
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                if (T.get(ist, r) > 0) {
                    R.set(ist, r, Q.get(ist, r) / T.get(ist, r));
                }
            }
        }

        Matrix STeff = new Matrix(M, K);
        STeff.zero();
        double runtime = (System.nanoTime() - tStart) / 1.0e9;
        return new SolverNC.SolverNCReturn(Q, U, R, T, sn.nchains, X, lG, STeff, 1, runtime, "is");
    }

    // ======================================================================
    @SuppressWarnings("unchecked")
    private static SerializableFunction<Matrix, Double>[] getSvc(NetworkStruct sn, int M) {
        SerializableFunction<Matrix, Double>[] svc = new SerializableFunction[M];
        for (int ist = 0; ist < M; ist++) {
            Node node = sn.nodes.get((int) sn.stationToNode.get(ist));
            NodeParam np = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
            svc[ist] = (np instanceof QueueNodeParam) ? ((QueueNodeParam) np).svcRateFun : null;
        }
        return svc;
    }

    private static ToDoubleFunction<int[]> makeRankRate(final SerializableFunction<Matrix, Double> fun) {
        // OI rank rate on a per-class count/support vector n. svcRateFun(c) takes
        // an ordered microstate list; permutation-invariant for an OI station, so
        // evaluate on a canonical 0-based microstate with n_r copies of class r.
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                int len = 0;
                for (int v : n) {
                    if (v > 0) {
                        len += v;
                    }
                }
                Matrix c = new Matrix(1, len);
                int col = 0;
                for (int r = 0; r < n.length; r++) {
                    for (int t = 0; t < n[r]; t++) {
                        c.set(0, col++, r);
                    }
                }
                return fun.apply(c);
            }
        };
    }

}
