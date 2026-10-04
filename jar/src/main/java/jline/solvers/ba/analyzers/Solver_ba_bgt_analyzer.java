/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.api.npfqn.Npfqn_bnd_bgt;
import jline.api.sn.SnRtStations;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
 * multitype open Markovian network, valid for EVERY work-conserving Markovian
 * policy.
 *
 * <p>Port of matlab/src/solvers/BA/solver_ba_bgt_analyzer.m. The polyhedron and
 * the bound are {@link Npfqn_bnd_bgt}; this analyzer maps the LINE model onto
 * them and reads the bound back per station and class.</p>
 *
 * <p>CLASS SPACE. The reference's network is a MULTITYPE one: each type follows
 * a FIXED sequence of stages, and stage k of type i is its own buffer. LINE's
 * (station, job class) pair is that buffer, so the analyzer walks the routing
 * matrix from the Source and turns each open class into one type whose stages
 * are the pairs it visits. Two gates follow from the model and are enforced by
 * name rather than approximated: routing must be DETERMINISTIC (a pair sends
 * everything to one successor, or everything to the Sink), and routes must NOT
 * MERGE (a pair belongs to exactly one type, else the reference's class index
 * (i,k) is not defined). A re-entrant line is expressible by giving the
 * revisits distinct LINE classes.</p>
 *
 * <p>THE BOUND IS LOOSE, and knowingly so: the exception parameter of the
 * smoothed Lyapunov function carries (Lmax+gamma)^3/gamma^2 and dominates as
 * soon as there is more than one station. What is sharp is the STABILITY
 * CERTIFICATE and the geometric tail RATE.</p>
 *
 * <p>Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance
 * of multiclass Markovian queueing networks via piecewise linear Lyapunov
 * functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1.</p>
 */
public final class Solver_ba_bgt_analyzer {

    private Solver_ba_bgt_analyzer() {}

    public static MVAResult solver_ba_bgt_analyzer(NetworkStruct sn, SolverOptions options) {
        long t0 = System.nanoTime();
        final int M = sn.nstations;
        final int K = sn.nclasses;

        MVAResult ret = new MVAResult();
        ret.QN = new Matrix(M, K);
        ret.UN = new Matrix(M, K);
        ret.RN = new Matrix(M, K);
        ret.TN = new Matrix(M, K);
        ret.CN = new Matrix(1, K);
        ret.XN = new Matrix(1, K);
        ret.logNormConstAggr = Double.NaN;
        ret.iter = 1;

        // ---- model gates ----
        for (int r = 0; r < sn.njobs.length(); r++) {
            if (Double.isFinite(sn.njobs.get(r))) {
                throw new RuntimeException(
                        "Method 'bgt.upper' supports fully open networks only (no closed classes).");
            }
        }
        List<Integer> srcList = new ArrayList<Integer>();
        List<Integer> qstat = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
                srcList.add(i);
            } else {
                qstat.add(i);
            }
        }
        if (srcList.isEmpty()) {
            throw new RuntimeException(
                    "Method 'bgt.upper' requires an open network with a Source station.");
        }
        for (int a = 0; a < qstat.size(); a++) {
            int i = qstat.get(a);
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                throw new RuntimeException("Method 'bgt.upper' does not support delay "
                        + "(infinite-server) stations: the reference's network has one server "
                        + "per station.");
            }
            if (sn.nservers.get(i, 0) > 1) {
                throw new RuntimeException(
                        "Method 'bgt.upper' does not support multi-server stations.");
            }
        }

        Pair<Matrix, Matrix> rtstPair = SnRtStations.snRtStations(sn);
        Matrix rtst = rtstPair.getLeft();

        int np = qstat.size() * K;
        int[] pairStation = new int[np];
        int[] pairClass = new int[np];
        int[] pairFlat = new int[np];
        int p = 0;
        for (int a = 0; a < qstat.size(); a++) {
            int i = qstat.get(a);
            for (int r = 0; r < K; r++) {
                pairStation[p] = i;
                pairClass[p] = r;
                pairFlat[p] = i * K + r;
                p++;
            }
        }

        // ---- walk one deterministic route per source class ----
        List<Double> lambdaL = new ArrayList<Double>();
        List<int[]> routes = new ArrayList<int[]>();
        boolean[] used = new boolean[np];
        for (int si = 0; si < srcList.size(); si++) {
            int s = srcList.get(si);
            for (int r0 = 0; r0 < K; r0++) {
                double arr = sn.rates.get(s, r0);
                if (!Double.isFinite(arr) || arr <= 0) {
                    continue;
                }
                int cur = singleSuccessor(rtst, s * K + r0, pairFlat,
                        "the Source for class " + (r0 + 1));
                List<Integer> route = new ArrayList<Integer>();
                while (cur >= 0) {
                    if (used[cur]) {
                        throw new RuntimeException("Method 'bgt.upper' needs routes that do not "
                                + "merge: station " + (pairStation[cur] + 1) + " class "
                                + (pairClass[cur] + 1) + " is visited by more than one type. "
                                + "Give the visits distinct job classes.");
                    }
                    used[cur] = true;
                    route.add(cur);
                    cur = singleSuccessor(rtst, pairFlat[cur], pairFlat, "station "
                            + (pairStation[route.get(route.size() - 1)] + 1) + " class "
                            + (pairClass[route.get(route.size() - 1)] + 1));
                }
                if (route.isEmpty()) {
                    throw new RuntimeException("Class " + (r0 + 1)
                            + " leaves the Source and reaches no station.");
                }
                int[] rv = new int[route.size()];
                for (int k = 0; k < rv.length; k++) {
                    rv[k] = route.get(k);
                }
                routes.add(rv);
                lambdaL.add(arr);
            }
        }
        if (routes.isEmpty()) {
            throw new RuntimeException("The model carries no open traffic.");
        }

        final int I = routes.size();
        double[] lambda = new double[I];
        double[][] mu = new double[I][];
        int[][] sigma = new int[I][];
        List<Integer> ustat = new ArrayList<Integer>();
        for (int i = 0; i < I; i++) {
            lambda[i] = lambdaL.get(i);
            int[] route = routes.get(i);
            mu[i] = new double[route.length];
            sigma[i] = new int[route.length];
            for (int k = 0; k < route.length; k++) {
                int q = route[k];
                mu[i][k] = sn.rates.get(pairStation[q], pairClass[q]);
                if (!Double.isFinite(mu[i][k]) || mu[i][k] <= 0) {
                    throw new RuntimeException("Station " + (pairStation[q] + 1)
                            + " has no service rate for class " + (pairClass[q] + 1)
                            + " but carries its traffic.");
                }
                ProcessType pt = sn.procid.get(sn.stations.get(pairStation[q]))
                        .get(sn.jobclasses.get(pairClass[q]));
                if (pt != ProcessType.EXP) {
                    throw new RuntimeException("Method 'bgt.upper' requires exponential service: "
                            + "station " + (pairStation[q] + 1) + " class " + (pairClass[q] + 1)
                            + " is " + pt + ".");
                }
                if (!ustat.contains(pairStation[q])) {
                    ustat.add(pairStation[q]);
                }
            }
        }
        // Dense station index space for the LP.
        for (int i = 0; i < I; i++) {
            int[] route = routes.get(i);
            for (int k = 0; k < route.length; k++) {
                sigma[i][k] = ustat.indexOf(pairStation[route[k]]);
            }
        }

        Npfqn_bnd_bgt.Result info = Npfqn_bnd_bgt.npfqn_bnd_bgt(lambda, mu, sigma, ustat.size());

        // ---- read the bound back per station and class ----
        for (int i = 0; i < I; i++) {
            int[] route = routes.get(i);
            for (int k = 0; k < route.length; k++) {
                int q = route[k];
                int ist = pairStation[q], r = pairClass[q];
                ret.QN.set(ist, r, ret.QN.get(ist, r) + info.Qub[i][k]);
                ret.TN.set(ist, r, ret.TN.get(ist, r) + lambda[i]);
                ret.UN.set(ist, r, ret.UN.get(ist, r) + lambda[i] / sn.rates.get(ist, r));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                if (ret.TN.get(i, r) > 0) {
                    ret.RN.set(i, r, ret.QN.get(i, r) / ret.TN.get(i, r));
                }
            }
        }

        // ---- exact open-network quantities ----
        for (int si = 0; si < srcList.size(); si++) {
            int s = srcList.get(si);
            for (int r = 0; r < K; r++) {
                double arr = sn.rates.get(s, r);
                if (Double.isFinite(arr) && arr > 0) {
                    ret.TN.set(s, r, ret.TN.get(s, r) + arr);
                    ret.XN.set(0, r, ret.XN.get(0, r) + arr);
                }
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

        ret.runtime = (System.nanoTime() - t0) / 1.0e9;
        return ret;
    }

    /**
     * The single successor of a routing row, as a pair index, or -1 when
     * everything leaves the network. A probabilistic split is refused by name:
     * the reference's network has deterministic routing and a split is a
     * different model, not an approximation of this one.
     */
    private static int singleSuccessor(Matrix rtst, int row, int[] pairFlat, String who) {
        final double tol = 1e-9;
        double mass = 0, best = 0;
        int p = -1;
        for (int q = 0; q < pairFlat.length; q++) {
            double v = rtst.get(row, pairFlat[q]);
            if (v > tol) {
                mass += v;
                if (v > best) {
                    best = v;
                    p = q;
                }
            }
        }
        if (mass <= tol) {
            return -1;
        }
        if (Math.abs(mass - 1) > tol || Math.abs(best - 1) > tol) {
            throw new RuntimeException("Method 'bgt.upper' needs deterministic routing: " + who
                    + " splits its departures (the largest branch carries " + best + " of them). "
                    + "The reference's network routes each type along a fixed sequence of stages.");
        }
        return p;
    }
}
