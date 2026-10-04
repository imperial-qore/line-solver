/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.api.npfqn.Npfqn_bnd_bpt;
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
 * Achievable-region LOWER bound on the mean response times of a multiclass
 * open Markovian network, valid for EVERY non-idling scheduling policy at
 * every station.
 *
 * <p>Port of matlab/src/solvers/BA/solver_ba_bpt_analyzer.m. The polyhedron is
 * the first-order linear-programming relaxation of the achievable region
 * ({@link Npfqn_bnd_bpt}); this analyzer maps the LINE model onto it and reads
 * the bound back per station and class.</p>
 *
 * <p>CLASS SPACE. The reference's "class" is a buffer: one exponential service
 * rate, one Markovian routing law. LINE's (station, job class) pair is exactly
 * that, so a pair carrying traffic becomes one LP class, the Source is absorbed
 * into the external arrival vector, and class switching needs no special
 * treatment because sn.rt already carries it.</p>
 *
 * <p>BOUND CONVENTION. R(i,r) is obtained by minimizing x over the polyhedron
 * with the objective set to the unit vector of that pair, so each entry is a
 * valid lower bound on its own. Q follows by Little's law from the bounded R
 * and the EXACT throughput T (an open network's per-class rates are fixed by
 * the traffic equations, not by the policy), and so does C. U is exact for the
 * same reason.</p>
 *
 * <p>TIGHTNESS. The relaxation is exact on M/M/1 and tight on the externally
 * fed classes, but weak on a class whose arrivals are all internal: the only
 * term coupling x_r to the second-moment block carries the factor lambda0_r,
 * so an internally fed class can fall back to its own mean service time.</p>
 *
 * <p>Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994).
 * Optimization of multiclass queueing networks: polyhedral and nonlinear
 * characterizations of achievable performance. Annals of Applied Probability
 * 4(1), 43-75.</p>
 */
public final class Solver_ba_bpt_analyzer {

    private Solver_ba_bpt_analyzer() {}

    public static MVAResult solver_ba_bpt_analyzer(NetworkStruct sn, SolverOptions options) {
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
                        "Method 'bpt.lower' supports fully open networks only (no closed classes).");
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
                    "Method 'bpt.lower' requires an open network with a Source station.");
        }
        for (int a = 0; a < qstat.size(); a++) {
            int i = qstat.get(a);
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                throw new RuntimeException("Method 'bpt.lower' does not support delay "
                        + "(infinite-server) stations: the achievable region is derived for one "
                        + "server per station.");
            }
            if (sn.nservers.get(i, 0) > 1) {
                throw new RuntimeException(
                        "Method 'bpt.lower' does not support multi-server stations.");
            }
        }

        // ---- station-space routing, with the Source absorbed into lambda0 ----
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

        double[] lambda0 = new double[np];
        for (int si = 0; si < srcList.size(); si++) {
            int s = srcList.get(si);
            for (int r0 = 0; r0 < K; r0++) {
                double arr = sn.rates.get(s, r0);
                if (!Double.isFinite(arr) || arr <= 0) {
                    continue;
                }
                int srow = s * K + r0;
                for (int p = 0; p < np; p++) {
                    lambda0[p] += arr * rtst.get(srow, pairFlat[p]);
                }
            }
        }

        // Pair-to-pair routing. Flow to the Sink or back to a Source is the
        // exit probability, i.e. the row deficit, and needs no column.
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
        double[] lam0k = new double[nk];
        double[] muk = new double[nk];
        int[] statk = new int[nk];
        int[] clsk = new int[nk];
        Matrix Pk = new Matrix(nk, nk);
        List<Integer> ustat = new ArrayList<Integer>();
        for (int a = 0; a < nk; a++) {
            int p = keep.get(a);
            lam0k[a] = lambda0[p];
            statk[a] = pairStation[p];
            clsk[a] = pairClass[p];
            for (int b = 0; b < nk; b++) {
                double v = Pfull.get(p, keep.get(b));
                if (v != 0) {
                    Pk.set(a, b, v);
                }
            }
            muk[a] = sn.rates.get(statk[a], clsk[a]);
            if (!Double.isFinite(muk[a]) || muk[a] <= 0) {
                throw new RuntimeException("Station " + (statk[a] + 1) + " has no service rate for "
                        + "class " + (clsk[a] + 1) + " but carries its traffic.");
            }
            ProcessType pt = sn.procid.get(sn.stations.get(statk[a])).get(sn.jobclasses.get(clsk[a]));
            if (pt != ProcessType.EXP) {
                throw new RuntimeException("Method 'bpt.lower' requires exponential service: "
                        + "station " + (statk[a] + 1) + " class " + (clsk[a] + 1) + " is " + pt + ".");
            }
            if (!ustat.contains(statk[a])) {
                ustat.add(statk[a]);
            }
        }
        // Dense station index space for the LP.
        int[] stationOf = new int[nk];
        for (int a = 0; a < nk; a++) {
            stationOf[a] = ustat.indexOf(statk[a]);
        }

        // ---- one LP per pair, objective = that pair's unit vector ----
        Npfqn_bnd_bpt.Result info = null;
        for (int a = 0; a < nk; a++) {
            double[] e = new double[nk];
            e[a] = 1.0;
            info = Npfqn_bnd_bpt.npfqn_bnd_bpt(lam0k, muk, Pk, stationOf, e);
            ret.RN.set(statk[a], clsk[a], info.zlb);
            ret.TN.set(statk[a], clsk[a], info.lambda[a]);
            ret.UN.set(statk[a], clsk[a], info.rho[a]);
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

        ret.runtime = (System.nanoTime() - t0) / 1.0e9;
        return ret;
    }
}
