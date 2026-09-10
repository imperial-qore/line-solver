/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.Pfqn_sdr;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

/**
 * Exact product-form analysis of a closed multiclass network with the
 * state-dependent routing of Krzesinski (1987), "Multiclass Queueing Networks
 * with State-Dependent Routing", Performance Evaluation 7(2):125-143.
 *
 * <p>The joint distribution is eq. (16); the coefficients xi are those of
 * Section 3.2, obtained from the state-independent part of the routing
 * matrix.</p>
 */
public class Solver_nc_sdr_analyzer {

    /**
     * Analyzes a network whose entry center routes by state-dependent routing.
     *
     * @param sn the network structure, carrying a non-null sdr
     * @param options solver options
     * @return the mean performance measures and the log normalizing constant
     */
    public static NCResult solver_nc_sdr_analyzer(NetworkStruct sn, SolverOptions options) {
        long tstart = System.nanoTime();
        int M = sn.nstations;
        int R = sn.nclasses;

        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                throw new RuntimeException("State-dependent routing is defined for closed networks only; "
                        + "the model has an open class.");
            }
        }
        if (sn.nchains != R) {
            throw new RuntimeException("State-dependent routing does not support class switching: the product form "
                    + "of Krzesinski (1987) is stated over closed chains whose customers keep their class. Merge the "
                    + "switching classes into one class.");
        }
        if (sn.nstateful != M) {
            throw new RuntimeException("State-dependent routing requires every stateful node to be a station: the "
                    + "product form is over queue lengths, and a stateless node holds none.");
        }

        Matrix S = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(i, r);
                if (rate > 0 && !Double.isInfinite(rate) && !Double.isNaN(rate)) {
                    S.set(i, r, 1.0 / rate);
                }
            }
        }

        int Ntot = 0;
        for (int r = 0; r < R; r++) {
            Ntot += (int) sn.njobs.get(r);
        }
        Matrix alpha = new Matrix(M, Math.max(1, Ntot));
        for (int i = 0; i < M; i++) {
            for (int k = 1; k <= alpha.getNumCols(); k++) {
                double a = 1.0;
                if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                    a = k;
                } else {
                    double c = sn.nservers.get(i);
                    if (!Double.isInfinite(c) && c > 1) {
                        a = Math.min(k, c);
                    }
                }
                if (sn.lldscaling != null && !sn.lldscaling.isEmpty()
                        && i < sn.lldscaling.getNumRows() && (k - 1) < sn.lldscaling.getNumCols()) {
                    double lld = sn.lldscaling.get(i, k - 1);
                    if (lld > 0) {
                        a *= lld;
                    }
                }
                alpha.set(i, k - 1, a);
            }
        }

        // A BCMP center served FCFS must hold one rate for every chain; eq. (16)
        // admits chain-dependent rates only at the symmetric disciplines.
        for (int i = 0; i < M; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                double ref = -1;
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) > 0 && S.get(i, r) > 0) {
                        if (ref < 0) {
                            ref = S.get(i, r);
                        } else if (Math.abs(ref - S.get(i, r)) > 1e-12) {
                            throw new RuntimeException("Station " + sn.stations.get(i).getName()
                                    + " is FCFS with chain-dependent service times, which has no BCMP product form. "
                                    + "Use PS, LCFSPR or INF, or equalize the service times.");
                        }
                    }
                }
            }
        }

        List<Matrix> P = new ArrayList<Matrix>();
        for (int r = 0; r < R; r++) {
            Matrix Pr = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                int isf = (int) sn.stationToStateful.get(i);
                for (int j = 0; j < M; j++) {
                    int jsf = (int) sn.stationToStateful.get(j);
                    Pr.set(i, j, sn.rt.get(isf * R + r, jsf * R + r));
                }
            }
            P.add(Pr);
        }
        Matrix xi = Pfqn_sdr.pfqn_sdrvisits(sn.sdr, P);

        Matrix N = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            N.set(0, r, sn.njobs.get(r));
        }
        // 'sdr' evaluates the product form (16) exactly by state enumeration,
        // which is general in the branch topology; 'sdr.mva' runs the paper's
        // Section 4 MVA and convolution, which costs O(J T M (V_1...V_J)^2)
        // instead of the state-space size but requires single-centre branches.
        // Both are exact.
        String method = "sdr";
        if (options != null && "sdr.mva".equalsIgnoreCase(options.method)) {
            method = "sdr.mva";
        }
        Pfqn_sdr.Result pf = "sdr.mva".equals(method)
                ? Pfqn_sdr.pfqn_sdrmva(S, xi, N, sn.sdr, alpha)
                : Pfqn_sdr.pfqn_sdr(S, xi, N, sn.sdr, alpha);

        NCResult res = new NCResult();
        res.QN = pf.Q;
        res.UN = pf.U;
        res.RN = pf.R;
        res.TN = pf.X;
        res.lG = pf.lG;
        res.CN = new Matrix(1, R);
        res.XN = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            int ref = (int) sn.refstat.get(r);
            double x = pf.X.get(ref, r);
            res.XN.set(0, r, x);
            if (x > 0) {
                res.CN.set(0, r, sn.njobs.get(r) / x);
            }
        }
        res.it = 1;
        res.method = method;
        res.runtime = (System.nanoTime() - tstart) / 1.0e9;
        return res;
    }
}
