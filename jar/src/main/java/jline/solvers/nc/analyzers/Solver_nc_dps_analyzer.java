/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc.analyzers;

import jline.api.npfqn.Npfqn_dps_morrison;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

/**
 * Heavy-usage asymptotic analysis of the closed two-station network with one think
 * (infinite-server) station and one discriminatory processor-sharing station, by the
 * generating-function expansion of J.A. Morrison, "Asymptotic analysis of a large closed queueing
 * network with discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.
 *
 * <p>Admitted only on the exact shape {@link #nc_is_dps_model} tests for. The kernel is
 * {@link Npfqn_dps_morrison}; this class maps the model struct onto it and lifts the per-class DPS
 * results into the station-by-class matrices the NC analyzers return.</p>
 *
 * <p>THERE IS NO NORMALIZING CONSTANT HERE. A DPS station is not product-form -- that is the
 * premise of the paper -- so lG is NaN, as on the maximum-entropy route. NC hosts this method
 * because NC is where LINE keeps the asymptotic expansions of generating functions and
 * normalizing-constant integrals (pana, mmint2, le, ble, gleint, rayint), which is the family
 * Morrison's expansion belongs to, not because a constant is being computed.</p>
 *
 * <p>Response times come from Little's law on the queue-length result rather than from the expanded
 * RESULT 2 (eq. 4.17), so that Q = R*T holds exactly in the returned table; the two agree to the
 * order of the approximation, since Morrison derives (4.17) as the ratio (4.11)/(4.15).</p>
 */
public class Solver_nc_dps_analyzer {

    /**
     * True when the model is the closed two-station network Morrison's expansion is derived for:
     * one infinite-server (think) station and one single-server DPS station, exponential service,
     * every class alternating between the two. The shape is checked exactly, not approximately:
     * outside it the expansion has no derivation behind it.
     *
     * @param sn the network structure
     * @return true when the Morrison route applies
     */
    public static boolean nc_is_dps_model(NetworkStruct sn) {
        if (sn == null || sn.njobs == null) {
            return false;
        }
        boolean anyPositive = false;
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                return false;
            }
            if (sn.njobs.get(r) > 0) {
                anyPositive = true;
            }
        }
        if (!anyPositive || sn.nstations != 2) {
            return false;
        }
        int iInf = -1;
        int iDps = -1;
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.INF) {
                iInf = (iInf == -1) ? ist : -2;
            } else if (s == SchedStrategy.DPS) {
                iDps = (iDps == -1) ? ist : -2;
            } else {
                return false;
            }
        }
        if (iInf < 0 || iDps < 0) {
            return false;
        }
        double c = sn.nservers.get(iDps);
        if (Double.isFinite(c) && c != 1) {
            return false;   // multi-server DPS: the min(n,c) share is not Morrison's
        }
        int K = sn.nclasses;
        if (sn.nchains != K) {
            return false;   // class switching
        }
        boolean lld = sn.lldscaling != null && sn.lldscaling.getNumElements() > 0;
        boolean cd = sn.cdscaling != null && !sn.cdscaling.isEmpty();
        boolean jd = sn.jdscaling != null && !sn.jdscaling.isEmpty();
        if (lld || cd || jd) {
            return false;
        }
        for (int r = 0; r < K; r++) {
            if (sn.njobs.get(r) <= 0) {
                return false;
            }
            for (int ist : new int[]{iInf, iDps}) {
                double rate = sn.rates.get(ist, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    return false;
                }
                double scv = sn.scv.get(ist, r);
                if (Double.isFinite(scv) && Math.abs(scv - 1) > 1e-6) {
                    return false;
                }
            }
            double wgt = sn.schedparam.get(iDps, r);
            if (!Double.isFinite(wgt) || wgt <= 0) {
                return false;
            }
        }
        Matrix V = visits(sn);
        if (V == null) {
            return false;
        }
        for (int r = 0; r < K; r++) {
            if (Math.abs(V.get(iInf, r) - V.get(iDps, r)) > 1e-9 * Math.max(1.0, V.get(iInf, r))) {
                return false;   // unequal visits: not the alternating cycle of the paper
            }
        }
        return true;
    }

    /**
     * Analyzes the closed think+DPS network.
     *
     * @param sn the network structure, of the shape {@link #nc_is_dps_model} accepts
     * @param options solver options
     * @return the mean performance measures, with lG = NaN
     */
    public static NCResult solver_nc_dps_analyzer(NetworkStruct sn, SolverOptions options) {
        long tstart = System.nanoTime();
        // The shape is re-checked here, not assumed from the caller: this analyzer is
        // reachable from SolverNC, the dispatch chain and directly from user code, and
        // every quantity below -- think time, DPS service time, weights, visit ratios --
        // is meaningless off the shape the expansion was derived for.
        if (!nc_is_dps_model(sn)) {
            throw new RuntimeException("solver_nc_dps_analyzer applies only to a CLOSED network of exactly two "
                    + "stations, one infinite-server (think) station and one single-server DPS station with "
                    + "exponential service and one visit each per cycle (see nc_is_dps_model).");
        }
        int M = sn.nstations;
        int K = sn.nclasses;
        int iInf = -1;
        int iDps = -1;
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.INF) {
                iInf = ist;
            } else if (s == SchedStrategy.DPS) {
                iDps = ist;
            }
        }
        if (iInf < 0 || iDps < 0) {
            throw new RuntimeException("Solver_nc_dps_analyzer requires one INF and one DPS station.");
        }

        Matrix Npop = new Matrix(1, K);
        Matrix Z = new Matrix(1, K);
        Matrix S = new Matrix(1, K);
        Matrix w = new Matrix(1, K);
        for (int r = 0; r < K; r++) {
            Npop.set(0, r, sn.njobs.get(r));
            Z.set(0, r, 1.0 / sn.rates.get(iInf, r));
            S.set(0, r, 1.0 / sn.rates.get(iDps, r));
            w.set(0, r, sn.schedparam.get(iDps, r));
        }

        Npfqn_dps_morrison.Result mor = Npfqn_dps_morrison.npfqn_dps_morrison(Npop, Z, S, w);

        Matrix V = visits(sn);
        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix XN = new Matrix(1, K);
        Matrix CN = new Matrix(1, K);

        double c = sn.nservers.get(iDps);
        if (!Double.isFinite(c) || c <= 0) {
            c = 1;
        }
        for (int r = 0; r < K; r++) {
            double q = mor.Q.get(0, r);
            double np = Npop.get(0, r);
            if (!(q >= 0) || q > np) {
                // outside the moderately-heavy regime the expansion can leave [0,N]
                q = Math.min(Math.max(q, 0.0), np);
            }
            double x = (np - q) / Z.get(0, r);
            QN.set(iDps, r, q);
            QN.set(iInf, r, np - q);          // population conservation (exact, closed)
            XN.set(0, r, x);
            for (int ist = 0; ist < M; ist++) {
                TN.set(ist, r, x * V.get(ist, r));
            }
            UN.set(iInf, r, QN.get(iInf, r));  // INF utilization convention
            UN.set(iDps, r, x * V.get(iDps, r) * S.get(0, r) / c);
            if (x > 0) {
                CN.set(0, r, np / x);
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                if (TN.get(ist, r) > 0) {
                    RN.set(ist, r, QN.get(ist, r) / TN.get(ist, r));
                }
            }
        }

        NCResult res = new NCResult();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.XN = XN;
        res.CN = CN;
        res.lG = Double.NaN;                   // no product form: no normalizing constant
        res.it = 1;
        res.method = "morrison";
        res.runtime = (System.nanoTime() - tstart) / 1.0e9;
        return res;
    }

    /** Per-class visit matrix, normalized at the reference station; null when unavailable. */
    private static Matrix visits(NetworkStruct sn) {
        if (sn.visits == null || sn.chains == null || sn.stationToStateful == null || sn.refstat == null) {
            return null;
        }
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix V = new Matrix(M, K);
        for (int r = 0; r < K; r++) {
            int chain = -1;
            for (int c = 0; c < sn.chains.getNumRows(); c++) {
                if (sn.chains.get(c, r) != 0) {
                    if (chain != -1) {
                        return null;
                    }
                    chain = c;
                }
            }
            if (chain < 0) {
                return null;
            }
            Matrix vc = sn.visits.get(chain);
            if (vc == null) {
                return null;
            }
            for (int ist = 0; ist < M; ist++) {
                V.set(ist, r, vc.get((int) sn.stationToStateful.get(ist), r));
            }
            double vref = V.get((int) sn.refstat.get(r), r);
            if (!(vref > 0)) {
                return null;
            }
            for (int ist = 0; ist < M; ist++) {
                V.set(ist, r, V.get(ist, r) / vref);
            }
        }
        return V;
    }
}
