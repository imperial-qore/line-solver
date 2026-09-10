/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba.analyzers;

import jline.api.spn.Spn_lpbnd;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

/**
 * Linear-programming bounds on the mean marking and the throughputs of a
 * stochastic timed Petri net.
 *
 * <p>Port of matlab/src/solvers/BA/solver_ba_spnlp_analyzer.m. The polytope and
 * the LP are {@link Spn_lpbnd}; this analyzer maps the LINE model onto them and
 * reads one side of the bracket back per place and class.</p>
 *
 * <p>METHOD NAMES. Four, in two families: {@code spnlp.upper} and
 * {@code spnlp.lower} are the Markovian LP and need exponential firing times;
 * {@code spnlp.op.upper} and {@code spnlp.op.lower} drop the second-moment,
 * covariance and Little's-law families and the whole E[X_p e_t] block with
 * them, which is what removes the exponential requirement and admits any
 * phase-type law. The operational pair is much looser, and is the reference's
 * own "without Markovian assumption" column.</p>
 *
 * <p>BOUND CONVENTION. QN(i,r) is the reported side of the bracket on the mean
 * number of class-r tokens in place i. TN(i,r) is the same side of the bracket
 * on the token throughput of that place, and RN follows by Little's law from
 * the two. UN(i,r) = QN(i,r) DELIBERATELY: a Place is an INF station and LINE
 * reports U = Q at an infinite server, which is what SolverCTMC and
 * {@code Solver_nc_spn_analyzer} both do on the same net. The reference's place
 * utilization 1 - P(m = 0) is a different quantity and is not this column.</p>
 *
 * <p>A TRANSITION GETS NO ROW. It is a StatefulNode and not a Station, so it
 * has no station index; the mode throughputs and enabling probabilities the LP
 * also brackets stay inside {@link Spn_lpbnd}'s return value, the same way
 * {@code Spn_metrics} keeps modeTput and modeUtil off the table.</p>
 *
 * <p>HOW TIGHT. The reference's own Table 2 measures it on a four-server
 * production line: the upper side lands 2% to 11% above simulation and the
 * lower side 30% to 40% below it, both comfortably inside the operational
 * bounds it also reports. Expect a usable upper bound and a weak lower one.</p>
 *
 * <p>Reference: Z. Liu (1998). Performance analysis of stochastic timed Petri
 * nets using linear programming approach. IEEE Transactions on Software
 * Engineering 24(11), 1014-1030.</p>
 */
public final class Solver_ba_spnlp_analyzer {

    private Solver_ba_spnlp_analyzer() {}

    public static MVAResult solver_ba_spnlp_analyzer(NetworkStruct sn, SolverOptions options) {
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

        String method = options != null && options.method != null ? options.method : "";
        boolean markovian;
        int side;
        if ("spnlp.upper".equals(method)) {
            markovian = true;
            side = 1;
        } else if ("spnlp.lower".equals(method)) {
            markovian = true;
            side = 0;
        } else if ("spnlp.op.upper".equals(method)) {
            markovian = false;
            side = 1;
        } else if ("spnlp.op.lower".equals(method)) {
            markovian = false;
            side = 0;
        } else {
            throw new RuntimeException("solver_ba_spnlp_analyzer: unknown SPN bound method '"
                    + method + "'. Valid: spnlp.upper, spnlp.lower, spnlp.op.upper, "
                    + "spnlp.op.lower.");
        }

        // ---- model gates ----
        // Every other check belongs to Spn_lpbnd, which refuses by name on the
        // mode it cannot represent. What must be decided here is only whether
        // this is a Petri net at all, and whether the places carry an embedded
        // queue the relaxation has no variable for.
        boolean hasTransition = false;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Transition) {
                hasTransition = true;
                break;
            }
        }
        if (!hasTransition) {
            throw new RuntimeException("solver_ba_spnlp_analyzer: method '" + method
                    + "' bounds a stochastic Petri net; this model has no Transition node. Use "
                    + "the queueing-network bound families, or SolverMVA/SolverNC.");
        }
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) != NodeType.Place) {
                continue;
            }
            int ist = (int) sn.nodeToStation.get(i);
            if (ist >= 0 && sn.sched != null
                    && sn.sched.get(sn.stations.get(ist)) != SchedStrategy.INF) {
                throw new RuntimeException("solver_ba_spnlp_analyzer: method '" + method
                        + "' does not support queueing places: place " + sn.nodenames.get(i)
                        + " serves under a non-INF discipline, and the relaxation carries one "
                        + "variable per (place, class) marking with no notion of an embedded "
                        + "queue.");
            }
        }

        Spn_lpbnd.SpnLpOptions lpopt = new Spn_lpbnd.SpnLpOptions();
        lpopt.markovian = markovian;
        Spn_lpbnd.SpnLpBounds bnd = Spn_lpbnd.spn_lpbnd(sn, lpopt);

        // ---- read the reported side back per place and class ----
        final int rn = bnd.nclasses;
        for (int pp = 0; pp < bnd.places.length; pp++) {
            int ist = (int) sn.nodeToStation.get(bnd.places[pp]);
            if (ist < 0 || ist >= M) {
                continue;
            }
            for (int k = 0; k < rn && k < K; k++) {
                int l = pp * rn + k;
                double q = bnd.tokens[side][l];
                double t = bnd.placeTput[side][l];
                ret.QN.set(ist, k, q);
                ret.UN.set(ist, k, q);
                ret.TN.set(ist, k, t);
                if (t > 0) {
                    ret.RN.set(ist, k, q / t);
                }
            }
        }

        for (int k = 0; k < K; k++) {
            int ref = (int) sn.refstat.get(k);
            if (ref >= 0 && ref < M) {
                ret.XN.set(0, k, ret.TN.get(ref, k));
            }
            double nk = 0;
            for (int i = 0; i < M; i++) {
                nk += ret.QN.get(i, k);
            }
            if (ret.XN.get(0, k) > 0 && nk > 0) {
                ret.CN.set(0, k, nk / ret.XN.get(0, k));
            }
        }

        ret.runtime = (System.nanoTime() - t0) / 1.0e9;
        ret.method = method;
        return ret;
    }
}
