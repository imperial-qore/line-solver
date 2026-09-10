/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import jline.api.spn.Spn_metrics;
import jline.api.spn.Spn_pf;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

/**
 * Stationary analysis of a PRODUCT-FORM stochastic Petri net by MDD-rec: the
 * normalising constant is obtained from one memoised walk of the decision
 * diagram holding the reachable set, and every reported measure is a masked walk
 * of the same diagram.
 *
 * <p>This is the "rec" method of SolverNC, and the first analytical route LINE
 * offers for a Petri net -- CTMC solves the explicit generator, SSA and LDES
 * simulate, FLD fluidises. Three functions do the work and each is the subject
 * of its own reference. {@link Spn_pf} decides the product form and derives the
 * per-place factors g_l (Coleman-Henderson-Taylor complex balance).
 * {@link jline.api.mdd.Mdd_rec} evaluates G = sum_S prod_l g_l(s_l) in
 * O(sum_l nodes_l * |S_l|) rather than O(|S|) (Balsamo-Marin-Stojic).
 * {@link Spn_metrics} reads the mean tokens, the place and mode utilisation and
 * the throughputs off masked walks of the same diagram.</p>
 *
 * <p><b>What this reaches that the explicit generator does not.</b> The diagram
 * stores the reachable set, never the generator, so the cost is set by the number
 * of diagram nodes and not by |S|. It also does not need the marking to be a
 * conserved job population: a mode may consume two tokens and produce one, or
 * consume one and produce two, which is the fork-join and batch case that the
 * MDD-rec paper exists to serve.</p>
 *
 * <p><b>UN follows LINE, not the paper.</b> A Place is an INF station, and LINE
 * reports U = Q at an infinite server, which is what SolverCTMC returns for the
 * same net. The paper's place utilisation u(P_j) = 1 - P(m_j = 0) is a different
 * quantity and is reported separately, on the certificate's metrics block.</p>
 */
public final class Solver_nc_spn_analyzer {

    private Solver_nc_spn_analyzer() {
    }

    /**
     * Analyse a product-form stochastic Petri net.
     *
     * @param model   the Network holding the Places and Transitions
     * @param sn      its network structure
     * @param options solver options; tol and verbose are read
     * @return QN, UN, RN and TN per (Place station, class), CN and XN per class,
     *         and lG the log normalising constant
     */
    public static NCResult solver_nc_spn_analyzer(Network model, NetworkStruct sn,
                                                  SolverOptions options) {
        long Tstart = System.nanoTime();

        Spn_pf.SpnPfOptions pfopt = new Spn_pf.SpnPfOptions();
        if (options != null) {
            pfopt.verbose = options.verbose != null && options.verbose.ordinal() > 1;
            if (options.tol > 0) {
                pfopt.tol = Math.max(options.tol, 1e-12);
            }
        }
        Spn_pf.SpnPfResult pf = Spn_pf.spn_pf(model, pfopt);
        Spn_metrics.SpnMetricsResult met =
                Spn_metrics.spn_metrics(pf.spn.mdds, pf.g, pf.spn.info);
        pf.metrics = met;

        int M = sn.nstations;
        int R = sn.nclasses;
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix Rt = new Matrix(M, R);
        Matrix T = new Matrix(M, R);
        Q.fill(0.0);
        U.fill(0.0);
        Rt.fill(0.0);
        T.fill(0.0);

        int[] places = pf.spn.info.places;
        for (int pp = 0; pp < places.length; pp++) {
            int ist = (int) sn.nodeToStation.get(places[pp]);
            if (ist < 0) {
                continue;
            }
            for (int k = 0; k < R; k++) {
                int l = pp * R + k;
                Q.set(ist, k, met.tokens[l]);
                // INF station: LINE charges one server per resident token, so U = Q.
                // The paper's 1 - P(m = 0) is met.placeUtil, on the certificate.
                U.set(ist, k, met.tokens[l]);
                T.set(ist, k, met.placeTput[l]);
                if (met.placeTput[l] > 0) {
                    Rt.set(ist, k, met.tokens[l] / met.placeTput[l]);   // Little at the place
                }
            }
        }

        // System throughput at the reference station of each class, and the
        // response time Little's law then fixes. A net whose class population is
        // not conserved has no meaningful N/X, so CN stays zero there rather than
        // reporting a ratio against a moving population.
        Matrix X = new Matrix(1, R);
        Matrix CN = new Matrix(1, R);
        X.fill(0.0);
        CN.fill(0.0);
        for (int k = 0; k < R; k++) {
            int ref = (int) sn.refstat.get(k);
            if (ref >= 0 && ref < M) {
                X.set(0, k, T.get(ref, k));
            }
            double Nk = 0;
            for (int i = 0; i < M; i++) {
                Nk += Q.get(i, k);
            }
            if (X.get(0, k) > 0 && Nk > 0) {
                CN.set(0, k, Nk / X.get(0, k));
            }
        }

        NCResult res = new NCResult();
        res.QN = Q;
        res.UN = U;
        res.RN = Rt;
        res.TN = T;
        res.XN = X;
        res.CN = CN;
        res.lG = Math.log(met.G);
        res.it = 1;
        res.iter = 1;
        res.method = "rec";
        res.spnpf = pf;
        res.runtime = (double) (System.nanoTime() - Tstart) / 1000000000.0;
        return res;
    }
}
