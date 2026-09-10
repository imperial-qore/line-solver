/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.npfqn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;

/**
 * Near-immediate feedback elimination for the robust queueing network analyzer.
 *
 * <p>WHY FEEDBACK BREAKS DECOMPOSITION. A parametric decomposition treats the
 * arrival stream at each station as renewal. Feedback destroys that badly: a
 * customer that leaves a busy station and comes straight back arrives exactly
 * when the station is busy, so the flow is strongly correlated with the queue it
 * feeds. The fix is not to model the correlation but to REMOVE the feedback, by
 * folding the repeated visits into the service time:
 *
 * <ul>
 *   <li>effective mean service E[S]/(1-p),</li>
 *   <li>effective service SCV p + (1-p)cs^2 (eq. 37),</li>
 *   <li>fresh arrival rate lambda(1-p),</li>
 *   <li>per-visit waiting time (1-p) times the wait in the modified system.</li>
 * </ul>
 *
 * <p>The modified system has the SAME heavy-traffic limits for queue length,
 * workload, waiting time and external departures, so this is asymptotically
 * exact rather than merely plausible.
 *
 * <p>NEAR-IMMEDIATE, NOT JUST IMMEDIATE. What matters is whether the customer
 * returns WITHOUT PASSING A BUSIER STATION: a detour through a station of lower
 * traffic intensity is fast on the time scale of the busy station. The
 * probability computed here is the probability of returning to station i through
 * stations of strictly smaller rho only.
 *
 * <p>Port of MATLAB npfqn_feedback_elim.m.
 *
 * <p>Reference: W. Whitt, W. You (2022). A robust queueing network analyzer
 * based on indices of dispersion. Naval Research Logistics 69(1), 36-56.
 *
 * @since LINE 3.1.0
 */
public final class Npfqn_feedback_elim {

    private Npfqn_feedback_elim() {
    }

    /**
     * @param P             routing matrix, substochastic
     * @param rho           traffic intensity of each station
     * @param cs2           service SCV of each station, or null
     * @param lambda        arrival rate of each station, or null
     * @param immediateOnly keep only the self-loops, i.e. Section 4.1 feedback
     * @return map with feedbackProb, visitInflation, modifiedScv, modifiedRates
     *         (as double[]), modifiedRouting (a Matrix) and reductionExact
     */
    public static Map<String, Object> npfqn_feedback_elim(Matrix P, double[] rho, double[] cs2,
                                                           double[] lambda, boolean immediateOnly) {
        int m = P.getNumRows();
        if (P.getNumCols() != m) {
            throw new RuntimeException("npfqn_feedback_elim: the routing matrix must be square");
        }
        for (int i = 0; i < m; i++) {
            double row = 0.0;
            for (int j = 0; j < m; j++) {
                if (P.get(i, j) < -1e-12) {
                    throw new RuntimeException(
                            "npfqn_feedback_elim: the routing matrix must be non-negative");
                }
                row += P.get(i, j);
            }
            if (row > 1 + 1e-9) {
                throw new RuntimeException(
                        "npfqn_feedback_elim: the routing matrix must be substochastic");
            }
        }
        if (rho.length != m) {
            throw new RuntimeException(
                    "npfqn_feedback_elim: one traffic intensity per station is required");
        }

        double[] phat = new double[m];
        double[] inflation = new double[m];
        for (int i = 0; i < m; i++) {
            double ret = P.get(i, i);
            if (!immediateOnly) {
                // Stations a customer may pass through on a near-immediate
                // return: those NOT MORE loaded than i. A detour through a
                // busier station is not fast on the time scale of station i, so
                // it is not near-immediate; one through a station of equal load
                // is, which is why the test is <= and not <. This is the cloud
                // of eqs. (3.8)-(3.9) with H = {i}, and the same one
                // Solver_rqna applies -- the two must not drift.
                List<Integer> idx = new ArrayList<Integer>();
                for (int j = 0; j < m; j++) {
                    if (j != i && rho[j] <= rho[i] + 1e-9) {
                        idx.add(j);
                    }
                }
                if (!idx.isEmpty()) {
                    int k = idx.size();
                    Matrix A = new Matrix(k, k);
                    Matrix b = new Matrix(k, 1);
                    for (int a = 0; a < k; a++) {
                        for (int c = 0; c < k; c++) {
                            A.set(a, c, (a == c ? 1.0 : 0.0) - P.get(idx.get(a), idx.get(c)));
                        }
                        b.set(a, 0, P.get(idx.get(a), i));
                    }
                    // (I-Q)^-1 r: the probability of eventually reaching i from
                    // each allowed station without leaving the allowed set.
                    Matrix reach = new Matrix(k, 1);
                    Matrix.solve(A, b, reach);
                    for (int a = 0; a < k; a++) {
                        ret += P.get(i, idx.get(a)) * reach.get(a, 0);
                    }
                }
            }
            phat[i] = Math.min(Math.max(ret, 0.0), 1.0 - 1e-12);
            inflation[i] = 1.0 / (1.0 - phat[i]);
        }

        Map<String, Object> res = new HashMap<String, Object>();
        res.put("feedbackProb", phat);
        res.put("visitInflation", inflation);
        if (cs2 != null) {
            if (cs2.length != m) {
                throw new RuntimeException(
                        "npfqn_feedback_elim: one service SCV per station is required");
            }
            double[] mod = new double[m];
            for (int i = 0; i < m; i++) {
                mod[i] = phat[i] + (1.0 - phat[i]) * cs2[i];         // eq. (37)
            }
            res.put("modifiedScv", mod);
        }
        if (lambda != null) {
            if (lambda.length != m) {
                throw new RuntimeException(
                        "npfqn_feedback_elim: one arrival rate per station is required");
            }
            double[] mod = new double[m];
            for (int i = 0; i < m; i++) {
                mod[i] = lambda[i] * (1.0 - phat[i]);
            }
            res.put("modifiedRates", mod);
        }

        // The immediate-feedback reduction: drop the self-loop and renormalize
        // the rest of the row. For near-immediate feedback the return path runs
        // through other stations, so no row-local reduction exists and the
        // elimination applies to the service description instead.
        Matrix Pmod = P.copy();
        boolean exact = true;
        for (int i = 0; i < m; i++) {
            if (Math.abs(phat[i] - P.get(i, i)) > 1e-12) {
                exact = false;
            }
            double loop = P.get(i, i);
            if (loop <= 0) {
                continue;
            }
            Pmod.set(i, i, 0.0);
            double rest = 0.0;
            double total = 0.0;
            for (int j = 0; j < m; j++) {
                rest += Pmod.get(i, j);
                total += P.get(i, j);
            }
            if (rest > 0) {
                for (int j = 0; j < m; j++) {
                    Pmod.set(i, j, Pmod.get(i, j) * (total - loop) / rest);
                }
            }
        }
        res.put("modifiedRouting", Pmod);
        res.put("reductionExact", immediateOnly || exact);
        return res;
    }
}
