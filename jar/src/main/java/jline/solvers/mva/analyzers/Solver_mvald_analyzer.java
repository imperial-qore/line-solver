/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.handlers.Solver_amva;
import jline.solvers.mva.handlers.Solver_mvald;

public final class Solver_mvald_analyzer {
    private Solver_mvald_analyzer() {}

    /**
     * MVALD Analyzer.
     */
    public static MVAResult solver_mvald_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        String method = options.method;
        MVAResult ret = null;
        if ("exact".equals(method) || "mva".equals(method)) {
            if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
                throw new RuntimeException("Exact class-dependent solver not available in MVA.");
            }
            if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
                throw new RuntimeException("Exact joint-dependent solver not available in MVA.");
            }
            ret = Solver_mvald.solver_mvald(sn, options);
        } else if ("default".equals(method) || "amva".equals(method) || "qd".equals(method)
                || "lin".equals(method) || "qdlin".equals(method)) {
            if ("default".equals(method) && qualifiesForExactLd(sn)) {
                // see _kb/06-solver-catalog.md for rationale
                ret = Solver_mvald.solver_mvald(sn, options);
                method = "exact";
            } else {
                ret = Solver_amva.solver_amva(sn, options);
            }
        } else {
            throw new RuntimeException("The " + method + " method is not supported by the load-dependent MVA solver.");
        }
        long endTime = System.nanoTime();
        res.QN = ret.QN;
        res.UN = ret.UN;
        res.RN = ret.RN;
        res.TN = ret.TN;
        res.CN = ret.CN;
        res.XN = ret.XN;
        res.AN = ret.AN;
        res.WN = ret.WN;
        res.logNormConstAggr = ret.logNormConstAggr;
        res.runtime = (endTime - startTime) / 1000000000.0;
        res.iter = ret.iter;
        res.method = method;
        return res;
    }

    /**
     * True when the LD model qualifies for the exact load-dependent MVA under
     * the default method: closed-only, no class-dependent scaling, at most 4
     * chains and 20 jobs, product-form, integral populations.
     */
    private static boolean qualifiesForExactLd(NetworkStruct sn) {
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            return false;
        }
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            return false;
        }
        double totJobs = 0.0;
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            double n = sn.njobs.get(r);
            if (!Double.isFinite(n) || n != Math.floor(n)) {
                return false;
            }
            totJobs += n;
        }
        if (sn.nchains > 4 || totJobs > 20) {
            return false;
        }
        return jline.api.sn.SnHasProductForm.snHasProductForm(sn);
    }
}
