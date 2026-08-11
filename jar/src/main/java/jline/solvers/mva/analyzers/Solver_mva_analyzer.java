/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import java.util.Map;

import jline.VerboseLevel;
import jline.api.sn.SnHasFractionalPopulations;
import jline.api.sn.SnHasProductForm;
import jline.api.sn.SnHasBurstyArrival;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.handlers.Solver_amva;
import jline.solvers.mva.handlers.Solver_sqd;
import jline.solvers.mva.handlers.Solver_mva;
import jline.solvers.mva.handlers.Solver_mva_sum;
import jline.solvers.mva.handlers.Solver_qna;
import jline.solvers.mva.handlers.Solver_rqna;

public final class Solver_mva_analyzer {
    private Solver_mva_analyzer() {}

    /**
     * Detects a closed single-chain network with Blocking-After-Service (BAS) finite-buffer
     * blocking, which the {@link Solver_sqd} approximation handles but exact/AMVA MVA does not.
     */
    public static boolean isBasModel(NetworkStruct sn) {
        if (sn.nchains != 1 || sn.nclosedjobs <= 0 || sn.droprule == null) {
            return false;
        }
        for (int r = 0; r < sn.nclasses; r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                return false; // open class present
            }
        }
        for (Station st : sn.stations) {
            Map<JobClass, DropStrategy> rules = sn.droprule.get(st);
            if (rules != null) {
                for (DropStrategy ds : rules.values()) {
                    if (ds == DropStrategy.BlockingAfterService) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * MVA Analyzer.
     */
    public static MVAResult solver_mva_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        String method = options.method;
        method = method.replaceFirst("^amva\\.", "");

        long startTime = System.nanoTime();
        MVAResult ret = null;
        if ("exact".equals(method) || "mva".equals(method)) {
            ret = Solver_mva.solver_mva(sn, options);
            ret.iter = 0;
        } else if ("mvac".equals(method)) {
            ret = Solver_mva.solver_mvac(sn, options);
            ret.iter = 0;
        } else if ("sqd".equals(method)) {
            ret = Solver_sqd.solver_sqd(sn, options);
        } else if ("sum".equals(method) || "esum".equals(method)) {
            ret = Solver_mva_sum.solver_mva_sum(sn, options);
        } else if ("qna".equals(method)) {
            ret = Solver_qna.solver_qna(sn, options);
        } else if ("rqna".equals(method)) {
            ret = Solver_rqna.solver_rqna(sn, options);
        } else if ("default".equals(method)) {
            // see _kb/06-solver-catalog.md for rationale
            boolean allOpen = true;
            for (int r = 0; r < sn.nclasses; r++) {
                if (!Double.isInfinite(sn.njobs.get(r))) { allOpen = false; break; }
            }
            if (sn.nclasses == 1 && allOpen && SnHasBurstyArrival.snHasBurstyArrival(sn)) {
                ret = Solver_rqna.solver_rqna(sn, options);
                method = "rqna";
            } else {
            // see _kb/06-solver-catalog.md for rationale
            boolean hasOpenClass = false, hasClosedClass = false;
            for (int r = 0; r < sn.nclasses; r++) {
                double nj = sn.njobs.get(r);
                if (Double.isInfinite(nj)) hasOpenClass = true;
                else if (nj > 0) hasClosedClass = true;
            }
            // see _kb/06-solver-catalog.md for rationale
            double maxFiniteServers = Double.NaN;
            boolean hasFiniteServer = false;
            for (int i = 0; i < sn.nstations; i++) {
                double si = sn.nservers.get(i);
                if (!Double.isInfinite(si)) {
                    hasFiniteServer = true;
                    if (Double.isNaN(maxFiniteServers) || si > maxFiniteServers) maxFiniteServers = si;
                }
            }
            // see _kb/06-solver-catalog.md for rationale
            boolean closedPopsIntegral = true;
            for (int r = 0; r < sn.nclasses; r++) {
                double nj = sn.njobs.get(r);
                if (!Double.isInfinite(nj) && nj != Math.floor(nj)) {
                    closedPopsIntegral = false;
                    break;
                }
            }
            if (isBasModel(sn)) {
                ret = Solver_sqd.solver_sqd(sn, options);
                method = "sqd";
            } else if (hasOpenClass && hasClosedClass && hasFiniteServer && maxFiniteServers == 1.0
                    && SnHasProductForm.snHasProductForm(sn)
                    && closedPopsIntegral) {
                ret = Solver_mva.solver_mva(sn, options);
                method = "exact";
            } else if (sn.nchains <= 4 && sn.njobs.sumRows().toDouble() <= 20
                    && SnHasProductForm.snHasProductForm(sn)
                    && !SnHasFractionalPopulations.snHasFractionalPopulations(sn)) {
                ret = Solver_mva.solver_mva(sn, options);
                method = "exact";
            } else {
                ret = Solver_amva.solver_amva(sn, options);
                method = ret.method;
            }
            }
        } else if ("amva".equals(method) || "bs".equals(method) || "qd".equals(method)
                || "qli".equals(method) || "fli".equals(method) || "lin".equals(method)
                || "qdlin".equals(method) || "sqni".equals(method) || "gflin".equals(method)
                || "egflin".equals(method) || "schmidt".equals(method) || "schmidt-ext".equals(method)
                || "ab".equals(method)) {
            ret = Solver_amva.solver_amva(sn, options);
            method = ret.method;
        } else {
            // An unsupported method must FAIL by name. Returning a result whose
            // every metric is null left the caller to dereference them, so the
            // failure surfaced as a NullPointerException naming neither the
            // method nor this analyzer, and with verbose SILENT there was no
            // message at all. Matches the line_error in
            // matlab/src/solvers/MVA/solver_mva_analyzer.m.
            throw new RuntimeException("solver_mva_analyzer: the '" + method
                    + "' method is not dispatched by solver_mva_analyzer. Supported: default, "
                    + "exact, mva, mvac, amva, bs, qd, qli, fli, lin, qdlin, sqni, egflin, gflin, ab, "
                    + "schmidt, schmidt-ext.");
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
        res.iter = ret.iter;
        res.method = method;
        res.runtime = (endTime - startTime) / 1000000000.0;

        // see _kb/06-solver-catalog.md for rationale
        res.converged = ret.converged;
        boolean nonConverged;
        if (ret.converged != null) {
            nonConverged = !ret.converged.booleanValue();
        } else {
            // Single-loop handlers report no flag; there the count is sound.
            nonConverged = options.iter_max > 0 && res.iter >= options.iter_max;
        }
        if (nonConverged) {
            InputOutput.line_warning_always("solver_mva_analyzer",
                    String.format("AMVA method '%s' did not meet the convergence tolerance %g after %d iterations; the returned metrics may not be converged. Try another method (e.g. 'qd' or 'bs'), raise options.iter_max, or loosen options.iter_tol.",
                            method, options.iter_tol, res.iter));
        }
        return res;
    }
}
