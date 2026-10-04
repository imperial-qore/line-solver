/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag.analyzers;

import jline.api.sn.SnNonmarkovToPh;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.ag.AGResult;
import jline.solvers.ag.handlers.Solver_ag;
import jline.util.matrix.Matrix;

/**
 * Analyzer of the agent-based (RCAT) solver.
 *
 * <p>Short by design: the whole analysis is the reversed-rate fixed point in
 * {@link Solver_ag}, and what remains here is the conversion every RCAT method
 * needs before it and the tidy-up every solver needs after it.</p>
 */
public final class Solver_ag_analyzer {

    private Solver_ag_analyzer() {}

    public static AGResult solver_ag_analyzer(NetworkStruct snInput, SolverOptions options) {
        double start = System.nanoTime();

        // RCAT builds a CTMC per agent out of (D0,D1): a preserved Det would
        // reach it with no matrix at all and be read back as its mean rate, and a
        // concentrated matrix exponential is not a generator at all. Unlike the
        // MAM analyzer, which decides this per method, EVERY AG method needs the
        // genuine phase-type, so the choice is not method-dependent here.
        options.config.preserveDet = false;
        options.config.phfit = "ph";
        NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(snInput, options);

        String method = (options.method == null || options.method.isEmpty())
                ? "default" : options.method;
        InputOutput.line_debug(options.verbose, "AG analyzer starting: method=" + method);

        AGResult result = Solver_ag.solver_ag(sn, options);

        // 'default' resolves to inap, and 'exact' falls back to it with a warning
        // (AutoCAT is unreachable). Report what actually ran: 'exact' is
        // classified globally as an exact method, so leaving the name in place
        // would banner an iterative approximation as exact.
        if ("default".equals(method) || "exact".equals(method)) {
            result.method = "inap";
        } else {
            result.method = method;
        }

        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                for (int r = 0; r < sn.nclasses; r++) {
                    result.TN.set(i, r, sn.rates.get(i, r));
                }
            }
        }

        zeroNaN(result.QN);
        zeroNaN(result.UN);
        zeroNaN(result.RN);
        zeroNaN(result.TN);
        zeroNaN(result.CN);
        zeroNaN(result.XN);

        result.runtime = (System.nanoTime() - start) / 1.0e9;
        return result;
    }

    private static void zeroNaN(Matrix m) {
        if (m == null) return;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (Double.isNaN(m.get(i, j))) m.set(i, j, 0.0);
            }
        }
    }
}
