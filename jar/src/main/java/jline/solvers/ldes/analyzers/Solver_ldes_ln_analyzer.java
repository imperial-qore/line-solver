/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ldes.analyzers;

import java.util.List;

import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LNLDESResult;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.ldes.handlers.Solver_ssj_ln;

public final class Solver_ldes_ln_analyzer {
    private Solver_ldes_ln_analyzer() {}

    /**
     * LN LDES analyzer - validates LQN structure and dispatches to SSJ backend.
     */
    public static LNLDESResult solver_ldes_ln_analyzer(
            LayeredNetworkStruct lsn,
            SolverOptions options,
            SolverLDES solverLDES) {
        long Tstart = System.nanoTime();

        validateLQNStructure(lsn);

        String backendMethod = options.method;
        if (backendMethod == null || "default".equals(backendMethod)) {
            backendMethod = "ssj";
        }

        LNLDESResult result;
        if ("ssj".equals(backendMethod)) {
            result = Solver_ssj_ln.solver_ssj_ln(lsn, options);
        } else {
            throw new RuntimeException("solver_ldes_ln_analyzer: Unknown method '" + backendMethod + "'");
        }

        result.method = options.method != null ? options.method : "default";
        result.runtime = (double) (System.nanoTime() - Tstart) / 1e9;

        return result;
    }

    private static void validateLQNStructure(LayeredNetworkStruct lsn) {
        if (lsn.nhosts == 0) {
            throw new RuntimeException("solver_ldes_ln_analyzer: LQN must have at least one Host/Processor");
        }
        if (lsn.ntasks == 0) {
            throw new RuntimeException("solver_ldes_ln_analyzer: LQN must have at least one Task");
        }

        for (int t = 1; t <= lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            List<Integer> entries = lsn.entriesof.get(tidx);
            if (entries == null || entries.isEmpty()) {
                String taskName = lsn.names.get(tidx);
                if (taskName == null) taskName = "Task" + tidx;
                throw new RuntimeException("solver_ldes_ln_analyzer: Task '" + taskName + "' has no entries");
            }
        }

        for (int e = 1; e <= lsn.nentries; e++) {
            int eidx = lsn.eshift + e;
            int tidx = (int) lsn.parent.get(0, eidx);
            List<Integer> activities = lsn.actsof.get(tidx);
            if (activities == null || activities.isEmpty()) {
                String entryName = lsn.names.get(eidx);
                if (entryName == null) entryName = "Entry" + eidx;
                throw new RuntimeException("solver_ldes_ln_analyzer: Entry '" + entryName + "' has no activities in its task");
            }
        }
    }
}
