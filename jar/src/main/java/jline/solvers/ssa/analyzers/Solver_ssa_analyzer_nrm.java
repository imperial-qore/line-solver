package jline.solvers.ssa.analyzers;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.StatefulNode;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SSAResult;
import jline.solvers.ssa.handlers.Solver_ssa_nrm;
import jline.solvers.ssa.handlers.SolverSSAResultNRM;
import jline.util.matrix.Matrix;

public final class Solver_ssa_analyzer_nrm {
    private Solver_ssa_analyzer_nrm() {}

    public static SSAResult solver_ssa_analyzer_nrm(
            NetworkStruct sn,
            Map<StatefulNode, Matrix> init_state,
            SolverOptions options) {
        // see _kb/06-solver-catalog.md for rationale
        Set<SchedStrategy> allowedSched = new HashSet<SchedStrategy>(Arrays.asList(
                SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
                SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS,
                SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO,
                SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT,
                SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING));
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy schedPolicy = sn.sched.get(sn.stations.get(i));
            if (!allowedSched.contains(schedPolicy)) {
                throw new RuntimeException("solver_ssa_analyzer_nrm:UnsupportedPolicy - the NRM method does not support the scheduling policy: " + schedPolicy);
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        if (!Solver_ssa_analyzer.phaseNrmOK(sn)) {
            throw new RuntimeException("solver_ssa_analyzer_nrm:UnsupportedPhaseService - the NRM expands phase-type service only at INF/PS and non-preemptive buffered stations; use method='serial' for phase-type service at LCFSPR or POLLING stations.");
        }

        // Check state space generation configuration
        String stateSpaceGen = "default";
        if (options.config != null && options.config.state_space_gen != null) {
            stateSpaceGen = options.config.state_space_gen;
        }

        if ("full".equals(stateSpaceGen)) {
            // Use the explicit state space generation version
            return Solver_ssa_analyzer_nrm_space.solver_ssa_analyzer_nrm_space(sn, options);
        } else {
            // Use the default direct metric computation version
            SolverSSAResultNRM result = Solver_ssa_nrm.solver_ssa_nrm(sn, options);
            return new SSAResult(result.getQN(), result.getUN(), result.getRN(), result.getTN(),
                    result.getCN(), result.getXN(), null, null, sn);
        }
    }
}
