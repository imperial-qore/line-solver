package jline.solvers.ssa.analyzers

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.StatefulNode
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAResult
import jline.solvers.ssa.handlers.solver_ssa_nrm
import jline.util.matrix.Matrix

fun solver_ssa_analyzer_nrm(
    sn: NetworkStruct,
    init_state: MutableMap<StatefulNode?, Matrix?>,
    options: SolverOptions,
): SSAResult {
    val allowedSched = setOf(SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
        SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO)
    for (i in 0 until sn.nstations) {
        val schedPolicy = sn.sched.get(sn.stations.get(i))
        if (schedPolicy !in allowedSched) {
            throw RuntimeException("solver_ssa_analyzer_nrm:UnsupportedPolicy - Only INF, EXT, PS, DPS, GPS, PSPRIO, DPSPRIO, GPSPRIO are supported.")
        }
    }

    // Check state space generation configuration
    val stateSpaceGen = options.config?.state_space_gen ?: "default"

    return if (stateSpaceGen == "full") {
        // Use the explicit state space generation version
        solver_ssa_analyzer_nrm_space(sn, options)
    } else {
        // Use the default direct metric computation version
        val result = solver_ssa_nrm(sn, options)
        SSAResult(result.QN, result.UN, result.RN, result.TN, result.CN, result.XN, null, null, sn)
    }
}
