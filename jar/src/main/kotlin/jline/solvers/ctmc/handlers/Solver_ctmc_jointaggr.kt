package jline.solvers.ctmc.handlers

import jline.api.mc.ctmc_solve_reducible
import jline.GlobalConstants
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.state.ToMarginal
import jline.solvers.SolverOptions
import jline.solvers.ctmc.ResultCTMC
import jline.solvers.ctmc.SolverCTMC
import jline.util.MatFileUtils
import jline.util.matrix.Matrix

class Solver_ctmc_jointaggr {
    companion object {
        @JvmStatic
        fun solver_ctmc_jointaggr(sn: NetworkStruct, options: SolverOptions): SolverCTMC.SolverCtmcJointResult {
            var sn = sn
            val M = sn.nstations
            val K = sn.nclasses
            var fname = ""
            val Tstart = System.nanoTime()

            val solverCTMCResult: ResultCTMC = Solver_ctmc.solver_ctmc(sn, options)
            val Q = solverCTMCResult.getQ()
            val SS = solverCTMCResult.getStateSpace()
            val SSq = solverCTMCResult.getStateSpaceAggr()
            sn = solverCTMCResult.getSn()

            if (options.keep) {
                try {
                    MatFileUtils.ensureWorkspaceDirectoryExists()
                    fname = MatFileUtils.genFilename("workspace")
                    MatFileUtils.saveCTMCWorkspace(SS, Q, SSq, fname)
                } catch (e: Exception) {
                    line_warning("solver_ctmc_jointaggr", "Could not save workspace to .mat file: %s", e.message)
                    fname = ""
                }
            }

            // Use ctmc_solve_reducible to handle reducible CTMCs (matches MATLAB)
            val (pi, _) = ctmc_solve_reducible(Q)
            for (row in 0 until pi.getNumRows()) {
                for (col in 0 until pi.getNumCols()) {
                    if (pi.get(row, col) < GlobalConstants.Zero) {
                        pi.set(row, col, 0.0)
                    }
                }
            }

            val state = sn.state
            val nvec = mutableListOf<Double>()

            // Build nvec using ToMarginal to get nir (marginal job counts), matching MATLAB
            for (i in 0 until sn.nstations) {
                val nodeIdx = sn.stationToNode.get(i).toInt()
                if (sn.isstateful.get(nodeIdx) != 0.0) {
                    val isf = sn.stationToStateful.get(i).toInt()
                    val nodeState = state[sn.stateful[isf]]
                    if (nodeState != null) {
                        // Use ToMarginal to get nir (marginal job counts per class)
                        // This matches MATLAB: [~,nir,~,~] = State.toMarginal(sn, isf, state{isf})
                        val margStats = ToMarginal.toMarginal(sn, nodeIdx, nodeState, null, null, null, null, null)
                        val nir = margStats.nir

                        // Add nir values to nvec (flattened, as in MATLAB: nvec = [nvec, nir(:)'])
                        for (col in 0 until nir.getNumCols()) {
                            nvec.add(nir.get(0, col))
                        }
                    }
                }
            }

            // Convert nvec to matrix for comparison
            val nvecMatrix = Matrix(1, nvec.size)
            for (i in nvec.indices) {
                nvecMatrix.set(0, i, nvec[i])
            }

            // Find matching rows in SSq (aggregate state space) using findrows logic
            var Pnir = 0.0
            for (row in 0 until SSq.getNumRows()) {
                val rowMatrix = SSq.getRow(row)
                // Check if nvecMatrix matches this row of SSq
                var matches = true
                if (nvecMatrix.getNumCols() != rowMatrix.getNumCols()) {
                    matches = false
                } else {
                    for (col in 0 until nvecMatrix.getNumCols()) {
                        if (Math.abs(nvecMatrix.get(0, col) - rowMatrix.get(0, col)) > 1e-10) {
                            matches = false
                            break
                        }
                    }
                }
                if (matches) {
                    Pnir += pi.get(row)
                }
            }

            val Tstop = System.nanoTime()
            val runtime = (Tstop - Tstart) / 1000000000.0

            val result = Matrix(1, 1)
            result.set(0, 0, Pnir)
            return SolverCTMC.SolverCtmcJointResult(result, runtime, fname)
        }
    }
}