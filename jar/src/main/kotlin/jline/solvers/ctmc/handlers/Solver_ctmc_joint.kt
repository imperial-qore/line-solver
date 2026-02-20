package jline.solvers.ctmc.handlers

import jline.api.mc.ctmc_solve
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.ctmc.ResultCTMC
import jline.solvers.ctmc.SolverCTMC
import jline.util.MatFileUtils
import jline.util.matrix.Matrix

class Solver_ctmc_joint(private val solverCTMC: SolverCTMC?) {
    companion object {
        @JvmStatic
        fun solver_ctmc_joint(sn: NetworkStruct, options: SolverOptions): SolverCTMC.SolverCtmcJointResult {
            val M = sn.nstations
            val K = sn.nclasses
            var fname = ""
            val T0 = System.nanoTime()
            val solverCTMCResult: ResultCTMC = Solver_ctmc.Companion.solver_ctmc(sn, options)
            val SS = solverCTMCResult.getStateSpace()
            val Q = solverCTMCResult.getQ()
            if (options.keep) {
                try {
                    MatFileUtils.ensureWorkspaceDirectoryExists()
                    fname = MatFileUtils.genFilename("workspace")
                    MatFileUtils.saveCTMCWorkspace(SS, Q, null, fname)
                } catch (e: Exception) {
                    line_warning("solver_ctmc_joint", "Could not save workspace to .mat file: %s", e.message)
                    fname = ""
                }
            }
            val pi = ctmc_solve(Q)
            for (row in 0..<pi.getNumRows()) {
                for (col in 0..<pi.getNumCols()) {
                    if (pi.get(row, col) < 0) {
                        pi.set(row, col, 0)
                    }
                }
            }
            val statevecList: MutableList<Matrix> = ArrayList<Matrix>()
            val state = sn.state
            // Iterate over nodes (not stations) to match MATLAB solver_ctmc_joint.m
            for (ind in 0..<sn.nnodes) {
                if (sn.isstateful.get(ind) != 0.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    val requiredlength = sn.space.get(sn.stateful.get(isf))!!.getNumCols()
                    // Get the state matrix for this specific stateful node
                    val stateMatrix = sn.state.get(sn.stateful.get(isf))
                    val currentlength = stateMatrix?.length() ?: 0
                    val state_i = Matrix(1, requiredlength)
                    val numZeros = requiredlength - currentlength
                    for (col in 0..<numZeros) {
                        state_i.set(0, col, 0)
                    }
                    for (col in 0..<currentlength) {
                        state_i.set(0, numZeros + col, stateMatrix!!.get(0, col))
                    }
                    statevecList.add(state_i)
                }
            }
            val totalCols = statevecList.stream().mapToInt { obj: Matrix? -> obj!!.getNumCols() }.sum()
            val statevec = Matrix(1, totalCols)
            var currentCol = 0
            for (matrix in statevecList) {
                for (col in 0..<matrix.getNumCols()) {
                    statevec.set(0, currentCol++, matrix.get(0, col))
                }
            }
            val PnirIndex = Matrix.findRows(SS, statevec)
            // pi is a row vector (1 x n), so we index by column using the state space row indices
            val Pnir = Matrix(PnirIndex.size, 1)
            for (i in PnirIndex.indices) {
                val stateIdx: Int = PnirIndex.get(i)!!
                // Get probability from pi at the state index (pi is 1 x n, so use column index)
                Pnir.set(i, 0, pi.get(0, stateIdx))
            }

            val runtime = (System.nanoTime() - T0) / 1000000000.0
            return SolverCTMC.SolverCtmcJointResult(Pnir, runtime, fname)
        }
    }
}