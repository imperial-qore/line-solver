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
            for (i in 0..<sn.nstations) {
                if (sn.isstateful.get(0, i) != 0.0) {
                    val isf = sn.nodeToStateful.get(0, i).toInt()
                    val requiredlength = sn.space.get(sn.stateful.get(isf))!!.getNumCols()
                    val currentlength = sn.state.size
                    val state_i = Matrix(1, requiredlength)
                    val numZeros = requiredlength - currentlength
                    for (col in 0..<numZeros) {
                        state_i.set(0, col, 0)
                    }
                    for (col in 0..<currentlength) {
                        state_i.set(0, numZeros + col, sn.state.get(sn.stateful.get(isf))!!.get(col))
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
            val cols: MutableList<Int?> = ArrayList<Int?>()
            cols.add(0)
            val Pnir = Matrix(PnirIndex.size, cols.size)
            for (row in PnirIndex.indices) {
                val rowIndex: Int = PnirIndex.get(row)!!
                for (col in cols.indices) {
                    val colIndex: Int = cols.get(col)!!
                    Matrix.extract(pi, rowIndex, rowIndex + 1, colIndex, colIndex + 1, Pnir, row, col)
                }
            }

            val runtime = (System.nanoTime() - T0) / 1000000000.0
            return SolverCTMC.SolverCtmcJointResult(Pnir, runtime, fname)
        }
    }
}