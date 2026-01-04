package jline.solvers.ctmc.handlers

import jline.api.mc.ctmc_solve
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.ctmc.ResultCTMC
import jline.solvers.ctmc.SolverCTMC
import jline.util.MatFileUtils
import jline.util.matrix.Matrix

class Solver_ctmc_marg(private val solverCTMC: SolverCTMC?) {
    companion object {
        @JvmStatic
        fun solver_ctmc_marg(sn: NetworkStruct, options: SolverOptions): Matrix {
            var sn = sn
            val M = sn.nstations
            val K = sn.nclasses
            val state = sn.state
            var fname = ""
            val Tstart = System.nanoTime()
            val solverCTMCResult: ResultCTMC = Solver_ctmc.Companion.solver_ctmc(sn, options)
            val Q = solverCTMCResult.getQ()
            val SS = solverCTMCResult.getStateSpace()
            val SSq = solverCTMCResult.getStateSpaceAggr()
            sn = solverCTMCResult.getSn()
            if (options.keep) {
                try {
                    MatFileUtils.ensureWorkspaceDirectoryExists()
                    fname = MatFileUtils.genFilename("workspace")
                    MatFileUtils.saveCTMCWorkspace(SS, Q, null, fname)
                } catch (e: Exception) {
                    line_warning("solver_ctmc_marg", "Could not save workspace to .mat file: %s", e.message)
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
            val pisum = pi.sumSubMatrix(0, pi.getNumRows() - 1, 0, pi.getNumCols() - 1)
            for (row in 0..<pi.getNumRows()) {
                for (col in 0..<pi.getNumCols()) {
                    pi.set(row, col, pi.get(row, col) / pisum)
                }
            }
            val statesz: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
            for (ind in 0..<sn.nnodes) {
                if (sn.isstateful.get(ind) != 0.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    statesz.putIfAbsent(isf, sn.space.get(sn.stateful.get(isf))!!.getNumCols())
                }
            }
            val cstatesz: MutableList<Int?> = ArrayList<Int?>()
            cstatesz.add(0)
            var cumulativeSum = 0
            for (i in 0..<statesz.size) {
                cumulativeSum += statesz.get(i)!!
                cstatesz.add(cumulativeSum)
            }

            val Pnir = Matrix(1, sn.nstations)
            Pnir.zero()
            for (ind in 0..<sn.nnodes) {
                if (sn.isstateful.get(ind) != 0.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    val ist = sn.nodeToStation.get(ind).toInt()
                    val requiredlength = sn.space.get(sn.stateful.get(isf))!!.getNumCols()
                    val currentlength = sn.state.size
                    val state_i = Matrix(1, requiredlength)
                    val numZeros = requiredlength - currentlength
                    for (col in 0..<numZeros) {
                        state_i.set(0, col, 0)
                    }
                    for (col in 0..<currentlength) {
                        state_i.set(0, numZeros + col, sn.state.get(sn.stateful.get(isf))!!.get(0, col))
                    }
                    val colStart = cstatesz.get(isf)!! + 1
                    val colStop = cstatesz.get(isf)!! + state_i.getNumCols() - 1
                    val SS_slice = Matrix(SS.getNumRows(), colStop - colStart + 1)
                    Matrix.extract(SS, 0, SS.getNumRows(), colStart, colStop + 1, SS_slice, 0, 0)

                    val rowmatched = Matrix.findRows(SS, state_i)
                    var sum = 0.0
                    for (row in rowmatched) {
                        sum += pi.sumRows(row!!)
                    }
                    Pnir.set(0, ist, sum)
                }
            }
            val Tstop = System.nanoTime()
            val runtime = (Tstop - Tstart) / 1000000000.0
            return Pnir
        }
    }
}