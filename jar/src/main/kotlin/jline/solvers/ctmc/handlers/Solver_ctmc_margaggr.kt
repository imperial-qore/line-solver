package jline.solvers.ctmc.handlers

import jline.api.mc.ctmc_solve_reducible
import jline.GlobalConstants
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.state.State
import jline.lang.state.ToMarginal
import jline.solvers.SolverOptions
import jline.solvers.ctmc.ResultCTMC
import jline.solvers.ctmc.ResultCTMCMargAggr
import jline.util.MatFileUtils
import jline.util.matrix.Matrix

class Solver_ctmc_margaggr {
    companion object {
        @JvmStatic
        fun solver_ctmc_margaggr(sn: NetworkStruct, options: SolverOptions): ResultCTMCMargAggr {
            var sn = sn
            val M = sn.nstations
            val K = sn.nclasses
            val state = sn.state
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
                    println("\nCTMC generator and state space saved in: ${fname}.mat")
                } catch (e: Exception) {
                    line_warning("solver_ctmc_margaggr", "Could not save workspace to .mat file: %s", e.message)
                    fname = ""
                }
            }

            // Use ctmc_solve_reducible to handle reducible CTMCs
            val (pi, _) = ctmc_solve_reducible(Q)
            for (row in 0 until pi.getNumRows()) {
                for (col in 0 until pi.getNumCols()) {
                    if (pi.get(row, col) < GlobalConstants.Zero) {
                        pi.set(row, col, 0.0)
                    }
                }
            }

            val statesz = mutableListOf<Int>()
            for (ind in 0 until sn.nnodes) {
                if (sn.isstateful.get(ind) != 0.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    statesz.add(sn.space[sn.stateful[isf]]!!.getNumCols())
                }
            }

            val cstatesz = mutableListOf<Int>()
            cstatesz.add(0)
            var cumulativeSum = 0
            for (size in statesz) {
                cumulativeSum += size
                cstatesz.add(cumulativeSum)
            }

            val Pnir = Matrix(1, sn.nstations)
            Pnir.zero()

            for (ind in 0 until sn.nnodes) {
                if (sn.isstateful.get(ind) != 0.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    val ist = sn.nodeToStation.get(ind).toInt()
                    
                    val nodeState = state[sn.stateful[isf]]
                    if (nodeState != null) {
                        val stateLength = sn.space[sn.stateful[isf]]!!.getNumCols()
                        val state_i = Matrix(1, stateLength)
                        
                        // Pad with zeros if needed (matches MATLAB line 39)
                        val numZeros = stateLength - nodeState.getNumCols()
                        for (col in 0 until numZeros) {
                            state_i.set(0, col, 0.0)
                        }
                        for (col in 0 until nodeState.getNumCols()) {
                            state_i.set(0, numZeros + col, nodeState.get(0, col))
                        }

                        // Use State.toMarginal to get the marginal statistics (matches MATLAB line 40)
                        // MATLAB passes state{isf} directly, not the padded state_i
                        val margStats = ToMarginal.toMarginal(sn, ind, nodeState, null, null, null, null, null)
                        val nivec = margStats.nir

                        Pnir.set(0, ist, 0.0)
                        for (s in 0 until SS.getNumRows()) {
                            val colStart = cstatesz[isf]
                            val colEnd = cstatesz[isf] + stateLength
                            val stateSlice = Matrix.extract(SS, s, s + 1, colStart, colEnd)
                            
                            // Use State.toMarginal for the state from SS (matches MATLAB line 43)
                            val sivecStats = ToMarginal.toMarginal(sn, ind, stateSlice, null, null, null, null, null)
                            val sivec = sivecStats.nir
                            
                            // Check if all elements of sivec match nivec (matches MATLAB line 44)
                            var matches = true
                            if (sivec.getNumCols() != nivec.getNumCols()) {
                                matches = false
                            } else {
                                for (col in 0 until sivec.getNumCols()) {
                                    if (Math.abs(sivec.get(0, col) - nivec.get(0, col)) > 1e-10) {
                                        matches = false
                                        break
                                    }
                                }
                            }
                            
                            if (matches) {
                                Pnir.set(0, ist, Pnir.get(0, ist) + pi.get(s))
                            }
                        }
                    }
                }
            }

            val Tstop = System.nanoTime()
            val runtime = (Tstop - Tstart) / 1000000000.0
            
            return ResultCTMCMargAggr(Pnir, pi, runtime, fname)
        }
    }
}