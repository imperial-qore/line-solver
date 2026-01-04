package jline.api.mc

import jline.lang.NetworkStruct
import jline.lang.state.State
import jline.lang.state.ToMarginal
import jline.solvers.SolverOptions
import jline.solvers.ctmc.SolverCTMC
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

fun ctmc_ssg(sn: NetworkStruct, options: SolverOptions): SolverCTMC.CtmcSsgResult {
    // Get properly dimensioned cutoff matrix
    val cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses)
    val spaceResult = State.spaceGenerator(sn, cutoffMatrix, options)
    val nodeStateSpace = spaceResult.sn.space
    sn.space = nodeStateSpace
    val stateSpace = spaceResult.SS
    val stateSpaceHashed = spaceResult.SSh

    val nstateful = sn.nstateful
    val nclasses = sn.nclasses
    val sync = sn.sync
    val A = sync.size
    val stateSpaceAggr = Matrix(stateSpaceHashed.getNumRows(), stateSpaceHashed.getNumCols())
    for (a in 0..<A) {
        val stateCell = MatrixCell()
        for (s in 0..<stateSpaceHashed.getNumRows()) {
            val state = stateSpaceHashed.getRow(s)
            for (ind in 0..<sn.nnodes) {
                if (sn.isstateful.get(ind) == 1.0) {
                    val isf = sn.nodeToStateful.get(ind).toInt()
                    stateCell.set(isf, sn.space.get(sn.stateful.get(isf))!!.getRow(state.get(isf).toInt()))
                    if (sn.isstation.get(ind) == 1.0) {
                        val ist = sn.nodeToStation.get(ind).toInt()
                        val marginalResult =
                            ToMarginal.toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null)
                        val nir = marginalResult.nir

                        val startCol = ist * nclasses
                        stateSpaceAggr.expandMatrix(stateSpaceAggr.getNumRows(),
                            (ist + 1) * nclasses,
                            stateSpaceAggr.getNumNonZeros())
                        for (i in 0..<nir.getNumCols()) {
                            stateSpaceAggr.set(s, startCol + i, nir.get(i))
                        }
                    }
                }
            }
        }
    }
    return SolverCTMC.CtmcSsgResult(stateSpace, stateSpaceAggr, stateSpaceHashed, nodeStateSpace, sn)
}
/**
 * CTMC ssg algorithms
 */
@Suppress("unused")
class CtmcSsgAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}