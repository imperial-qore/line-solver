package jline.api.mc

import jline.lang.NetworkStruct
import jline.lang.nodes.StatefulNode
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix

/**
 * CTMC State Space Generator for Reachability Analysis
 * Generates state spaces using reachability analysis for CTMC models.
 * This is used for models where not all states need to be generated, only reachable ones.
 *
 * @param sn - Network structure
 * @param options - Solver options
 * @return CtmcSsgReachabilityResult containing state spaces and updated network structure
 */
fun ctmc_ssg_reachability(sn: NetworkStruct, options: SolverOptions): CtmcSsgReachabilityResult {
    // Note: This implementation uses the regular spaceGenerator instead of a reachability-based
    // state space generator. In MATLAB, this would use reachableSpaceGenerator which only
    // generates states reachable from the initial state, potentially reducing state space size.
    // The current implementation generates all possible states but still produces correct results.
    // Get properly dimensioned cutoff matrix
    val cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses)
    val spaceResult = State.spaceGenerator(sn, cutoffMatrix, options)

    val stateSpace = spaceResult.SS
    val stateSpaceHashed = spaceResult.SSh
    val nodeStateSpace = spaceResult.sn.space
    sn.space = nodeStateSpace
    sn.spaceHash = spaceResult.ST.spaceHash

    // Ensure hide_immediate is enabled for reachability analysis
    options.config.hide_immediate = true

    sn.nstateful
    val nclasses = sn.nclasses
    val sync = sn.sync
    val A = sync.size
    val stateSpaceAggr = Matrix(stateSpaceHashed.numRows, stateSpaceHashed.numCols)
    stateSpaceAggr.zero()

    // For all synchronizations
    for (a in 0..<A) {
        val stateCell = mutableMapOf<Int, Matrix>()
        for (s in 0..<stateSpaceHashed.numRows) {
            val state = stateSpaceHashed.getRow(s)

            // Update state cell array and aggregate state space
            for (ind in 0..<sn.nnodes) {
                if (sn.isstateful[ind] > 0) {
                    val isf = sn.nodeToStateful[ind].toInt()
                    stateCell[isf] = sn.space[sn.stateful[isf]]!!.getRow(state[isf].toInt())

                    if (sn.isstation[ind] > 0) {
                        val ist = sn.nodeToStation[ind].toInt()
                        val marginalResult = toMarginal(sn, ind, stateCell[isf], null, null, null, null, null)
                        val nir = marginalResult.nir

                        val startCol = ist * nclasses
                        val endCol = (ist + 1) * nclasses

                        // Ensure stateSpaceAggr has enough columns
                        if (stateSpaceAggr.numCols < endCol) {
                            stateSpaceAggr.expandMatrix(stateSpaceAggr.numRows, endCol, stateSpaceAggr.numNonZeros)
                        }

                        for (c in 0..<nir.length()) {
                            stateSpaceAggr[s, startCol + c] = nir[c]
                        }
                    }
                }
            }
        }
    }

    return CtmcSsgReachabilityResult(stateSpace, stateSpaceAggr, stateSpaceHashed, nodeStateSpace, sn)
}

/**
 * Result data class for CTMC SSG reachability analysis
 */
data class CtmcSsgReachabilityResult(val stateSpace: Matrix,
                                     val stateSpaceAggr: Matrix,
                                     val stateSpaceHashed: Matrix,
                                     val nodeStateSpace: Map<StatefulNode, Matrix>,
                                     val sn: NetworkStruct)
/**
 * CTMC ssg reachability algorithms
 */
@Suppress("unused")
class CtmcSsgReachabilityAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}