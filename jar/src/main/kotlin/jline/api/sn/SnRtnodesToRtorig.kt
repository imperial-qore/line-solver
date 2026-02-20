package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix
import jline.api.mc.dtmc_stochcomp

/**
 * Converts routing matrices from nodes to original format, specifically handling class switching nodes
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return a pair containing:
 *         - rtorigcell: a 2D list of matrices indexed by [class_from][class_to] containing routing probabilities
 *         - rtorig: the full routing matrix after stochastic complement
 */
fun snRtnodesToRtorig(sn: NetworkStruct): Pair<List<List<Matrix>>, Matrix> {
    val K = sn.nclasses
    val rtnodes = sn.rtnodes
    
    // Find CS shift point - the index before the first class switching node
    var csshift = sn.nnodes
    for (ind in 0..<sn.nnodes) {
        if (sn.nodenames[ind].startsWith("CS_")) {
            csshift = ind
            break
        }
    }
    
    // Build columns to keep for stochastic complement
    val colToKeep = mutableListOf<Int?>()
    for (ind in 0..<csshift) {
        for (k in 0..<K) {
            colToKeep.add(ind * K + k)
        }
    }
    
    // Compute stochastic complement to remove class switching nodes
    val rtorig = dtmc_stochcomp(rtnodes, colToKeep)
    
    // Initialize the cell array equivalent - K x K matrices, each of size csshift x csshift
    val rtorigcell = List(K) { List(K) { Matrix(csshift, csshift) } }
    
    // Replace NaNs with 0 for cache routing probabilities
    for (i in 0..<rtorig.numRows) {
        for (j in 0..<rtorig.numCols) {
            if (rtorig[i, j].isNaN()) {
                rtorig[i, j] = 0.0
            }
        }
    }
    
    // Populate cell array with routing probabilities between classes
    for (ind in 0..<csshift) {
        if (sn.nodetype[ind] != NodeType.Sink) {
            for (jnd in 0..<csshift) {
                for (r in 0..<K) {
                    for (s in 0..<K) {
                        rtorigcell[r][s][ind, jnd] = rtorig[(ind * K) + r, (jnd * K) + s]
                    }
                }
            }
        }
    }
    
    return Pair(rtorigcell, rtorig)
}
/**
 * Stochastic network RtnodesToRtorig algorithms
 */
@Suppress("unused")
class SnrtnodestortorigAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}