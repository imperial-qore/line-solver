/**
 * @file Load Sharing Network Maximum Multiplicity Computation
 * 
 * Computes maximum multiplicity constraints for load sharing network (LSN)
 * analysis. Essential for determining feasible population bounds and
 * capacity constraints in layered queueing network models.
 * 
 * @since LINE 3.0
 */
package jline.api.lsn

import jline.GlobalConstants.Inf

import jline.lang.layered.LayeredNetworkElement
import jline.lang.layered.LayeredNetworkStruct
import jline.util.graph.DirectedGraph
import jline.util.matrix.Matrix
import kotlin.math.min

/**
 * Computes the maximum multiplicity (number of concurrent instances) that can be 
 * sustained by each task in a layered network while respecting capacity constraints.
 * 
 * This function uses flow analysis based on Kahn's topological sorting algorithm
 * to determine the maximum sustainable throughput for each task, considering
 * both the incoming flow and the multiplicity constraints.
 *
 * @param lsn The layered network structure containing task dependencies and constraints
 * @return Matrix of maximum multiplicities for each task in the network
 */
fun lsnMaxMultiplicity(lsn: LayeredNetworkStruct): Matrix {
    val ag = Matrix(lsn.dag)
    val n = lsn.dag.numRows
    for (ist in 0..<n) {
        for (jst in 0..<n) {
            if (lsn.dag[ist, jst] > 0) {
                ag[ist, jst] = 1.0
            }
        }
    }
    // Make a copy of mult to avoid modifying lsn.mult directly
    val mult = Matrix(lsn.mult)
    val type = lsn.type
    val isref = lsn.isref

    val order = DirectedGraph.kahn(ag)

    val inflow = Matrix(n, 1)
    for (ist in 0..<n) {
        if (type[ist] == LayeredNetworkElement.TASK.toDouble() && isref[ist] != 0.0) {
            inflow[ist, 0] = mult[ist]
        }
    }

    val outflow = Matrix(n, 1)

    if (mult.length() < n) {
        val cur_end = mult.length()
        mult.expandMatrix(1, n, n)
        for (ist in cur_end..<n) {
            mult[ist] = Inf
        }
    }

    for (k in 0..<n) {
        val ist = order[k].toInt()
        val inflowVal = inflow[ist]
        val multVal = mult[ist]
        outflow[ist, 0] = min(inflowVal, multVal)

        for (jst in 0..<n) {
            if (jst != ist && ag[ist, jst] != 0.0) {
                inflow[jst] = inflow[jst] + outflow[ist]
            }
        }
    }

    for (ist in 0..<n) {
        if (type[ist] == LayeredNetworkElement.TASK.toDouble() && mult[ist] == Inf && isref[ist] == 0.0) {
            outflow[ist] = Inf
        }
    }

    return outflow
}