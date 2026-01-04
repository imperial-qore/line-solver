/**
 * Stochastic Network State Validation Utility
 * 
 * Validates the consistency and feasibility of queueing network states by checking
 * job conservation laws, capacity constraints, and scheduling discipline requirements.
 * Essential for ensuring model correctness before numerical analysis.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.io.line_warning
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.util.matrix.Matrix

/**
 * Checks if the network state is valid.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return true if the state is valid, false otherwise
 */

fun snIsStateValid(sn: NetworkStruct): Boolean {
    val snTmp: NetworkStruct = sn.copy()
    val nir = Matrix(snTmp.nstations, snTmp.nclasses)
    val sir = Matrix(snTmp.nstations, snTmp.nclasses)

    for (ist in 0..<snTmp.nstations) {
        val statefulIdx = snTmp.stationToStateful.get(ist).toInt()
        if (snTmp.state[snTmp.stations[ist]]!!.numRows > 1) {
            if (snTmp.stateprior[snTmp.stations[ist]]!!.elementMax() < 1 - GlobalConstants.FineTol) {
                line_warning(mfilename(object : Any() {}),
                    "isStateValid will ignore some states of station %d, define a unique initial state to address this problem.\n",
                    ist)
            }
            val initialState = Matrix(1, snTmp.state[snTmp.stations[ist]]!!.numCols)
            Matrix.extractRows(snTmp.state[snTmp.stations[ist]], 0, 1, initialState)
            snTmp.state[snTmp.stations[ist]] = initialState
        }

        val stats = toMarginal(snTmp,
            snTmp.stationToNode[0, ist].toInt(),
            snTmp.state[snTmp.stations[ist]],
            null,
            null,
            null,
            null,
            null)

        for (i in 0..<snTmp.nclasses) {
            nir[ist, i] = stats.nir[0, i]
            sir[ist, i] = stats.sir[0, i]
        }
    }

    return State.isValid(snTmp, nir, sir)
}
/**
 * Stochastic network IsStateValid algorithms
 */
@Suppress("unused")
class SnisstatevalidAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}