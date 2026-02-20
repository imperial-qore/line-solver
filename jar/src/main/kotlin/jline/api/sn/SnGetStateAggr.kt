package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.nodes.StatefulNode
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginalAggr
import jline.util.matrix.Matrix


/**
 * Aggregates the state of the network.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return a map of stateful nodes to their aggregated state matrices
 */
fun snGetStateAggr(sn: NetworkStruct): Map<StatefulNode, Matrix> {
    val initialState = sn.state
    val initialStateAggr: MutableMap<StatefulNode, Matrix> = HashMap()
    for (isf in 0..<initialState.size) {
        val ind = sn.statefulToNode[isf].toInt()
        initialState[sn.stateful[isf]]
        val aggrState = toMarginalAggr(sn, ind, initialState[sn.stateful[isf]], null, null, null, null, null)
        initialStateAggr[sn.stateful[isf]] = aggrState.nir
    }
    return initialStateAggr
}
/**
 * Stochastic network GetStateAggr algorithms
 */
@Suppress("unused")
class SngetstateaggrAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}