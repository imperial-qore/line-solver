/**
 * @file Stochastic network state-dependent routing checker
 *
 * Determines if a queueing network has state-dependent routing strategies that
 * violate product-form assumptions. State-dependent routing includes Round-Robin,
 * Weighted Round-Robin, Join Shortest Queue, Power of K Choices, and Reinforcement Learning.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.RoutingStrategy

/**
 * Checks if the network has state-dependent routing strategies
 *
 * Product-form requires state-independent (Markovian) routing.
 * PROB and RAND are product-form compatible.
 * RROBIN, WRROBIN, JSQ, KCHOICES, RL are state-dependent and violate product-form.
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return true if the network has state-dependent routing, false otherwise
 */
fun snHasSDRouting(sn: NetworkStruct): Boolean {
    if (sn.routing == null || sn.routing.isEmpty()) {
        return false
    }

    for (nodeEntry in sn.routing.entries) {
        val classMap = nodeEntry.value
        if (classMap != null) {
            for (classEntry in classMap.entries) {
                val strategy = classEntry.value
                if (strategy == RoutingStrategy.RROBIN ||
                    strategy == RoutingStrategy.WRROBIN ||
                    strategy == RoutingStrategy.JSQ ||
                    strategy == RoutingStrategy.KCHOICES ||
                    strategy == RoutingStrategy.RL) {
                    return true
                }
            }
        }
    }
    return false
}

/**
 * Stochastic network HasSDRouting algorithms
 */
@Suppress("unused")
class SnHasSDRoutingAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
