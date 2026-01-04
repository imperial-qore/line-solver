/**
 * @file Stochastic network model type classifier for closed models
 * 
 * Identifies closed queueing network models with finite job populations and no external
 * arrivals. Closed models are fundamental in capacity planning and system design,
 * enabling analysis of systems with constrained resources and finite user populations.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network model is closed.
 *
 * @param sn the NetworkStruct object representing the network
 * @return true if the network has finite jobs and no infinite jobs, false otherwise
 */
fun snIsClosedModel(sn: NetworkStruct): Boolean {
    return !sn.njobs.hasInfinite() && sn.njobs.hasFinite()
}
/**
 * Stochastic network IsClosedModel algorithms
 */
@Suppress("unused")
class SnisclosedmodelAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}