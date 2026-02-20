/**
 * @file Stochastic network model type classifier for open models
 * 
 * Identifies open queueing network models with external arrivals and infinite 
 * job populations. Open models are essential for analyzing systems with external 
 * traffic sources and unlimited capacity for job creation.
 * 
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network is an open model.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return true if the network has infinite jobs and no finite jobs, false otherwise
 */

fun snIsOpenModel(sn: NetworkStruct): Boolean {
    return !sn.njobs.hasFinite() && sn.njobs.hasInfinite()
}
/**
 * Stochastic network IsOpenModel algorithms
 */
@Suppress("unused")
class SnisopenmodelAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}