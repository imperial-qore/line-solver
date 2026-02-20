/**
 * @file Stochastic network multi-class detection utility
 * 
 * Identifies queueing networks with multiple job classes, which require specialized
 * analysis algorithms that account for class-dependent service parameters, routing
 * probabilities, and scheduling policies in multi-class queueing systems.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct

fun snHasMultiClass(sn: NetworkStruct): Boolean {
    return sn.nclasses > 1
}
/**
 * Stochastic network HasMultiClass algorithms
 */
@Suppress("unused")
class SnhasmulticlassAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}