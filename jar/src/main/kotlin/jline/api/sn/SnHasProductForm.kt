/**
 * @file Stochastic network product-form solution checker
 * 
 * Determines if a queueing network has a known product-form solution by validating
 * scheduling disciplines, class structures, and network topology. Product-form networks
 * enable efficient exact analysis using algorithms like MVA and convolution.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

/**
 * Checks if the network has a known product-form solution
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasProductForm(sn: NetworkStruct): Boolean {
    var ret = true
    var hasLcfs = false
    var hasLcfspr = false
    var lcfsCount = 0
    var lcfsprCount = 0

    for (i in 0..<sn.sched.size) {
        val sched = sn.sched[sn.stations[i]]
        when (sched) {
            SchedStrategy.LCFS -> {
                hasLcfs = true
                lcfsCount++
            }
            SchedStrategy.LCFSPR -> {
                hasLcfspr = true
                lcfsprCount++
            }
            SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.EXT -> {}
            else -> ret = false
        }
    }

    // Special case: LCFS + LCFSPR 2-station network has product-form solution
    // Reference: G. Casale, "A family of multiclass LCFS queueing networks with
    // order-dependent product-form solutions", QUESTA 2026.
    if (hasLcfs && hasLcfspr && lcfsCount == 1 && lcfsprCount == 1) {
        // This is a valid product-form configuration
    } else if (hasLcfs) {
        // LCFS alone (without LCFSPR pairing) is not product-form
        ret = false
    }

    ret = ret && !snHasMultiClassHeterFCFS(sn)
    ret = ret && !snHasPriorities(sn)
    ret = ret && !snHasForkJoin(sn)
    ret = ret && !snHasSDRouting(sn)
    return ret
}
/**
 * Stochastic network HasProductForm algorithms
 */
@Suppress("unused")
class SnhasproductformAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}