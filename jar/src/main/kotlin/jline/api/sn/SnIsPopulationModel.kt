package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

/**
 * Checks if the model is a population model (only specific scheduling strategies without priorities or fork-join)
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return true if the model is a population model, false otherwise
 */
fun snIsPopulationModel(sn: NetworkStruct): Boolean {
    // Check if all scheduling strategies are allowed ones
    val allSchedulingValid = sn.sched.values.all { schedStrategy ->
        schedStrategy == SchedStrategy.INF || 
        schedStrategy == SchedStrategy.PS ||
        schedStrategy == SchedStrategy.PSPRIO ||
        schedStrategy == SchedStrategy.DPS ||
        schedStrategy == SchedStrategy.GPS ||
        schedStrategy == SchedStrategy.GPSPRIO ||
        schedStrategy == SchedStrategy.DPSPRIO ||
        schedStrategy == SchedStrategy.EXT
    }
    
    // Must have valid scheduling AND no priorities AND no fork-join
    return allSchedulingValid && !snHasPriorities(sn) && !snHasForkJoin(sn)
}
/**
 * Stochastic network IsPopulationModel algorithms
 */
@Suppress("unused")
class SnispopulationmodelAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}