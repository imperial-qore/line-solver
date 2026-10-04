package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnIsPopulationModel {
    private SnIsPopulationModel() {}

    /**
     * Checks if the model is a population model (only specific scheduling strategies without priorities or fork-join)
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if the model is a population model, false otherwise
     */
    public static boolean snIsPopulationModel(NetworkStruct sn) {
        // Check if all scheduling strategies are allowed ones
        boolean allSchedulingValid = true;
        for (SchedStrategy schedStrategy : sn.sched.values()) {
            if (!(schedStrategy == SchedStrategy.INF
                    || schedStrategy == SchedStrategy.PS
                    || schedStrategy == SchedStrategy.PSPRIO
                    || schedStrategy == SchedStrategy.DPS
                    || schedStrategy == SchedStrategy.GPS
                    || schedStrategy == SchedStrategy.GPSPRIO
                    || schedStrategy == SchedStrategy.DPSPRIO
                    || schedStrategy == SchedStrategy.EXT)) {
                allSchedulingValid = false;
                break;
            }
        }

        // Must have valid scheduling AND no priorities AND no fork-join
        return allSchedulingValid && !SnHasPriorities.snHasPriorities(sn) && !SnHasForkJoin.snHasForkJoin(sn);
    }
}
