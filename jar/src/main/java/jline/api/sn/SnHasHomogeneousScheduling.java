package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasHomogeneousScheduling {
    private SnHasHomogeneousScheduling() {}

    /**
     * Checks if the network uses an identical scheduling strategy at every station
     *
     * @param sn       - NetworkStruct object for the queueing network model
     * @param strategy - Scheduling strategy
     * @return boolean
     */
    public static boolean snHasHomogeneousScheduling(NetworkStruct sn, SchedStrategy strategy) {
        int stratCount = 0;
        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == strategy) stratCount++;
        }
        return stratCount == sn.sched.size();
    }
}
