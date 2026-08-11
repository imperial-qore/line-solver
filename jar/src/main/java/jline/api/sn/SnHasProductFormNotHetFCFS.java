package jline.api.sn;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasProductFormNotHetFCFS {
    private SnHasProductFormNotHetFCFS() {}

    /**
     * Checks if the network satisfies product-form assumptions (does not have heterogeneous FCFS)
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return true if the network satisfies the product-form assumptions
     */
    public static boolean snHasProductFormNotHetFCFS(NetworkStruct sn) {
        boolean ret = true;
        for (int i = 0; i < sn.sched.size(); i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            ret = ret && (s == SchedStrategy.INF
                    || s == SchedStrategy.PS
                    || s == SchedStrategy.FCFS
                    || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.EXT);
        }
        ret = ret && !SnHasPriorities.snHasPriorities(sn);
        ret = ret && !SnHasForkJoin.snHasForkJoin(sn);
        ret = ret && !SnHasSDRouting.snHasSDRouting(sn);

        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                for (int j = 0; j < sn.scv.getNumCols(); j++) {
                    double scv = sn.scv.get(i, j);
                    if (Double.isFinite(scv) && scv > 0) {
                        ret = ret && (scv > 1 - GlobalConstants.FineTol) && (scv < 1 + GlobalConstants.FineTol);
                    }
                }
            }
        }

        return ret;
    }
}
