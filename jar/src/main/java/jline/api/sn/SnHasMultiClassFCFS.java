package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;

import java.util.Map;

public final class SnHasMultiClassFCFS {
    private SnHasMultiClassFCFS() {}

    public static boolean snHasMultiClassFCFS(NetworkStruct sn) {
        for (Map.Entry<Station, SchedStrategy> entry : sn.sched.entrySet()) {
            Station key = entry.getKey();
            SchedStrategy value = entry.getValue();
            if (value == SchedStrategy.FCFS) {
                int nnz_rates = 0;
                for (int j = 0; j < sn.nclasses; j++) {
                    if (sn.rates.get(key.getStationIdx(), j) > 0) {
                        nnz_rates++;
                    }
                }
                if (nnz_rates > 1) {
                    return true;
                }
            }
        }
        return false;
    }
}
