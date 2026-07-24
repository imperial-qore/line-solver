package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasPolling {
    private SnHasPolling() {}

    public static boolean snHasPolling(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.POLLING);
    }
}
