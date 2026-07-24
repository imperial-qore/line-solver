package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasLCFS {
    private SnHasLCFS() {}

    public static boolean snHasLCFS(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.LCFS);
    }
}
