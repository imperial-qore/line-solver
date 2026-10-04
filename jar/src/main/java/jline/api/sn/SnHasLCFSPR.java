package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasLCFSPR {
    private SnHasLCFSPR() {}

    public static boolean snHasLCFSPR(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.LCFSPR);
    }
}
