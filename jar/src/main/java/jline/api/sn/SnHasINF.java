package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasINF {
    private SnHasINF() {}

    public static boolean snHasINF(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.INF);
    }
}
