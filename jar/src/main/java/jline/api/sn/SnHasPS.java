package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasPS {
    private SnHasPS() {}

    public static boolean snHasPS(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.PS);
    }
}
