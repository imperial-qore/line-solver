package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasDPS {
    private SnHasDPS() {}

    public static boolean snHasDPS(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.DPS);
    }
}
