package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasGPS {
    private SnHasGPS() {}

    public static boolean snHasGPS(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.GPS);
    }
}
