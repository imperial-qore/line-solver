package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasDPSPRIO {
    private SnHasDPSPRIO() {}

    public static boolean snHasDPSPRIO(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.DPSPRIO);
    }
}
