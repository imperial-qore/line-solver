package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasPSPRIO {
    private SnHasPSPRIO() {}

    public static boolean snHasPSPRIO(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.PSPRIO);
    }
}
