package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasGPSPRIO {
    private SnHasGPSPRIO() {}

    public static boolean snHasGPSPRIO(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.GPSPRIO);
    }
}
