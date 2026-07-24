package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasLEPT {
    private SnHasLEPT() {}

    public static boolean snHasLEPT(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.LEPT);
    }
}
