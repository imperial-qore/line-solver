package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasSEPT {
    private SnHasSEPT() {}

    public static boolean snHasSEPT(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.SEPT);
    }
}
