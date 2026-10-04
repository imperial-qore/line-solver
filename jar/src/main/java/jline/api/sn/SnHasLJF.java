package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasLJF {
    private SnHasLJF() {}

    public static boolean snHasLJF(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.LJF);
    }
}
