package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasSJF {
    private SnHasSJF() {}

    public static boolean snHasSJF(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.SJF);
    }
}
