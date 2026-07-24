package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasSIRO {
    private SnHasSIRO() {}

    public static boolean snHasSIRO(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.SIRO);
    }
}
