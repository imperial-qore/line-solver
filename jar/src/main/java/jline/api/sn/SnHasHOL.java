package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasHOL {
    private SnHasHOL() {}

    public static boolean snHasHOL(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.HOL) || sn.sched.containsValue(SchedStrategy.FCFSPRIO);
    }
}
