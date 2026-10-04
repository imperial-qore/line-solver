package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasSrpt {
    private SnHasSrpt() {}

    public static boolean sn_has_srpt(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.SRPT);
    }
}
