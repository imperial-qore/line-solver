package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasLCFSPI {
    private SnHasLCFSPI() {}

    public static boolean snHasLCFSPI(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.LCFSPI);
    }
}
