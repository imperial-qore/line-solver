package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

/** True when some station schedules by non-preemptive shortest job next. */
public final class SnHasSjn {
    private SnHasSjn() {}

    public static boolean snHasSjn(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.SJF);
    }
}
