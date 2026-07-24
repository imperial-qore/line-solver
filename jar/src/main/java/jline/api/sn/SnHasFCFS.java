/**
 * @file Stochastic network FCFS scheduling detector
 *
 * Identifies queueing networks using First-Come-First-Served scheduling disciplines.
 * FCFS detection is crucial for algorithm selection and performance analysis since
 * FCFS networks often admit efficient solution techniques.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

public final class SnHasFCFS {
    private SnHasFCFS() {}

    public static boolean snHasFCFS(NetworkStruct sn) {
        return sn.sched.containsValue(SchedStrategy.FCFS);
    }
}
