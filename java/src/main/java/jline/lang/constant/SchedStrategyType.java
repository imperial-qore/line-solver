package jline.lang.constant;

import java.util.List;

/**
 *  Constants for specifying a scheduling strategy type at stations
 */
public enum SchedStrategyType {
    PR, // preemptive resume
    PNR, // preemptive non-resume
    NP, // non-preemptive
    NPPrio // non-preemptive priority
}
