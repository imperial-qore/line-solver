/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for specifying a scheduling strategy type at stations
 */
public enum SchedStrategyType {
    PR, // preemptive resume
    PNR, // preemptive non-resume
    NP, // non-preemptive
    NPPrio // non-preemptive priority
}
