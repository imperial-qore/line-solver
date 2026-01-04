/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for specifying service strategies at stations
 */
public enum ServiceStrategy {
    LI, // load-independent
    LD, // load-dependent
    CD, // class-dependent
    SD, // state-dependent
    JD  // joint-dependent
}
