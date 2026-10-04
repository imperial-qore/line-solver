/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying timing strategies at Petri net transitions
 */
public enum TimingStrategy implements Serializable {
    TIMED,
    IMMEDIATE
}