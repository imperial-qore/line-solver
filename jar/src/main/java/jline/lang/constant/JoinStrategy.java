/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying a join strategy
 */
public enum JoinStrategy implements Serializable {
    STD,
    PARTIAL,
    Quorum,
    Guard
}
