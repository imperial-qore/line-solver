/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

/**
 * Alias for the Host class, i.e., a processor that can run Tasks
 */
public class Processor extends Host {
    public Processor(LayeredNetwork myLN, String name, int multiplicity, SchedStrategy scheduling, double quantum, double speedFactor) {
        super(myLN, name, multiplicity, scheduling, quantum, speedFactor);
    }

    public Processor(LayeredNetwork myLN, String name, int multiplicity, SchedStrategy scheduling, double quantum) {
        super(myLN, name, multiplicity, scheduling, quantum, 1);
    }

    public Processor(LayeredNetwork myLN, String name, int multiplicity, SchedStrategy scheduling) {
        super(myLN, name, multiplicity, scheduling, 0.001, 1);
    }

    public Processor(LayeredNetwork myLN, String name, int multiplicity) {
        super(myLN, name, multiplicity, SchedStrategy.PS, 0.001, 1);
    }

    public Processor(LayeredNetwork myLN, String name) {
        super(myLN, name, 1, SchedStrategy.PS, 0.001, 1);
    }
}
