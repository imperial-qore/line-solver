/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import java.io.Serializable;

/**
 * A service section of a Transition in a stochastic Petri net model
 */
public class Timing extends ServiceSection implements Serializable {
    public Timing(String name) {
        super(name);

        this.numberOfServers = 1;
    }

    public Timing() {
        this("Timing");
    }
}
