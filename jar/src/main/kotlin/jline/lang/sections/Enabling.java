/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.constant.SchedStrategyType;

import java.io.Serializable;

/**
 * A section that models enabling conditions in a stochastic Petri net transition
 */
public class Enabling extends InputSection implements Serializable {
    protected int size;

    public Enabling(String name) {
        super(name);

        this.size = -1;
        this.schedPolicy = SchedStrategyType.NP;
    }

    public Enabling() {
        this("Enabling");
    }
}
