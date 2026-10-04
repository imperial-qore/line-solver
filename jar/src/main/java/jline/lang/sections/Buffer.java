/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.constant.SchedStrategyType;

import java.io.Serializable;
import java.util.List;

/**
 * Input buffer of a station
 */
public class Buffer extends InputSection implements Serializable {
    protected int size;

    public Buffer(List<JobClass> classes) {
        super("Buffer");

        this.size = -1;
        this.schedPolicy = SchedStrategyType.NP;
    }

    public void setSize(int size) {
        this.size = size;
    }
}
