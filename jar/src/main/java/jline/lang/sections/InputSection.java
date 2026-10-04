/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.ServiceBinding;
import jline.lang.constant.SchedStrategyType;

import java.io.Serializable;

/**
 * Input section of a station
 */
public class InputSection extends Section implements Serializable {
    protected SchedStrategyType schedPolicy;

    public InputSection(String className) {
        super(className);
    }

    public void setServiceProcess(ServiceBinding serviceProcess) {

    }
}
