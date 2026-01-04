/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;

import java.io.Serializable;
import java.util.List;

/**
 * A class switcher that depends on its local state
 */
public class StatefulClassSwitcher extends ClassSwitcher implements Serializable {
    public StatefulClassSwitcher(List<JobClass> jobClasses, String name) {
        super(jobClasses, name);
    }

    public StatefulClassSwitcher(List<JobClass> jobClasses) {
        this(jobClasses, "StatefulClassSwitcher");
    }
}
