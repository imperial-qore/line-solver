/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.util.SerializableFunction;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;

/**
 * A job class switcher based on a static probability table
 */
public class ClassSwitcher extends ServiceSection implements Serializable {
    protected SerializableFunction<CSFunInput, Double> csFun;
    protected List<JobClass> jobClasses;

    public ClassSwitcher(List<JobClass> jobClasses, String name) {
        super(name);

        this.jobClasses = jobClasses;
        this.numberOfServers = 1;
        this.serviceProcesses = new HashMap<JobClass, ServiceBinding>();
    }

    public double applyCsFun(int r, int s) {
        return this.csFun.apply(new CSFunInput(r, s, null, null));
    }

}
