/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.constant.RoutingStrategy;

import java.util.List;

/**
 * Output section that models the process of firing for a transition in a Stochastic Petri net model
 */
public class Firing extends OutputSection {
    protected List<JobClass> jobClasses;

    public Firing(List<JobClass> customerClasses) {
        super("Firing");
        this.jobClasses = customerClasses;
        this.initDispatcherJobClasses(customerClasses);
    }

    private void initDispatcherJobClasses(List<JobClass> customerClasses) {
        for (JobClass jobClass : customerClasses) {
            this.outputStrategies.add(new OutputStrategy(jobClass, RoutingStrategy.RAND));
        }
    }

}
