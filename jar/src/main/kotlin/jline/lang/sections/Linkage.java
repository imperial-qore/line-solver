/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.constant.RoutingStrategy;

import java.io.Serializable;
import java.util.List;

/**
 * Output section of a Place in a Stochastic Petri net model
 */
public class Linkage extends OutputSection implements Serializable {
    public Linkage(List<JobClass> customerClasses) {
        super("Linkage");
        this.initDispatcherJobClasses(customerClasses);
    }

    public void initDispatcherJobClasses(List<JobClass> customerClasses) {
        // Use setOutputStrategy to properly handle duplicate entries
        // The old approach used list index position which could create duplicates
        for (int r = 0; r < customerClasses.size(); r++) {
            JobClass jc = customerClasses.get(r);
            // Use the 2-param setOutputStrategy which properly removes all existing entries
            // for this class before adding a new one
            setOutputStrategy(jc, RoutingStrategy.RAND);
        }
    }
}
