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
 * Output section that routes jobs to nodes
 */
public class Dispatcher extends OutputSection implements Serializable {
    public Dispatcher(List<JobClass> customerClasses) {
        super("Dispatcher");
        initDispatcherJobClasses(customerClasses);
    }

    public void initDispatcherJobClasses(List<JobClass> customerClasses) {
        // Use setOutputStrategy to properly handle duplicate entries
        // The old approach used list index position which could create duplicates
        for (int r = 0; r < customerClasses.size(); r++) {
            JobClass jc = customerClasses.get(r);
            // Use the 2-param setOutputStrategy which properly removes all existing entries
            // for this class before adding a new one
            setOutputStrategy(jc, RoutingStrategy.DISABLED);
        }
    }
}
