/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Exp;

import java.io.Serializable;
import java.util.List;

/**
 * A server shared by multiple jobs simultaneously
 */
public class SharedServer extends ServiceSection implements Serializable {

    /**
     * Creates a new shared server section for the specified customer classes.
     * A shared server allows multiple jobs to be processed simultaneously
     * by a single server resource.
     * 
     * @param customerClasses the list of customer classes this shared server will handle
     */
    public SharedServer(List<JobClass> customerClasses) {
        super("SharedServer");
        this.numberOfServers = 1;
    }

    /**
     * Initializes service processes for all customer classes.
     * Sets exponential(0) as default service distribution with load-independent strategy.
     * 
     * @param customerClasses the list of customer classes to initialize
     */
    private void initServers(List<JobClass> customerClasses) {
        for (JobClass jobClass : customerClasses) {
            this.serviceProcesses.put(jobClass, new ServiceBinding(jobClass, ServiceStrategy.LI, new Exp(0)));
        }
    }
}