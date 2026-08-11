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
 * A preemptive service section that can interrupt lower priority jobs
 * to serve higher priority jobs.
 * <p>
 * This server implements preemptive scheduling where an arriving high-priority
 * job can interrupt the service of a lower-priority job currently being served.
 * The interrupted job may resume service later depending on the preemption policy.
 */
public class PreemptiveServer extends ServiceSection implements Serializable {
    /**
     * Creates a new preemptive server section for the specified customer classes.
     * Initializes with a single server that can preempt lower priority jobs.
     * 
     * @param customerClasses the list of customer classes with different priorities
     */
    public PreemptiveServer(List<JobClass> customerClasses) {
        super("PreemptiveServer");
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
