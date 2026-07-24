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
 * A service section that processes jobs
 */
public class Server extends ServiceSection implements Serializable {
    /**
     * Creates a new server section for the specified job classes.
     * Initializes with a single server.
     * 
     * @param jobClasses the list of job classes this server will handle
     */
    public Server(List<JobClass> jobClasses) {
        super("Server");
        this.numberOfServers = 1;
    }

    /**
     * Initializes service processes for all job classes.
     * Sets exponential(0) as default service distribution with load-independent strategy.
     * 
     * @param jobClasses the list of job classes to initialize
     */
    private void initServers(List<JobClass> jobClasses) {
        for (JobClass jobClass : jobClasses) {
            this.serviceProcesses.put(jobClass, new ServiceBinding(jobClass, ServiceStrategy.LI, new Exp(0)));
        }
    }

}
