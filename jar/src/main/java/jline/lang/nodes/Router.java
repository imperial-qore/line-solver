/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.ServiceBinding;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Distribution;
import jline.lang.sections.Buffer;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.ServiceTunnel;

import java.io.Serializable;
import java.util.List;

/**
 * A node that routes jobs without imposing any delay
 */
public class Router extends StatefulNode implements Serializable {
    protected int cap;
    protected int numberOfServers;
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;

    /**
     * Creates a new router node that routes jobs without delay.
     * Uses FCFS scheduling and a single server by default.
     * 
     * @param model the network model to add this router to
     * @param name the name for this router node
     */
    public Router(Network model, String name) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        this.server = new ServiceTunnel();
        this.input = new Buffer(jobClasses);
        this.output = new Dispatcher(jobClasses);
        this.cap = Integer.MAX_VALUE;

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.numberOfServers = 1;
        this.setModel(model);
        this.model.addNode(this);
    }

    /**
     * Gets the scheduling strategy used by this router.
     * 
     * @return the scheduling strategy
     */
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    /**
     * Gets the service distribution for a specific job class.
     * 
     * @param jobClass the job class to query
     * @return the service distribution for the job class
     */
    public Distribution getServiceProcess(JobClass jobClass) {
        return this.server.getServiceDistribution(jobClass);
    }

    /**
     * Sets the scheduling policy type for this router.
     * 
     * @param schedPolicy the scheduling policy type (PR for preemptive, NP for non-preemptive)
     */
    public void setSchedPolicy(SchedStrategyType schedPolicy) {
        this.schedPolicy = schedPolicy;
    }

    /**
     * Sets the service distribution for a specific job class.
     * 
     * @param jobClass the job class to configure
     * @param distribution the service distribution to set
     */
    public void setService(JobClass jobClass, Distribution distribution) {
        this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, distribution));
    }

}
