/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.sections.Linkage;
import jline.lang.sections.ServiceTunnel;
import jline.lang.sections.Storage;

import java.io.Serializable;
import java.util.HashMap;

/**
 * A queueing station within a Network model
 */
public class Place extends Station implements Serializable {
    protected HashMap<JobClass, SchedStrategy> schedStrategies;

    public Place(Network model, String name) {
        this(model, name, SchedStrategy.FCFS);
    }

    public Place(Network model, String name, SchedStrategy schedStrategy) {
        super(name);
        this.schedStrategies = new HashMap<JobClass, SchedStrategy>();
        this.schedStrategyPar = new HashMap<JobClass, Double>();

        this.setModel(model);
        this.model.addNode(this);
        this.schedStrategy = schedStrategy;
        this.input = new Storage(model.getClasses());
        this.server = new ServiceTunnel();
        this.output = new Linkage(model.getClasses());
        this.schedStrategy = schedStrategy;
        this.numberOfServers = 1;
    }

    public void init() {
        for (JobClass jobclass : this.model.getClasses()) {
            this.classCap.put(jobclass, Integer.MAX_VALUE);
            this.cap = Integer.MAX_VALUE;
            this.setSchedStrategy(jobclass, SchedStrategy.FCFS);
            this.setDropRule(jobclass, DropStrategy.WaitingQueue);
        }
    }

    public void setClassCapacity(JobClass jobclass, int capacity) {
        this.classCap.put(jobclass, capacity);
    }

    public void setSchedStrategy(JobClass jobClass, SchedStrategy strategy) {
        this.schedStrategies.put(jobClass, strategy);
    }

}
