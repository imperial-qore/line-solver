/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.ServiceBinding;
import jline.lang.constant.DepartureDiscipline;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;
import jline.lang.sections.InfiniteServer;
import jline.lang.sections.Linkage;
import jline.lang.sections.Server;
import jline.lang.sections.ServiceTunnel;
import jline.lang.sections.SharedServer;
import jline.lang.sections.Storage;

import java.io.Serializable;
import java.util.HashMap;

/**
 * A place within a stochastic Petri net / queueing Petri net model.
 *
 * <p>By default a Place is an <em>ordinary place</em> (tokens are immediately available to the
 * output transitions, modelled by a {@link ServiceTunnel} pass-through server). When a service
 * process is assigned via {@link #setService(JobClass, Distribution)}, the place becomes a
 * <em>queueing place</em> (QPN semantics): incoming tokens are served in an embedded queue under
 * the place's scheduling strategy and, on completion, move to a depository from which they become
 * available to the output transitions according to the departure discipline.</p>
 */
public class Place extends Station implements Serializable {
    protected HashMap<JobClass, SchedStrategy> schedStrategies;

    /** True once a service process has been assigned, i.e. this is a queueing place. */
    protected boolean queueing;

    /** Per-class service (queue) processes of the embedded queue; empty for ordinary places. */
    protected HashMap<JobClass, Distribution> serviceProcesses;

    /** Per-class departure discipline of the depository. */
    protected HashMap<JobClass, DepartureDiscipline> departureDiscipline;

    public Place(Network model, String name) {
        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        this(model, name, SchedStrategy.INF);
    }

    public Place(Network model, String name, SchedStrategy schedStrategy) {
        super(name);
        this.schedStrategies = new HashMap<JobClass, SchedStrategy>();
        this.schedStrategyPar = new HashMap<JobClass, Double>();
        this.serviceProcesses = new HashMap<JobClass, Distribution>();
        this.departureDiscipline = new HashMap<JobClass, DepartureDiscipline>();
        this.queueing = false;

        this.setModel(model);
        this.model.addNode(this);
        this.schedStrategy = schedStrategy;
        this.input = new Storage(model.getClasses());
        this.server = new ServiceTunnel();
        this.output = new Linkage(model.getClasses());
        this.numberOfServers = 1;
    }

    public void init() {
        for (JobClass jobclass : this.model.getClasses()) {
            // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
            if (!this.classCap.containsKey(jobclass)) {
                this.classCap.put(jobclass, Integer.MAX_VALUE);
            }
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

    /**
     * Assigns a service process to a token color, turning this ordinary place into a queueing
     * place. The embedded queue serves tokens under the place's scheduling strategy; the first
     * such call installs the concrete server section for that strategy.
     *
     * @param jobClass     the token color (job class)
     * @param distribution the service-time distribution of the embedded queue for this color
     */
    public void setService(JobClass jobClass, Distribution distribution) {
        if (!this.queueing) {
            this.installQueueServer();
            this.queueing = true;
        }
        this.serviceProcesses.put(jobClass, distribution);
        if (!this.departureDiscipline.containsKey(jobClass)) {
            this.departureDiscipline.put(jobClass, DepartureDiscipline.Normal);
        }
        Distribution actual = distribution.isImmediate() ? Immediate.getInstance() : distribution;
        if (this.server != null) {
            this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, actual));
        }
    }

    /**
     * Installs the concrete server section matching the place's scheduling strategy, replacing the
     * default {@link ServiceTunnel} pass-through used by ordinary places.
     */
    protected void installQueueServer() {
        switch (this.schedStrategy) {
            case INF:
                this.server = new InfiniteServer(this.model.getClasses());
                this.numberOfServers = Integer.MAX_VALUE;
                break;
            case PS:
            case DPS:
            case GPS:
            case LPS:
                this.server = new SharedServer(this.model.getClasses());
                break;
            default:
                this.server = new Server(this.model.getClasses());
                break;
        }
    }

    /** @return true if this place is a queueing place (has an embedded queue). */
    public boolean isQueueing() {
        return this.queueing;
    }

    public Distribution getService(JobClass jobClass) {
        return this.serviceProcesses.get(jobClass);
    }

    public void setDepartureDiscipline(JobClass jobClass, DepartureDiscipline discipline) {
        this.departureDiscipline.put(jobClass, discipline);
    }

    public DepartureDiscipline getDepartureDiscipline(JobClass jobClass) {
        return this.departureDiscipline.getOrDefault(jobClass, DepartureDiscipline.Normal);
    }

    public void setNumberOfServers(int numberOfServers) {
        this.numberOfServers = numberOfServers;
    }

    /**
     * Alias for {@link StatefulNode#setState(int)} using Petri-net terminology:
     * sets the initial token marking of this place.
     *
     * @param marking initial number of tokens
     */
    public void setMarking(int marking) {
        this.setState(marking);
    }

    /**
     * Alias for {@link StatefulNode#setState(jline.util.matrix.Matrix)} using
     * Petri-net terminology: sets the initial per-class token marking of this place.
     *
     * @param marking initial token counts per class
     */
    public void setMarking(jline.util.matrix.Matrix marking) {
        this.setState(marking);
    }

}
