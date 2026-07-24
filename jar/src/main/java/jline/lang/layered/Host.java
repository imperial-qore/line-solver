/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

import java.util.ArrayList;
import java.util.List;

/**
 * A processor that can run Tasks
 */
public class Host extends LayeredNetworkElement {
    protected int multiplicity;
    protected int replication;
    protected SchedStrategy scheduling;
    protected double quantum;
    protected double speedFactor;
    protected List<Task> tasks;
    private int ID;

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, double quantum, double speedFactor) {
        super(name);
        this.multiplicity = multiplicity;
        this.replication = 1;
        this.scheduling = scheduling;
        this.quantum = quantum;
        this.speedFactor = speedFactor;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(), this);
        model.nodes.put(model.hosts.size(), this);
    }

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, double quantum) {
        this(model, name, multiplicity, scheduling, quantum, 1);
    }

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        this(model, name, multiplicity, scheduling, 0.01, 1);
    }

    public Host(LayeredNetwork model, String name, int multiplicity) {
        this(model, name, multiplicity, SchedStrategy.PS, 0.01, 1);
    }

    public Host(LayeredNetwork model, String name) {
        this(model, name, 1, SchedStrategy.PS, 0.01, 1);
    }

    public void addTask(Task newTask) {
        tasks.add(newTask);
    }

    public boolean removeTask(Task newTask) {
        return tasks.remove(newTask);
    }

    public void setReplication(int replication) {
        this.replication = replication;
    }

    /**
     * Returns the list of tasks running on this host.
     *
     * @return the list of tasks
     */
    public List<Task> getTasks() {
        return tasks;
    }

    public int getMultiplicity() {
        return multiplicity;
    }

    public void setMultiplicity(int multiplicity) {
        this.multiplicity = multiplicity;
    }

    public int getReplication() {
        return replication;
    }

    public SchedStrategy getScheduling() {
        return scheduling;
    }

    public void setScheduling(SchedStrategy scheduling) {
        this.scheduling = scheduling;
    }

    public double getQuantum() {
        return quantum;
    }

    public void setQuantum(double quantum) {
        this.quantum = quantum;
    }

    public double getSpeedFactor() {
        return speedFactor;
    }

    public void setSpeedFactor(double speedFactor) {
        this.speedFactor = speedFactor;
    }
}
