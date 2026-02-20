/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.DropStrategy;
import jline.lang.nodes.Node;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Collection of stations with constraints on the number of admitted jobs
 */
public class Region implements Serializable {

    private static final long serialVersionUID = 1L;

    public static final int UNBOUNDED = -1;

    private String name;
    public List<Node> nodes;
    public List<JobClass> classes;
    public Map<JobClass, Integer> classMaxJobs;
    public Map<JobClass, Integer> classMaxMemory;
    public Map<JobClass, DropStrategy> dropRule;
    public Map<JobClass, Integer> classSize;
    public Map<JobClass, Double> classWeight;
    private int globalMaxJobs;
    private int globalMaxMemory;

    public Region(List<Node> nodes, List<JobClass> classes) {
        this.nodes = nodes;
        this.classes = classes;
        this.globalMaxJobs = UNBOUNDED;
        this.globalMaxMemory = UNBOUNDED;

        // Initialize Maps for per-class properties
        this.classMaxJobs = new HashMap<>();
        this.classMaxMemory = new HashMap<>();
        this.dropRule = new HashMap<>();
        this.classSize = new HashMap<>();
        this.classWeight = new HashMap<>();

        // Initialize default values for each class
        for (JobClass jobClass : classes) {
            this.classMaxJobs.put(jobClass, UNBOUNDED);
            this.classMaxMemory.put(jobClass, UNBOUNDED);
            this.dropRule.put(jobClass, DropStrategy.WaitingQueue);
            this.classSize.put(jobClass, 1);
            this.classWeight.put(jobClass, 1.0);  // Default weight = 1.0
        }
    }

    // Global methods
    public int getGlobalMaxJobs() {
        return globalMaxJobs;
    }

    public void setGlobalMaxJobs(int njobs) {
        this.globalMaxJobs = njobs;
    }

    public int getGlobalMaxMemory() {
        return globalMaxMemory;
    }

    public void setGlobalMaxMemory(int memlim) {
        this.globalMaxMemory = memlim;
    }

    // Per-class methods
    public void setClassMaxJobs(JobClass jobClass, int njobs) {
        classMaxJobs.put(jobClass, njobs);
    }

    public int getClassMaxJobs(JobClass jobClass) {
        return classMaxJobs.getOrDefault(jobClass, UNBOUNDED);
    }

    public void setClassMaxMemory(JobClass jobClass, int memlim) {
        classMaxMemory.put(jobClass, memlim);
    }

    public int getClassMaxMemory(JobClass jobClass) {
        return classMaxMemory.getOrDefault(jobClass, UNBOUNDED);
    }

    /**
     * Sets the drop strategy for a specific job class.
     *
     * @param jobClass the job class to configure
     * @param strategy the drop strategy to apply (DROP, WaitingQueue, BlockingAfterService, etc.)
     */
    public void setDropRule(JobClass jobClass, DropStrategy strategy) {
        dropRule.put(jobClass, strategy);
    }

    /**
     * Sets the drop rule using a boolean for backwards compatibility.
     *
     * @param jobClass the job class to configure
     * @param isDropEnabled true for DROP strategy, false for WaitingQueue
     */
    public void setDropRule(JobClass jobClass, boolean isDropEnabled) {
        if (isDropEnabled) {
            dropRule.put(jobClass, DropStrategy.Drop);
        } else {
            dropRule.put(jobClass, DropStrategy.WaitingQueue);
        }
    }

    /**
     * Gets the drop strategy for a specific job class.
     *
     * @param jobClass the job class to query
     * @return the drop strategy for the job class
     */
    public DropStrategy getDropStrategy(JobClass jobClass) {
        return dropRule.getOrDefault(jobClass, DropStrategy.WaitingQueue);
    }

    /**
     * Gets the drop rule as a boolean for backwards compatibility.
     * Returns true if the strategy is DROP, false otherwise.
     *
     * @param jobClass the job class to query
     * @return true if drop is enabled, false otherwise
     */
    public boolean getDropRule(JobClass jobClass) {
        DropStrategy strategy = dropRule.getOrDefault(jobClass, DropStrategy.WaitingQueue);
        return strategy == DropStrategy.Drop;
    }

    public void setClassSize(JobClass jobClass, int size) {
        classSize.put(jobClass, size);
    }

    public int getClassSize(JobClass jobClass) {
        return classSize.getOrDefault(jobClass, 1);
    }

    public void setClassWeight(JobClass jobClass, double weight) {
        classWeight.put(jobClass, weight);
    }

    public double getClassWeight(JobClass jobClass) {
        return classWeight.getOrDefault(jobClass, 1.0);
    }

    public List<Node> getNodes() {
        return nodes;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}

