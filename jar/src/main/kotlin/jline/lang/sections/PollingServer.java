/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.PollingType;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A service section that processes jobs using Polling scheduling
 */
public class PollingServer extends ServiceSection implements Serializable {
    private PollingType pollingType;
    private Map<JobClass, Distribution> switchoverTimes;
    private int pollingK; // K value for K-LIMITED polling

    /**
     * Creates a new polling server for the specified job classes.
     * Initializes with GATED polling type, K=1 for K-LIMITED mode, and zero switchover times.
     * 
     * @param jobClasses the list of job classes this server will handle
     */
    public PollingServer(List<JobClass> jobClasses) {
        super("PollingServer");
        this.numberOfServers = 1;
        this.pollingType = PollingType.GATED; // Default polling type
        this.switchoverTimes = new HashMap<JobClass, Distribution>();
        this.pollingK = 1; // Default K value for K-LIMITED polling
        this.initServers(jobClasses);
    }

    /**
     * Initializes service processes and switchover times for all job classes.
     * Sets exponential(0) as default for both service and switchover times.
     * 
     * @param jobClasses the list of job classes to initialize
     */
    private void initServers(List<JobClass> jobClasses) {
        for (JobClass jobClass : jobClasses) {
            this.serviceProcesses.put(jobClass, new ServiceBinding(jobClass, ServiceStrategy.LI, new Exp(0)));
            // Initialize with zero switchover time by default
            this.switchoverTimes.put(jobClass, new Exp(0));
        }
    }

    /**
     * Sets the polling type for this polling server.
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     */
    public void setPollingType(PollingType pollingType) {
        this.pollingType = pollingType;
    }
    
    /**
     * Sets the polling type for this polling server with K value for K-LIMITED.
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     * @param k the K value for K-LIMITED polling (ignored for other types)
     */
    public void setPollingType(PollingType pollingType, int k) {
        this.pollingType = pollingType;
        if (pollingType == PollingType.KLIMITED) {
            this.pollingK = k;
        }
    }

    /**
     * Gets the current polling type.
     *
     * @return the polling type
     */
    public PollingType getPollingType() {
        return this.pollingType;
    }
    
    /**
     * Gets the K value for K-LIMITED polling.
     *
     * @return the K value
     */
    public int getPollingK() {
        return this.pollingK;
    }
    
    /**
     * Sets the K value for K-LIMITED polling.
     *
     * @param k the K value (must be greater than 0)
     */
    public void setPollingK(int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("K value must be greater than 0");
        }
        this.pollingK = k;
    }

    /**
     * Sets the switchover time for a job class.
     *
     * @param jobClass the job class
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchover(JobClass jobClass, Distribution switchoverTime) {
        this.switchoverTimes.put(jobClass, switchoverTime);
    }

    /**
     * Gets the switchover time for a job class.
     *
     * @param jobClass the job class
     * @return the switchover time distribution
     */
    public Distribution getSwitchover(JobClass jobClass) {
        return this.switchoverTimes.get(jobClass);
    }

}
