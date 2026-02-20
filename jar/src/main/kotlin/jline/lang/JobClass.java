/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import static jline.GlobalConstants.Inf;

import jline.lang.constant.BalkingStrategy;
import jline.lang.constant.BalkingThreshold;
import jline.lang.constant.JobClassType;
import jline.lang.constant.ImpatienceType;
import jline.lang.processes.Distribution;
import jline.lang.processes.BMAP;
import jline.lang.processes.MAP;
import jline.lang.processes.MMPP2;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Superclass representing a class of jobs
 */
public class JobClass extends NetworkElement implements Serializable {
    protected JobClassType type;
    protected int priority;
    protected boolean completes;
    protected Station refstat;
    protected boolean isrefclass;
    protected int index;
    protected double deadline;

    private Integer[] attribute;

    /**
     * Index of the reply signal class expected after processing jobs of this class.
     * When set to a non-negative value, the server will block after sending a job
     * forward until receiving a REPLY signal from the specified class.
     * Value of -1 indicates no reply is expected (normal processing).
     */
    private int replySignalClassIndex = -1;

    /**
     * Global patience distribution for this job class.
     * Customers of this class will abandon queues after waiting for a time
     * drawn from this distribution. This setting applies to all queues unless
     * overridden by a queue-specific patience setting.
     */
    private Distribution patience;

    /**
     * Global impatience type for this job class (RENEGING or BALKING).
     * This setting applies to all queues unless overridden by a queue-specific setting.
     */
    private ImpatienceType impatienceType;

    /**
     * Global balking strategy for this job class.
     * Determines how balking decisions are made (queue length, expected wait, or combined).
     */
    private BalkingStrategy balkingStrategy;

    /**
     * Global balking thresholds for queue-length based balking.
     * Each threshold maps a queue length range to a balking probability.
     */
    private List<BalkingThreshold> balkingThresholds;

    /**
     * Global retrial delay distribution for this job class.
     * When a customer is rejected due to capacity and retrial is enabled,
     * they wait for a time drawn from this distribution before retrying.
     */
    private Distribution retrialDelayDistribution;

    /**
     * Maximum number of retrial attempts for this job class.
     * Value of -1 indicates unlimited retries.
     * Value of 0 means no retries (effectively a drop).
     * Positive values limit the number of retry attempts.
     */
    private int maxRetrialAttempts = -1;

    /**
     * Global immediate feedback setting for this job class.
     * When true, jobs of this class that self-loop at any station stay in service
     * instead of rejoining the queue.
     */
    private boolean immediateFeedback = false;

    /**
     * Creates a new job class with the specified type and name.
     *
     * @param type the type of job class (OPEN, CLOSED, etc.)
     * @param name the name for this job class
     */
    public JobClass(JobClassType type, String name) {
        super(name);
        this.priority = 0;
        this.deadline = Double.POSITIVE_INFINITY;
        this.refstat = null;
        this.type = type;
        this.completes = true;
        this.isrefclass = false;
        this.attribute = new Integer[]{null, null};
    }

    public Integer[] getAttribute() {
        return attribute;
    }

    public void setAttribute(Integer[] attribute) {
        this.attribute = attribute;
    }

    /**
     * Checks if jobs in this class complete their service and leave the system.
     *
     * @return true if jobs complete, false if they remain in the system
     */
    public boolean getCompletes() {
        return this.completes;
    }

    public void setCompletes(boolean completes) {
        this.completes = completes;
    }

    public int getIndex() {
        return index;
    }

    /**
     * Returns the type of this job class.
     *
     * @return the job class type (OPEN, CLOSED, etc.)
     */
    public JobClassType getJobClassType() {
        return this.type;
    }

    /**
     * Returns the number of jobs in this class.
     * Default implementation returns infinite (for open classes).
     *
     * @return number of jobs in this class
     */
    public double getNumberOfJobs() {
        return Inf;
    }

    /**
     * Returns the priority level of this job class.
     * Lower values indicate higher priority (priority 0 is highest).
     *
     * @return the priority level
     */
    public int getPriority() {
        return this.priority;
    }

    /**
     * Sets the priority level for this job class.
     *
     * @param p the priority level to set
     */
    public void setPriority(int p) {
        this.priority = p;
    }

    /**
     * Returns the relative deadline for jobs in this class.
     * The deadline is specified as time units from job arrival.
     * A value of Double.POSITIVE_INFINITY indicates no deadline constraint.
     *
     * @return the relative deadline in time units, or Double.POSITIVE_INFINITY if no deadline
     */
    public double getDeadline() {
        return this.deadline;
    }

    /**
     * Sets the relative deadline for jobs in this class.
     * Jobs will be assigned an absolute deadline of arrivalTime + deadline
     * when they enter queues with EDD/EDF scheduling.
     *
     * @param deadline the relative deadline in time units (must be > 0), or Double.POSITIVE_INFINITY for no deadline
     */
    public void setDeadline(double deadline) {
        if (deadline <= 0.0 && !Double.isInfinite(deadline)) {
            throw new IllegalArgumentException("Deadline must be positive or Infinity (no deadline)");
        }
        this.deadline = deadline;
    }

    /**
     * Returns the reference station for this job class.
     * The reference station is used for normalization in closed networks.
     *
     * @return the reference station, or null if not set
     */
    public Station getReferenceStation() {
        return this.refstat;
    }

    /**
     * Sets the reference station for this job class.
     *
     * @param ref the station to use as reference
     * @throws Exception if the reference station cannot be set
     */
    public void setReferenceStation(Station ref) throws Exception {
        this.refstat = ref;
    }

    /**
     * Checks if this is a reference job class within its chain.
     *
     * @return true if this is a reference class
     */
    public boolean isReferenceClass() {
        return this.isrefclass;
    }

    /**
     * Sets whether this job class is a reference class within its chain.
     *
     * @param isrefclass true to make this a reference class
     */
    public void setReferenceClass(boolean isrefclass) {
        this.isrefclass = isrefclass;
    }

    /**
     * Checks if the specified node is the reference station for this class.
     *
     * @param node the node to check
     * @return true if the node is the reference station
     */
    public boolean isReferenceStation(Node node) {
        return name.equals(node.getName());
    }

    /**
     * Sets the reply signal class index for synchronous call semantics.
     * When a job of this class completes service, the server will block
     * until receiving a REPLY signal from the signal class at the specified index.
     *
     * @param signalClassIndex the index of the signal class that will unblock the server
     */
    public void setReplySignalClassIndex(int signalClassIndex) {
        this.replySignalClassIndex = signalClassIndex;
    }

    /**
     * Gets the index of the reply signal class.
     *
     * @return the index of the reply signal class, or -1 if no reply is expected
     */
    public int getReplySignalClassIndex() {
        return this.replySignalClassIndex;
    }

    /**
     * Checks if this job class expects a reply signal after processing.
     *
     * @return true if a reply signal is expected, false otherwise
     */
    public boolean expectsReply() {
        return this.replySignalClassIndex >= 0;
    }

    /**
     * Sets the global patience distribution for this job class.
     * This applies to all queues unless overridden by a queue-specific patience setting.
     * Defaults to RENEGING patience type for backwards compatibility.
     *
     * @param distribution the patience time distribution
     * @throws IllegalArgumentException if distribution is a modulated process or argument is null
     */
    public void setPatience(Distribution distribution) {
        setPatience(ImpatienceType.RENEGING, distribution);
    }

    /**
     * Sets the global impatience type and distribution for this job class.
     * This applies to all queues unless overridden by a queue-specific patience setting.
     *
     * @param impatienceType the type of impatience (RENEGING or BALKING)
     * @param distribution the patience time distribution
     * @throws IllegalArgumentException if distribution is a modulated process or arguments are null
     * @throws UnsupportedOperationException if BALKING type is used (not yet supported)
     */
    public void setPatience(ImpatienceType impatienceType, Distribution distribution) {
        if (impatienceType == null) {
            throw new IllegalArgumentException("impatienceType cannot be null");
        }
        if (distribution == null) {
            throw new IllegalArgumentException("distribution cannot be null");
        }

        // Validate impatience type
        if (impatienceType != ImpatienceType.RENEGING && impatienceType != ImpatienceType.BALKING) {
            throw new IllegalArgumentException("Invalid impatience type. Use ImpatienceType.RENEGING or ImpatienceType.BALKING.");
        }

        // Only RENEGING is currently supported
        if (impatienceType == ImpatienceType.BALKING) {
            throw new UnsupportedOperationException("BALKING impatience type is not yet supported. Use ImpatienceType.RENEGING.");
        }

        // Validate distribution type
        if (distribution instanceof BMAP || distribution instanceof MAP ||
            distribution instanceof MMPP2) {
            throw new IllegalArgumentException(
                "Modulated processes (BMAP, MAP, MMPP2) are not supported for patience distributions.");
        }

        this.patience = distribution;
        this.impatienceType = impatienceType;
    }

    /**
     * Gets the global patience distribution for this job class.
     *
     * @return the patience distribution, or null if not set
     */
    public Distribution getPatience() {
        return this.patience;
    }

    /**
     * Gets the global impatience type for this job class.
     *
     * @return the impatience type, or null if not set
     */
    public ImpatienceType getImpatienceType() {
        return this.impatienceType;
    }

    /**
     * Checks if this job class has a patience distribution configured.
     *
     * @return true if patience is configured, false otherwise
     */
    public boolean hasPatience() {
        return this.patience != null && !this.patience.isDisabled();
    }

    // =================== BALKING METHODS ===================

    /**
     * Sets the global balking configuration for this job class using queue-length thresholds.
     * This applies to all queues unless overridden by a queue-specific balking setting.
     *
     * @param strategy the balking strategy to use
     * @param thresholds list of balking thresholds mapping queue lengths to probabilities
     * @throws IllegalArgumentException if strategy or thresholds are null
     */
    public void setBalking(BalkingStrategy strategy, List<BalkingThreshold> thresholds) {
        if (strategy == null) {
            throw new IllegalArgumentException("BalkingStrategy cannot be null");
        }
        if (thresholds == null || thresholds.isEmpty()) {
            throw new IllegalArgumentException("Balking thresholds cannot be null or empty");
        }
        this.balkingStrategy = strategy;
        this.balkingThresholds = new ArrayList<BalkingThreshold>(thresholds);
    }

    /**
     * Sets the global balking configuration using a single threshold.
     * Convenience method for simple balking configurations.
     *
     * @param strategy the balking strategy to use
     * @param minJobs minimum queue length to trigger balking
     * @param probability balking probability when queue length >= minJobs
     */
    public void setBalking(BalkingStrategy strategy, int minJobs, double probability) {
        List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
        thresholds.add(new BalkingThreshold(minJobs, probability));
        setBalking(strategy, thresholds);
    }

    /**
     * Gets the global balking strategy for this job class.
     *
     * @return the balking strategy, or null if not configured
     */
    public BalkingStrategy getBalkingStrategy() {
        return this.balkingStrategy;
    }

    /**
     * Gets the global balking thresholds for this job class.
     *
     * @return list of balking thresholds, or null if not configured
     */
    public List<BalkingThreshold> getBalkingThresholds() {
        return this.balkingThresholds;
    }

    /**
     * Checks if this job class has balking configured.
     *
     * @return true if balking is configured, false otherwise
     */
    public boolean hasBalking() {
        return this.balkingStrategy != null && this.balkingThresholds != null && !this.balkingThresholds.isEmpty();
    }

    // =================== RETRIAL METHODS ===================

    /**
     * Sets the global retrial configuration for this job class.
     * When a customer is rejected due to capacity, they will enter an orbit
     * and retry after a delay drawn from the specified distribution.
     *
     * @param delayDistribution the distribution for retrial delay times
     * @param maxAttempts maximum retry attempts (-1 for unlimited, 0 for no retries)
     * @throws IllegalArgumentException if delayDistribution is null or a modulated process
     */
    public void setRetrial(Distribution delayDistribution, int maxAttempts) {
        if (delayDistribution == null) {
            throw new IllegalArgumentException("Retrial delay distribution cannot be null");
        }
        if (delayDistribution instanceof BMAP || delayDistribution instanceof MAP ||
            delayDistribution instanceof MMPP2) {
            throw new IllegalArgumentException(
                "Modulated processes (BMAP, MAP, MMPP2) are not supported for retrial distributions.");
        }
        this.retrialDelayDistribution = delayDistribution;
        this.maxRetrialAttempts = maxAttempts;
    }

    /**
     * Sets the global retrial configuration with unlimited retry attempts.
     *
     * @param delayDistribution the distribution for retrial delay times
     */
    public void setRetrial(Distribution delayDistribution) {
        setRetrial(delayDistribution, -1);
    }

    /**
     * Gets the global retrial delay distribution for this job class.
     *
     * @return the retrial delay distribution, or null if not configured
     */
    public Distribution getRetrialDelayDistribution() {
        return this.retrialDelayDistribution;
    }

    /**
     * Gets the maximum number of retrial attempts for this job class.
     *
     * @return maximum attempts (-1 for unlimited)
     */
    public int getMaxRetrialAttempts() {
        return this.maxRetrialAttempts;
    }

    /**
     * Checks if this job class has retrial configured.
     *
     * @return true if retrial is configured, false otherwise
     */
    public boolean hasRetrial() {
        return this.retrialDelayDistribution != null && !this.retrialDelayDistribution.isDisabled();
    }

    /**
     * Prints a summary of this job class configuration to the console.
     */
    public void printSummary() {
        System.out.format("Job Class: %s\n", this.getName());
        System.out.println();
    }

    /*
     * ===================================================================================
     * MISSING METHODS FROM MATLAB JOBCLASS IMPLEMENTATION - NOT YET MIGRATED
     * ===================================================================================
     *
     * Based on analysis of /matlab/src/lang/JobClass.m
     */

    // =================== MATLAB-SPECIFIC METHODS ===================
    // public int subsindex()  // MATLAB indexing method (0-based index for array subscripts)
    // public void summary()   // Enhanced summary with formatted output

    // =================== IMMEDIATE FEEDBACK METHODS ===================

    /**
     * Enables or disables immediate feedback for this job class globally.
     * When enabled, jobs of this class that self-loop at any station stay in service
     * instead of rejoining the queue.
     *
     * @param enabled true to enable immediate feedback, false to disable
     */
    public void setImmediateFeedback(boolean enabled) {
        this.immediateFeedback = enabled;
    }

    /**
     * Checks if immediate feedback is enabled for this job class.
     *
     * @return true if immediate feedback is enabled
     */
    public boolean hasImmediateFeedback() {
        return this.immediateFeedback;
    }

}
