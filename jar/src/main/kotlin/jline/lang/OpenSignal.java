/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SignalType;
import jline.lang.processes.DiscreteDistribution;

import java.io.Serializable;

/**
 * An open signal class representing special customers in G-networks and related models.
 *
 * <p>OpenSignal is a specialized OpenClass that can have special effects on queues
 * they visit. Unlike regular customers that simply receive service, signals can:
 * <ul>
 *   <li>Remove jobs from destination queues (NEGATIVE signals)</li>
 *   <li>Reset queue states (CATASTROPHE signals)</li>
 *   <li>Trigger reply actions (REPLY signals)</li>
 * </ul>
 *
 * <p>For closed networks, use {@link ClosedSignal} instead.
 *
 * <p>Reference: Gelenbe, E. (1991). "Product-form queueing networks with
 * negative and positive customers", Journal of Applied Probability
 *
 * @see OpenClass
 * @see ClosedSignal
 * @see SignalType
 */
public class OpenSignal extends OpenClass implements Serializable {

    private SignalType signalType;
    private JobClass targetJobClass;
    private DiscreteDistribution removalDistribution;
    private RemovalPolicy removalPolicy;

    /**
     * Creates a new open signal class with full configuration for batch removal.
     *
     * @param model               the network model to add this class to
     * @param name                the name for this signal class
     * @param signalType          the type of signal (NEGATIVE, CATASTROPHE, or REPLY)
     * @param priority            the priority level for signals in this class
     * @param removalDistribution the distribution for number of jobs to remove (null for exactly 1)
     * @param removalPolicy       the policy for selecting which jobs to remove
     */
    public OpenSignal(Network model, String name, SignalType signalType, int priority,
                      DiscreteDistribution removalDistribution, RemovalPolicy removalPolicy) {
        super(model, name, priority);
        this.signalType = signalType;
        this.targetJobClass = null;
        this.removalDistribution = removalDistribution;
        this.removalPolicy = removalPolicy != null ? removalPolicy : RemovalPolicy.RANDOM;
    }

    /**
     * Creates a new open signal class with the specified type and priority.
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE, CATASTROPHE, or REPLY)
     * @param priority   the priority level for signals in this class
     */
    public OpenSignal(Network model, String name, SignalType signalType, int priority) {
        this(model, name, signalType, priority, null, RemovalPolicy.RANDOM);
    }

    /**
     * Creates a new open signal class with the specified type and default priority (0).
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE, CATASTROPHE, or REPLY) - REQUIRED
     */
    public OpenSignal(Network model, String name, SignalType signalType) {
        this(model, name, signalType, 0, null, RemovalPolicy.RANDOM);
    }

    /**
     * Internal constructor for Signal resolution that skips model registration.
     * Package-private: only for use by Signal.resolve().
     *
     * @param model               the network model (class is NOT added to model)
     * @param name                the name for this signal class
     * @param signalType          the type of signal (NEGATIVE, CATASTROPHE, or REPLY)
     * @param priority            the priority level for signals in this class
     * @param removalDistribution the distribution for number of jobs to remove (null for exactly 1)
     * @param removalPolicy       the policy for selecting which jobs to remove
     * @param existingIndex       the existing index to preserve from the Signal being replaced
     */
    OpenSignal(Network model, String name, SignalType signalType, int priority,
               DiscreteDistribution removalDistribution, RemovalPolicy removalPolicy,
               int existingIndex) {
        super(model, name, priority, existingIndex);  // Uses internal OpenClass constructor
        this.signalType = signalType;
        this.targetJobClass = null;
        this.removalDistribution = removalDistribution;
        this.removalPolicy = removalPolicy != null ? removalPolicy : RemovalPolicy.RANDOM;
    }

    /**
     * Gets the signal type for this signal class.
     *
     * @return the signal type
     */
    public SignalType getSignalType() {
        return signalType;
    }

    /**
     * Sets the signal type for this signal class.
     *
     * @param signalType the new signal type
     */
    public void setSignalType(SignalType signalType) {
        this.signalType = signalType;
    }

    /**
     * Associates this signal with a job class.
     *
     * <p>For REPLY signals, this specifies which job class's servers
     * will be unblocked when this signal arrives.
     *
     * @param jobClass the JobClass to associate with this signal
     * @return this OpenSignal instance (for method chaining)
     */
    public OpenSignal forJobClass(JobClass jobClass) {
        this.targetJobClass = jobClass;
        if (jobClass != null) {
            jobClass.setReplySignalClassIndex(this.getIndex());
        }
        return this;
    }

    /**
     * Gets the associated job class.
     *
     * @return the JobClass associated with this signal, or null if none
     */
    public JobClass getTargetJobClass() {
        return targetJobClass;
    }

    /**
     * Gets the index of the associated job class.
     *
     * @return the index of the associated JobClass, or -1 if none
     */
    public int getTargetJobClassIndex() {
        if (targetJobClass == null) {
            return -1;
        }
        return targetJobClass.getIndex();
    }

    /**
     * Gets the removal distribution for this negative signal.
     *
     * <p>The removal distribution determines how many positive customers
     * are removed when this signal arrives at a queue. If null, exactly
     * one customer is removed (default G-network behavior).
     *
     * @return the removal distribution, or null for single removal
     */
    public DiscreteDistribution getRemovalDistribution() {
        return removalDistribution;
    }

    /**
     * Sets the removal distribution for this negative signal.
     *
     * @param removalDistribution the distribution for number of removals
     */
    public void setRemovalDistribution(DiscreteDistribution removalDistribution) {
        this.removalDistribution = removalDistribution;
    }

    /**
     * Gets the removal policy for this negative signal.
     *
     * <p>The removal policy determines which customers are selected
     * for removal when this signal arrives at a queue.
     *
     * @return the removal policy (RANDOM, FCFS, or LCFS)
     */
    public RemovalPolicy getRemovalPolicy() {
        return removalPolicy;
    }

    /**
     * Sets the removal policy for this negative signal.
     *
     * @param removalPolicy the policy for selecting customers to remove
     */
    public void setRemovalPolicy(RemovalPolicy removalPolicy) {
        this.removalPolicy = removalPolicy;
    }

    /**
     * Returns whether this signal is a catastrophe (removes all jobs).
     *
     * @return true if signalType is CATASTROPHE
     */
    public boolean isCatastrophe() {
        return signalType == SignalType.CATASTROPHE;
    }

    /**
     * Prints a summary of this open signal class configuration.
     */
    @Override
    public void printSummary() {
        String removalInfo = "";
        if (removalDistribution != null) {
            removalInfo = String.format(", removal: %s", removalDistribution.getName());
        }
        System.out.format("OpenSignal class: %s (type: %s, policy: %s%s)\n",
            this.getName(), SignalType.toText(this.signalType),
            RemovalPolicy.toText(this.removalPolicy), removalInfo);
    }
}
