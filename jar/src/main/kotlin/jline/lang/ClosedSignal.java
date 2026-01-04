/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Station;
import jline.lang.processes.DiscreteDistribution;

import java.io.Serializable;

/**
 * A closed signal class for modeling signals in closed queueing networks.
 *
 * <p>ClosedSignal is a specialized ClosedClass for signals that need to circulate
 * in closed networks. Unlike regular Signal (which extends OpenClass), ClosedSignal
 * can be used in networks without Source/Sink nodes.
 *
 * <p>ClosedSignal has zero population - signals are created dynamically through
 * class switching from the target job class.
 *
 * <p>Signal types:
 * <ul>
 *   <li>NEGATIVE: Removes a job from the destination queue</li>
 *   <li>REPLY: Unblocks servers waiting for a reply (LQN synchronous call semantics)</li>
 * </ul>
 *
 * @see ClosedClass
 * @see Signal
 * @see SignalType
 */
public class ClosedSignal extends ClosedClass implements Serializable {

    private SignalType signalType;
    private JobClass targetJobClass;
    private DiscreteDistribution removalDistribution;
    private RemovalPolicy removalPolicy;

    /**
     * Creates a new closed signal class with full configuration for batch removal.
     *
     * @param model               the network model to add this class to
     * @param name                the name for this signal class
     * @param signalType          the type of signal (NEGATIVE or REPLY)
     * @param refstat             the reference station (should match target job class)
     * @param priority            the priority level for signals in this class
     * @param removalDistribution the distribution for number of jobs to remove (null for exactly 1)
     * @param removalPolicy       the policy for selecting which jobs to remove
     */
    public ClosedSignal(Network model, String name, SignalType signalType, Station refstat, int priority,
                        DiscreteDistribution removalDistribution, RemovalPolicy removalPolicy) {
        super(model, name, 0, refstat, priority);  // 0 population - signals are created by class switching
        this.signalType = signalType;
        this.targetJobClass = null;
        this.removalDistribution = removalDistribution;
        this.removalPolicy = removalPolicy != null ? removalPolicy : RemovalPolicy.RANDOM;
    }

    /**
     * Creates a new closed signal class with the specified type, reference station, and priority.
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE or REPLY)
     * @param refstat    the reference station (should match target job class)
     * @param priority   the priority level for signals in this class
     */
    public ClosedSignal(Network model, String name, SignalType signalType, Station refstat, int priority) {
        this(model, name, signalType, refstat, priority, null, RemovalPolicy.RANDOM);
    }

    /**
     * Creates a new closed signal class with the specified type, reference station, and default priority (0).
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE or REPLY)
     * @param refstat    the reference station (should match target job class)
     */
    public ClosedSignal(Network model, String name, SignalType signalType, Station refstat) {
        this(model, name, signalType, refstat, 0, null, RemovalPolicy.RANDOM);
    }

    /**
     * Creates a new negative closed signal class with default priority (0).
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model   the network model to add this class to
     * @param name    the name for this signal class
     * @param refstat the reference station
     */
    public ClosedSignal(Network model, String name, Station refstat) {
        this(model, name, SignalType.NEGATIVE, refstat, 0, null, RemovalPolicy.RANDOM);
    }

    /**
     * Internal constructor for Signal resolution that skips model registration.
     * Package-private: only for use by Signal.resolve().
     *
     * @param model               the network model (class is NOT added to model)
     * @param name                the name for this signal class
     * @param signalType          the type of signal (NEGATIVE or REPLY)
     * @param refstat             the reference station (should match target job class)
     * @param priority            the priority level for signals in this class
     * @param removalDistribution the distribution for number of jobs to remove (null for exactly 1)
     * @param removalPolicy       the policy for selecting which jobs to remove
     * @param existingIndex       the existing index to preserve from the Signal being replaced
     */
    ClosedSignal(Network model, String name, SignalType signalType, Station refstat, int priority,
                 DiscreteDistribution removalDistribution, RemovalPolicy removalPolicy,
                 int existingIndex) {
        super(model, name, 0, refstat, priority, existingIndex, true);  // Uses internal ClosedClass constructor
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
     * @return this ClosedSignal instance (for method chaining)
     */
    public ClosedSignal forJobClass(JobClass jobClass) {
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
     * Prints a summary of this closed signal class configuration.
     */
    @Override
    public void printSummary() {
        String removalInfo = "";
        if (removalDistribution != null) {
            removalInfo = String.format(", removal: %s", removalDistribution.getName());
        }
        System.out.format("ClosedSignal class: %s (type: %s, target: %s, policy: %s%s)\n",
            this.getName(), SignalType.toText(this.signalType),
            targetJobClass != null ? targetJobClass.getName() : "none",
            RemovalPolicy.toText(this.removalPolicy), removalInfo);
    }
}
