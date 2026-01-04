/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Station;
import jline.lang.processes.DiscreteDistribution;

import java.io.Serializable;

/**
 * A signal placeholder class that automatically resolves to OpenSignal or ClosedSignal.
 *
 * <p>Signal is a placeholder class that users can use in both open and closed networks.
 * The resolution to the concrete type (OpenSignal or ClosedSignal) happens during model
 * finalization when getStruct() is called. This provides a simplified API where users
 * don't need to explicitly choose between OpenSignal and ClosedSignal.
 *
 * <p>Resolution rules:
 * <ul>
 *   <li>If the network has a Source node: resolves to OpenSignal</li>
 *   <li>If the network has no Source node: resolves to ClosedSignal</li>
 * </ul>
 *
 * <p>For explicit control, users can directly use {@link OpenSignal} or {@link ClosedSignal}.
 *
 * <p>Signal types:
 * <ul>
 *   <li>NEGATIVE: Removes a job from the destination queue</li>
 *   <li>CATASTROPHE: Resets the destination queue to empty state</li>
 *   <li>REPLY: Triggers a reply action (unblocks servers waiting for reply)</li>
 * </ul>
 *
 * <p>Reference: Gelenbe, E. (1991). "Product-form queueing networks with
 * negative and positive customers", Journal of Applied Probability
 *
 * @see OpenSignal
 * @see ClosedSignal
 * @see SignalType
 */
public class Signal extends JobClass implements Serializable {

    private SignalType signalType;
    private JobClass targetJobClass;
    private DiscreteDistribution removalDistribution;
    private RemovalPolicy removalPolicy;
    private Network model;

    /**
     * Creates a new signal placeholder with full configuration for batch removal.
     *
     * @param model               the network model to add this class to
     * @param name                the name for this signal class
     * @param signalType          the type of signal (NEGATIVE, CATASTROPHE, or REPLY)
     * @param priority            the priority level for signals in this class
     * @param removalDistribution the distribution for number of jobs to remove (null for exactly 1)
     * @param removalPolicy       the policy for selecting which jobs to remove
     */
    public Signal(Network model, String name, SignalType signalType, int priority,
                  DiscreteDistribution removalDistribution, RemovalPolicy removalPolicy) {
        super(JobClassType.OPEN, name);  // Default to OPEN, resolved later
        // Set index before adding to model (1-based, like OpenClass)
        this.index = model.getNumberOfClasses() + 1;
        this.setPriority(priority);
        this.model = model;
        this.signalType = signalType;
        this.targetJobClass = null;
        this.removalDistribution = removalDistribution;
        this.removalPolicy = removalPolicy != null ? removalPolicy : RemovalPolicy.RANDOM;
        model.addJobClass(this);
    }

    /**
     * Creates a new signal placeholder with the specified type and priority.
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE, CATASTROPHE, or REPLY)
     * @param priority   the priority level for signals in this class
     */
    public Signal(Network model, String name, SignalType signalType, int priority) {
        this(model, name, signalType, priority, null, RemovalPolicy.RANDOM);
    }

    /**
     * Creates a new signal placeholder with the specified type and default priority (0).
     * Uses default removal behavior (remove exactly 1 job, random selection).
     *
     * @param model      the network model to add this class to
     * @param name       the name for this signal class
     * @param signalType the type of signal (NEGATIVE, CATASTROPHE, or REPLY) - REQUIRED
     */
    public Signal(Network model, String name, SignalType signalType) {
        this(model, name, signalType, 0, null, RemovalPolicy.RANDOM);
    }

    /**
     * Resolves this Signal placeholder to OpenSignal or ClosedSignal.
     *
     * @param isOpen  true if the network is open (has Source node)
     * @param refstat reference station for closed networks (ignored for open)
     * @return the resolved OpenSignal or ClosedSignal instance
     */
    public JobClass resolve(boolean isOpen, Station refstat) {
        JobClass concrete;
        int existingIndex = this.getIndex();
        if (isOpen) {
            // Use internal constructor that skips model registration
            OpenSignal openSignal = new OpenSignal(model, getName(), signalType, getPriority(),
                    removalDistribution, removalPolicy, existingIndex);
            if (targetJobClass != null) {
                openSignal.forJobClass(targetJobClass);
            }
            concrete = openSignal;
        } else {
            // Use internal constructor that skips model registration
            ClosedSignal closedSignal = new ClosedSignal(model, getName(), signalType, refstat,
                    getPriority(), removalDistribution, removalPolicy, existingIndex);
            if (targetJobClass != null) {
                closedSignal.forJobClass(targetJobClass);
            }
            concrete = closedSignal;
        }
        return concrete;
    }

    /**
     * Gets the network model this signal belongs to.
     *
     * @return the network model
     */
    public Network getModel() {
        return model;
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
     * @return this Signal instance (for method chaining)
     */
    public Signal forJobClass(JobClass jobClass) {
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
     * Prints a summary of this signal class configuration.
     */
    @Override
    public void printSummary() {
        String removalInfo = "";
        if (removalDistribution != null) {
            removalInfo = String.format(", removal: %s", removalDistribution.getName());
        }
        System.out.format("Signal class (placeholder): %s (type: %s, policy: %s%s)\n",
            this.getName(), SignalType.toText(this.signalType),
            RemovalPolicy.toText(this.removalPolicy), removalInfo);
    }
}
