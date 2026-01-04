/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;

import java.io.Serializable;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A class of jobs that arrives from the external world to the Network and, after completion, leaves it
 */
public class OpenClass extends JobClass implements Serializable {
    protected Network model;

    /**
     * Creates a new open job class with the specified priority and deadline.
     * Open classes represent jobs that arrive from external sources
     * and leave the system after service completion.
     *
     * @param model    the network model to add this class to
     * @param name     the name for this open class
     * @param priority the priority level for jobs in this class
     * @param deadline the relative deadline in time units from arrival, or Double.POSITIVE_INFINITY for no deadline
     */
    public OpenClass(Network model, String name, int priority, double deadline) {
        super(JobClassType.OPEN, name);
        this.index = model.getNumberOfClasses() + 1;
        model.addJobClass(this);
        try {
            setReferenceStation(model.getSource());
        } catch (Exception e) {
            line_error(mfilename(this), "The model requires a Source prior to instantiating open classes.");
        }

        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node currentNode = model.getNodes().get(i);
            if (currentNode instanceof Join) {
                model.setJoinNodeStrategy(i, this, JoinStrategy.STD);
                model.setJoinNodeRequired(i, this, -1);
            }
            // Skip setting routing for Cache nodes - they must use setProbRouting
            if (currentNode != null && !(currentNode instanceof jline.lang.nodes.Cache)) {
                model.setNodeRouting(i, this, RoutingStrategy.RAND);
            }
        }
        this.model = model;
        this.setPriority(priority);
        this.setDeadline(deadline);
    }

    /**
     * Creates a new open job class with the specified priority.
     * Open classes represent jobs that arrive from external sources
     * and leave the system after service completion.
     *
     * @param model    the network model to add this class to
     * @param name     the name for this open class
     * @param priority the priority level for jobs in this class
     */
    public OpenClass(Network model, String name, int priority) {
        this(model, name, priority, Double.POSITIVE_INFINITY);
    }

    /**
     * Creates a new open job class with default priority (0).
     *
     * @param model the network model to add this class to
     * @param name  the name for this open class
     */
    public OpenClass(Network model, String name) {
        this(model, name, 0, Double.POSITIVE_INFINITY);
    }

    /**
     * Internal constructor for Signal resolution that skips model registration.
     * Package-private: only for use by OpenSignal during Signal.resolve().
     *
     * @param model         the network model (class is NOT added to model)
     * @param name          the name for this open class
     * @param priority      the priority level for jobs in this class
     * @param existingIndex the existing index to preserve from the Signal being replaced
     */
    OpenClass(Network model, String name, int priority, int existingIndex) {
        super(JobClassType.OPEN, name);
        this.index = existingIndex;
        // Skip addJobClass - this class is replacing an existing registered Signal
        try {
            setReferenceStation(model.getSource());
        } catch (Exception e) {
            // Signal classes may not need a Source reference during resolution
        }
        this.model = model;
        this.setPriority(priority);
        this.setDeadline(Double.POSITIVE_INFINITY);
    }

    /**
     * Prints a summary of this open class configuration.
     */
    @Override
    public void printSummary() {
        System.out.format("Open class: %s\n", this.getName());
    }

    /**
     * Sets the reference station for this open class.
     * For open classes, the reference station must be a Source node.
     *
     * @param source the source station to use as reference
     * @throws Exception if the station is not a Source
     */
    @Override
    public void setReferenceStation(Station source) throws Exception {
        if (!(source instanceof Source)) {
            throw new Exception("The reference station for an open class must be a jline.Source.");
        }
        super.setReferenceStation(source);
    }

}
