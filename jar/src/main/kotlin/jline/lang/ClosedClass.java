/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

import java.io.Serializable;

/**
 * Class where jobs perpetually loop without arriving or leaving (Closed class)
 */
public class ClosedClass extends JobClass implements Serializable {
    protected double population;
    protected Network model;

    /**
     * Creates a new closed job class with the specified population, reference station, priority, and deadline.
     * Closed classes have a fixed number of jobs that circulate within the network.
     *
     * @param model    the network model to add this class to
     * @param name     the name for this closed class
     * @param njobs    the fixed population of jobs in this class
     * @param refstat  the reference station for performance normalization
     * @param priority the priority level for jobs in this class
     * @param deadline the relative deadline in time units from arrival, or Double.POSITIVE_INFINITY for no deadline
     */
    public ClosedClass(Network model, String name, double njobs, Station refstat, int priority, double deadline) {
        super(JobClassType.CLOSED, name);
        this.index = model.getNumberOfClasses() + 1;
        model.addJobClass(this);
        this.population = njobs;
        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node currentNode = model.getNodes().get(i);
//            model.setNodeRouting(i, this, RoutingStrategy.RAND);
            if (currentNode instanceof Join) {
                model.setJoinNodeStrategy(i, this, JoinStrategy.STD);
                model.setJoinNodeRequired(i, this, -1);
            }
        }
        this.model = model;
        this.refstat = refstat;
        this.setPriority(priority);
        this.setDeadline(deadline);
    }

    public ClosedClass(Network model, String name, int njobs, Station refstat, int priority, double deadline) {
        this(model, name, (double) njobs, refstat, priority, deadline);
    }

    /**
     * Creates a new closed job class with the specified population and reference station.
     * Closed classes have a fixed number of jobs that circulate within the network.
     *
     * @param model    the network model to add this class to
     * @param name     the name for this closed class
     * @param njobs    the fixed population of jobs in this class
     * @param refstat  the reference station for performance normalization
     * @param priority the priority level for jobs in this class
     */
    public ClosedClass(Network model, String name, double njobs, Station refstat, int priority) {
        this(model, name, njobs, refstat, priority, Double.POSITIVE_INFINITY);
    }

    public ClosedClass(Network model, String name, int njobs, Station refstat, int priority) {
        this(model, name, (double) njobs, refstat, priority, Double.POSITIVE_INFINITY);
    }

    /**
     * Creates a new closed job class with default priority (0).
     *
     * @param model   the network model to add this class to
     * @param name    the name for this closed class
     * @param njobs   the fixed population of jobs in this class
     * @param refstat the reference station for performance normalization
     */
    public ClosedClass(Network model, String name, double njobs, Station refstat) {
        this(model, name, njobs, refstat, 0, Double.POSITIVE_INFINITY);
    }

    public ClosedClass(Network model, String name, int njobs, Station refstat) {
        this(model, name, (double) njobs, refstat, 0, Double.POSITIVE_INFINITY);
    }

    /**
     * Internal constructor for Signal resolution that skips model registration.
     * Package-private: only for use by ClosedSignal during Signal.resolve().
     *
     * @param model         the network model (class is NOT added to model)
     * @param name          the name for this closed class
     * @param njobs         the fixed population of jobs in this class
     * @param refstat       the reference station for performance normalization
     * @param priority      the priority level for jobs in this class
     * @param existingIndex the existing index to preserve from the Signal being replaced
     * @param internal      marker parameter to distinguish from public constructors (always true)
     */
    ClosedClass(Network model, String name, double njobs, Station refstat, int priority, int existingIndex, boolean internal) {
        super(JobClassType.CLOSED, name);
        this.index = existingIndex;
        // Skip addJobClass - this class is replacing an existing registered Signal
        this.population = njobs;
        this.model = model;
        this.refstat = refstat;
        this.setPriority(priority);
        this.setDeadline(Double.POSITIVE_INFINITY);
    }

    /**
     * Returns the fixed population of jobs in this closed class.
     *
     * @return the number of jobs in this class
     */
    @Override
    public double getNumberOfJobs() {
        return population;
    }

    /**
     * Returns the population size for this closed class.
     *
     * @return the population size
     */
    public double getPopulation() {
        return this.population;
    }

    /**
     * Sets the population size for this closed class.
     *
     * @param pop the new population size
     */
    public void setPopulation(double pop) {
        this.population = pop;
    }

    /**
     * Returns the reference station for this closed class.
     *
     * @return the reference station
     */
    public Station getReferenceStation() {
        return this.refstat;
    }

    /**
     * Prints a summary of this closed class configuration.
     */
    @Override
    public void printSummary() {
        System.out.format("Closed class: %s\n", this.getName());
    }

}
