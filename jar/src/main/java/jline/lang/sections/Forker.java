/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Node;

import java.util.List;

/**
 * Output section that forks incoming jobs into sibling tasks
 */
public class Forker extends OutputSection {
    /**
     * Nominal tasks emitted on each outgoing link, shared by every link and
     * every class. For a fork whose degree is random this is the MEAN, so a
     * consumer that only knows this field gets E[tasks per link] rather than a
     * number the fork never emits.
     */
    public double tasksPerLink;

    /**
     * One variable-forking-level override. The destination is kept by NAME
     * rather than by link ordinal, because the link order is an artefact of
     * connmatrix traversal and renumbers whenever the model is relinked; an
     * empty name means every connected link of that class.
     */
    public static class ForkOverride {
        public final String dest;
        public final int jobClass;      // 1-based, as JobClass.getIndex() returns
        public final double value;      // tasks per link, or branch probability
        public final jline.lang.processes.DiscreteSampler dist;  // null unless random

        public ForkOverride(String dest, int jobClass, double value) {
            this(dest, jobClass, value, null);
        }

        public ForkOverride(String dest, int jobClass, double value,
                            jline.lang.processes.DiscreteSampler dist) {
            this.dest = dest;
            this.jobClass = jobClass;
            this.value = value;
            this.dist = dist;
        }
    }

    /** Per-destination, per-class tasks-per-link overrides. */
    public List<ForkOverride> tasksPerLinkByDest = new java.util.ArrayList<ForkOverride>();
    /** Per-destination, per-class jobs-per-link distributions. */
    public List<ForkOverride> tasksPerLinkDist = new java.util.ArrayList<ForkOverride>();
    /** Per-destination, per-class branch activation probabilities. */
    public List<ForkOverride> branchProb = new java.util.ArrayList<ForkOverride>();

    protected List<JobClass> jobClasses;

    public Forker(List<JobClass> customerClasses) {
        super("Forker");
        this.jobClasses = customerClasses;
        this.tasksPerLink = 1.0;
        this.initDispatcherJobClasses(customerClasses);
    }

    private void initDispatcherJobClasses(List<JobClass> customerClasses) {
        for (JobClass jobClass : customerClasses) {
            this.outputStrategies.add(new OutputStrategy(jobClass, RoutingStrategy.RAND));
        }
    }

    @Override
    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if ((outputStrategy.getJobClass().getIndex() == jobClass.getIndex()) &&
                (outputStrategy.getRoutingStrategy() != RoutingStrategy.PROB ||
                 outputStrategy.getDestination() == null ||
                 outputStrategy.getDestination().getNodeIndex() == destination.getNodeIndex())) {
                outputStrategy.setRoutingStrategy(routingStrategy);
                outputStrategy.setDestination(destination);
                outputStrategy.setProbability(probability);
                this.probabilityUpdate();
                return;
            }
        }

        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy, destination, probability);
        outputStrategies.add(outputStrategy);
        this.probabilityUpdate();
    }
}
