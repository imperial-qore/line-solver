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
    public double tasksPerLink;
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
                 (outputStrategy.getDestination() != null && outputStrategy.getDestination().getNodeIndex() == destination.getNodeIndex()))) {
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
