/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Node;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Class modelling the output section of a Node
 */
public class OutputStrategy implements Serializable {
    public static List<RoutingStrategy> legalStrategies = Arrays.asList(RoutingStrategy.DISABLED, RoutingStrategy.PROB, RoutingStrategy.RAND, RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ);
    private final JobClass jobClass;
    private RoutingStrategy routingStrategy;
    private double probability;
    private Node destination;

    public OutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        this.jobClass = jobClass;
        this.routingStrategy = routingStrategy;
        this.destination = destination;
        this.probability = probability;

        if (!legalStrategies.contains(routingStrategy)) {
            // Note: legalStrategies currently excludes FIRING and KCHOICES routing strategies
            // These may require special handling or are not supported in OutputStrategy context
            throw new RuntimeException("Unsupported Routing Strategy: " + routingStrategy);
        }
    }

    public OutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy) {
        this(jobClass, routingStrategy, null, 1);
    }

    public Node getDestination() {
        return this.destination;
    }

    public void setDestination(Node destination) {
        this.destination = destination;
    }

    public JobClass getJobClass() {
        return this.jobClass;
    }

    public double getProbability() {
        return this.probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    public RoutingStrategy getRoutingStrategy() {
        return this.routingStrategy;
    }

    public void setRoutingStrategy(RoutingStrategy routingStrategy) {
        this.routingStrategy = routingStrategy;
    }
}
