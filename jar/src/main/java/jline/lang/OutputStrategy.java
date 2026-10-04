/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Class modelling the output section of a Node
 */
public class OutputStrategy implements Serializable {
    /**
     * Routing strategies an OutputStrategy may carry. FIRING is excluded: a
     * Transition's outgoing arcs are described by its firing outcomes, not by
     * an output strategy.
     */
    public static List<RoutingStrategy> legalStrategies = Arrays.asList(RoutingStrategy.DISABLED, RoutingStrategy.PROB, RoutingStrategy.RAND, RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.SQ, RoutingStrategy.SDR);
    private final JobClass jobClass;
    private RoutingStrategy routingStrategy;
    private double probability;
    private Node destination;
    /** SQ: d, number of randomly sampled destinations to compare (>= 1). */
    private int sqD = 2;

    public OutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        this.jobClass = jobClass;
        this.routingStrategy = routingStrategy;
        this.destination = destination;
        this.probability = probability;

        if (!legalStrategies.contains(routingStrategy)) {
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

    public int getSqD() {
        return this.sqD;
    }

    public void setSqD(int d) {
        this.sqD = d;
    }


}
