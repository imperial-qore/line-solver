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
     * Routing strategies an OutputStrategy may carry. RL belongs here: it is a
     * state-dependent routing strategy on the same footing as JSQ and SQ
     * (see Network.refreshRouting, Network.sub_rl, SnHasSDRouting and the
     * RoutingStrategy_RL feature), and Node.setRLRouting installs it through
     * this class. Omitting it made setRLRouting throw, so RL routing could not
     * be configured at all. FIRING is excluded: a Transition's outgoing arcs are
     * described by its firing outcomes, not by an output strategy.
     */
    public static List<RoutingStrategy> legalStrategies = Arrays.asList(RoutingStrategy.DISABLED, RoutingStrategy.PROB, RoutingStrategy.RAND, RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.SQ, RoutingStrategy.RL);
    private final JobClass jobClass;
    private RoutingStrategy routingStrategy;
    private double probability;
    private Node destination;
    /** SQ: d, number of randomly sampled destinations to compare (>= 1). */
    private int sqD = 2;
    /**
     * RL: value function. Flat tabular table when rlStateSize=0, coefficient row
     * vector when rlStateSize&gt;0. Held here, alongside the SQ parameters,
     * because the NetworkStruct is a derived cache that refreshStruct rebuilds:
     * parameters kept only there do not survive a refresh.
     */
    private Matrix rlValueFunction;
    /** RL: per-axis sizes of the tabular value function; used when rlStateSize=0. */
    private int[] rlValueFunctionShape;
    /** RL: node indices that consult the value function; others fall back to JSQ. */
    private int[] rlNodesNeedAction;
    /** RL: 0 = tabular, &gt;0 = linear approximation, &lt;0 = JSQ fallback. */
    private int rlStateSize = -1;

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

    /**
     * Returns the RL value function, or null when RL routing is not configured.
     *
     * @return the value function
     */
    public Matrix getRlValueFunction() {
        return this.rlValueFunction;
    }

    /**
     * Returns the per-axis shape of the tabular RL value function.
     *
     * @return the shape, or null when not configured
     */
    public int[] getRlValueFunctionShape() {
        return this.rlValueFunctionShape;
    }

    /**
     * Returns the node indices that consult the RL value function.
     *
     * @return the node indices, or null when not configured
     */
    public int[] getRlNodesNeedAction() {
        return this.rlNodesNeedAction;
    }

    /**
     * Returns the RL state-size hint.
     *
     * @return 0 for tabular, &gt;0 for linear approximation, &lt;0 for JSQ fallback
     */
    public int getRlStateSize() {
        return this.rlStateSize;
    }

    /**
     * Stores the RL routing parameters on this strategy.
     *
     * @param valueFunction   the value function
     * @param vfShape         per-axis sizes of the tabular value function
     * @param nodesNeedAction node indices that consult the value function
     * @param stateSize       0 = tabular, &gt;0 = linear approximation, &lt;0 = JSQ fallback
     */
    public void setRlParams(Matrix valueFunction, int[] vfShape, int[] nodesNeedAction,
                            int stateSize) {
        this.rlValueFunction = valueFunction;
        this.rlValueFunctionShape = vfShape;
        this.rlNodesNeedAction = nodesNeedAction;
        this.rlStateSize = stateSize;
    }

}
