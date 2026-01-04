/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import static jline.GlobalConstants.Inf;

import jline.lang.*;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.sections.InputSection;
import jline.lang.sections.OutputSection;
import jline.lang.sections.Section;
import jline.lang.sections.ServiceSection;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Superclass for a node element within a Network model
 */
public class Node extends NetworkElement implements Serializable {
    private final NodeAttribute attribute;
    public Network model;
    protected InputSection input;
    protected OutputSection output;
    protected ServiceSection server;
    protected DropStrategy dropStrategy;
    protected int statefulIdx;
    protected int nodeIndex;
    protected int stationIdx;
    protected Matrix state;

    /**
     * Creates a new node with the specified name.
     * Initializes default routing and service configurations.
     *
     * @param nodeName the name for this node
     */
    public Node(String nodeName) {
        super(nodeName);

        this.output = new OutputSection("Generic Output");
        this.input = new InputSection("Generic Input");
        this.dropStrategy = DropStrategy.Drop;
        this.statefulIdx = -1;
        this.nodeIndex = -1;
        this.stationIdx = -1;
        this.attribute = new NodeAttribute();
    }

    /**
     * Gets the attribute object containing additional metadata for this node.
     * 
     * @return the node attribute object
     */
    public NodeAttribute getAttribute() {
        return attribute;
    }

    /**
     * Returns the total capacity limit for this node.
     * Default implementation returns infinite capacity.
     *
     * @return the total capacity limit
     */
    public double getCap() {
        return Inf;
    }

    /**
     * Returns the capacity limit for a specific job class at this node.
     * Default implementation returns infinite capacity.
     *
     * @param jobClass the job class to query
     * @return the capacity limit for the job class
     */
    public double getClassCap(JobClass jobClass) {
        return Inf;
    }

    /**
     * Gets the drop strategy used by this node when capacity limits are exceeded.
     * 
     * @return the drop strategy
     */
    public DropStrategy getDropStrategy() {
        return this.dropStrategy;
    }

    /**
     * Gets the input section that handles incoming jobs to this node.
     * 
     * @return the input section
     */
    public InputSection getInput() {
        return this.input;
    }

    /**
     * Returns the network model containing this node.
     *
     * @return the parent network model
     */
    public Network getModel() {
        return this.model;
    }

    /**
     * Sets the network model containing this node.
     *
     * @param model the parent network model
     */
    public void setModel(Network model) {
        this.model = model;
        this.nodeIndex = -1; // Reset cached index when model changes
    }

    /**
     * Gets the index of this node in the network's node collection.
     * The index is lazily computed if not already set.
     * 
     * @return the node index
     */
    public int getNodeIndex() {
        if (this.nodeIndex == -1) {
            this.nodeIndex = this.model.getNodeIndex(this);
        }

        return this.nodeIndex;
    }

    /**
     * Sets the index of this node in the network's node collection.
     * 
     * @param index the node index to set
     */
    public void setNodeIdx(int index) {
        this.nodeIndex = index;
    }

    /**
     * Gets the output section that handles job routing from this node.
     * 
     * @return the output section
     */
    public OutputSection getOutput() {
        return this.output;
    }

    /**
     * Returns the list of output strategies configured for this node.
     *
     * @return list of output routing strategies
     */
    public List<OutputStrategy> getOutputStrategies() {
        return this.output.getOutputStrategies();
    }

    /**
     * Returns the routing strategy configured for a specific job class.
     *
     * @param jobClass the job class to query
     * @return the routing strategy, or RAND if none specified
     */
    public RoutingStrategy getRoutingStrategy(JobClass jobClass) {
        for (OutputStrategy outputStrategy : this.output.getOutputStrategies()) {
            if (outputStrategy.getDestination() != null) {
                continue;
            }

            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (outputStrategy.getJobClass().getIndex() == jobClass.getIndex()) {
                return outputStrategy.getRoutingStrategy();
            }
        }

        return RoutingStrategy.RAND;
    }

    /**
     * Gets all sections (input, server, output) that compose this node.
     * 
     * @return list containing all node sections
     */
    public List<Section> getSections() {
        List<Section> ret = new ArrayList<>();
        ret.add(this.input);
        ret.add(this.server);
        ret.add(this.output);
        return ret;
    }

    /**
     * Gets the service section that handles job processing at this node.
     * 
     * @return the service section
     */
    public ServiceSection getServer() {
        return this.server;
    }

    /**
     * Gets the index of this node in the network's stateful node collection.
     * Returns -1 if this node is not stateful.
     * 
     * @return the stateful node index, or -1 if not stateful
     */
    public int getStatefulIdx() {
        this.statefulIdx = this.model.getStatefulNodeIndex(this);

        return this.statefulIdx;
    }

    /**
     * Gets the index of this node in the network's station collection.
     * The index is lazily computed if not already set.
     * 
     * @return the station index
     */
    public int getStationIdx() {
        if (this.stationIdx == -1) {
            this.stationIdx = this.model.getStationIndex(this);
        }

        return this.stationIdx;
    }

    /**
     * Sets the index of this node in the network's station collection.
     * 
     * @param index the station index to set
     */
    public void setStationIdx(int index) {
        this.stationIdx = index;
    }

    /**
     * Checks if this node is a reference station.
     * Default implementation returns false; override in subclasses as needed.
     * 
     * @return true if this is a reference station, false otherwise
     */
    public boolean isReferenceStation() {
        return false;
    }

    /**
     * Checks if this node maintains state between job visits.
     * 
     * @return true if this node is stateful, false otherwise
     */
    public boolean isStateful() {
        return this.model.getStatefulNodeIndex(this) != -1; // call needed to ensure statefulIdx is generated
    }
    
    /**
     * Sets the state of this node.
     * 
     * @param state the state matrix to set
     */
    public void setState(Matrix state) {
        this.state = state;
    }
    
    /**
     * Gets the current state of this node.
     * 
     * @return the current state matrix
     */
    public Matrix getState() {
        return this.state;
    }

    /**
     * Prints a summary of this node's configuration to the console.
     */
    public void printSummary() {
        System.out.format("jline.Node: %s\n", this.getName());
        this.output.printSummary();
    }

    /**
     * Resets the internal state of this node.
     * Called when the network model is reset.
     * Default implementation is no-op; override in subclasses as needed.
     */
    public void reset() {
        // Reset internal data structures when the network model is
        // reset
        /* no-op */
    }

    /**
     * Resets all routing configurations for this node.
     */
    public void resetRouting() {
        this.output.resetRouting();
    }

    /**
     * Sets probabilistic routing for a job class to a specific destination.
     *
     * @param jobClass    the job class to configure routing for
     * @param destination the destination node
     * @param probability the routing probability (0.0 to 1.0)
     */
    public void setProbRouting(JobClass jobClass, Node destination, double probability) {
        this.output.setOutputStrategy(jobClass, RoutingStrategy.PROB, destination, probability);
    }

    /**
     * Sets the routing strategy for a specific job class.
     *
     * @param jobClass        the job class to configure routing for
     * @param routingStrategy the routing strategy to use
     */
    public void setRouting(JobClass jobClass, RoutingStrategy routingStrategy) {
        // Cache nodes can only have PROB routing strategy set directly
        if (this instanceof Cache && routingStrategy != RoutingStrategy.PROB) {
            // jline.io.InputOutputKt.line_error(jline.io.InputOutputKt.mfilename(this),
            //     "Cannot set routing strategy on Cache node. Cache nodes must route to a Router node using setProbRouting instead.");
            return;
        }
        this.output.setOutputStrategy(jobClass, routingStrategy);
    }

    /**
     * Sets the routing strategy with destination and probability for a job class.
     *
     * @param jobClass        the job class to configure routing for
     * @param routingStrategy the routing strategy to use
     * @param destination     the destination node
     * @param probability     the routing probability to this destination
     */
    public void setRouting(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        // Cache nodes can only have PROB routing strategy set directly
        if (this instanceof Cache && routingStrategy != RoutingStrategy.PROB) {
            // jline.io.InputOutputKt.line_error(jline.io.InputOutputKt.mfilename(this),
            //     "Cannot set routing strategy on Cache node. Cache nodes must route to a Router node using setProbRouting instead.");
            return;
        }
        this.output.setOutputStrategy(jobClass, routingStrategy, destination, probability);
    }

    /*
     * ===================================================================================
     * MISSING METHODS FROM MATLAB NODE IMPLEMENTATION - NOT YET MIGRATED
     * ===================================================================================
     *
     * Based on analysis of /matlab/src/lang/nodes/Node.m
     */

    // =================== ENHANCED ROUTING METHODS ===================
    // public void setRouting(JobClass jobClass, RoutingStrategy strategy, Object param)  // Single parameter routing
    // public void setRouting(JobClass jobClass, RoutingStrategy strategy, Node destination, double probability)  // Enhanced routing with destination

    // =================== NODE TYPE CHECKING METHODS ===================
    // public boolean hasClassSwitching()  // Check if node has class switching capability
    // public boolean isStateful()         // Check if node is stateful
    // public boolean isStation()          // Check if node is a station

    // =================== ENHANCED NODE MANAGEMENT ===================
    // public Node link(Node nodeTo)       // Link this node to another node
    // public void summary()               // Enhanced summary with formatted output

    // =================== DEEP COPY SUPPORT ===================
    // public Node copyElement()           // Create deep copy of node with all sections
    // public Node copy()                  // Public interface for copying

}
