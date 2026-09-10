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
public class Node extends NetworkElement implements Serializable, Cloneable {
    private final NodeAttribute attribute;
    public Network model;
    protected InputSection input;
    protected OutputSection output;
    protected ServiceSection server;
    protected DropStrategy dropStrategy;
    protected int statefulIdx;
    protected int nodeIndex;
    public int stationIdx;
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
     * Removes all per-class configuration referencing the given job class from
     * this node. The base implementation drops the class's output (routing)
     * strategies; subclasses extend it to drop service, capacity, arrival and
     * class-switching configuration. Called by {@link Network#removeClass}.
     *
     * @param jobClass the job class being removed from the model
     */
    public void removeJobClass(JobClass jobClass) {
        if (this.output != null) {
            this.output.getOutputStrategies().removeIf(
                os -> os.getJobClass() == jobClass);
        }
    }

    /**
     * Returns the routing strategy configured for a specific job class.
     *
     * @param jobClass the job class to query
     * @return the routing strategy, or RAND if none specified
     */
    public RoutingStrategy getRoutingStrategy(JobClass jobClass) {
        // First pass: look for global (no destination) strategy entries
        for (OutputStrategy outputStrategy : this.output.getOutputStrategies()) {
            if (outputStrategy.getDestination() != null) {
                continue;
            }

            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (outputStrategy.getJobClass().getIndex() == jobClass.getIndex()) {
                return outputStrategy.getRoutingStrategy();
            }
        }

        // Second pass: check destination-based entries (e.g., WRROBIN stores strategy per-destination)
        for (OutputStrategy outputStrategy : this.output.getOutputStrategies()) {
            if (outputStrategy.getDestination() == null) {
                continue;
            }
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
            // jline.io.InputOutput.line_error(jline.io.InputOutput.mfilename(this),
            //     "Cannot set routing strategy on Cache node. Cache nodes must route to a Router node using setProbRouting instead.");
            return;
        }
        this.output.setOutputStrategy(jobClass, routingStrategy);
    }

    /** Declared Krzesinski state-dependent routing structures, one per job class. */
    protected java.util.Map<JobClass, jline.lang.StateDepRouting> sdrDeclarations =
            new java.util.HashMap<JobClass, jline.lang.StateDepRouting>();

    /**
     * Declares this node the entry center e of a subnetwork Q(V,V) served by the
     * product-form state-dependent routing of Krzesinski (1987), "Multiclass
     * Queueing Networks with State-Dependent Routing", Performance Evaluation
     * 7(2):125-143.
     *
     * <p>The departure center may be this node itself in a central server model.
     * Branches follow the paper's own indexing: index 0 of the list is the
     * unused complement M-V and must be null or empty, and entry b lists the
     * nodes of branch b with its entry center first and its departure center
     * last. A single-center branch is a one-element list. level[b] is the index
     * t of the subnetwork with B_b in V_t - V_{t+1}, and level[0] is ignored.</p>
     *
     * <p>Negative C_t and positive d_tb make the routing prefer the least
     * congested branches and impose the population bounds m_b &lt;= d_tb/(-C_t)
     * and v_t &lt;= D_tt/(-C_t).</p>
     *
     * @param jobClass  the job class routed by this subnetwork
     * @param departure the departure center d of Q(V,V)
     * @param branches  branch node lists, index 0 unused
     * @param level     branch levels, index 0 unused
     * @param C         coefficients C_t, length T
     * @param d         coefficients d_tb, T rows by B columns
     */
    public void setStateDepRouting(JobClass jobClass, Node departure,
                                   java.util.List<java.util.List<Node>> branches,
                                   int[] level, double[] C, double[][] d) {
        if (branches == null || branches.size() < 2
                || (branches.get(0) != null && !branches.get(0).isEmpty())) {
            throw new RuntimeException("The branch list must start with an empty entry: "
                    + "branch index 1 denotes the complement M-V.");
        }
        int B = branches.size();
        int T = C.length;
        if (level.length != B) {
            throw new RuntimeException("level must have one entry per branch index, including the unused index 0.");
        }
        if (d.length < T) {
            throw new RuntimeException("d must be at least " + T + "x" + B + ".");
        }
        jline.lang.StateDepRouting sdr = new jline.lang.StateDepRouting();
        sdr.entryNode = this;
        sdr.departureNode = departure;
        sdr.branchNodes = new java.util.ArrayList<java.util.List<Node>>(branches);
        sdr.level = level.clone();
        sdr.C = C.clone();
        sdr.d = new double[d.length][];
        for (int t = 0; t < d.length; t++) {
            sdr.d[t] = d[t].clone();
        }
        for (int b = 1; b < B; b++) {
            if (branches.get(b) == null || branches.get(b).isEmpty()) {
                throw new RuntimeException("Branch " + (b + 1) + " is empty.");
            }
        }
        this.sdrDeclarations.put(jobClass, sdr);
        this.output.setOutputStrategy(jobClass, RoutingStrategy.SDR);
    }

    /**
     * Returns the state-dependent routing structure declared for a job class.
     *
     * @param jobClass the job class
     * @return the declared structure, or null when the class is not routed by SDR
     */
    public jline.lang.StateDepRouting getStateDepRouting(JobClass jobClass) {
        return this.sdrDeclarations.get(jobClass);
    }

    /**
     * Configures power-of-K-choices routing for a job class.
     *
     * @param jobClass    the job class
     * @param d           number of randomly sampled destinations (must be >= 1)
     */
    public void setSQRouting(JobClass jobClass, int d) {
        if (d < 1) {
            throw new IllegalArgumentException("d must be >= 1");
        }
        this.output.setSQParams(jobClass, d);
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
            // jline.io.InputOutput.line_error(jline.io.InputOutput.mfilename(this),
            //     "Cannot set routing strategy on Cache node. Cache nodes must route to a Router node using setProbRouting instead.");
            return;
        }
        this.output.setOutputStrategy(jobClass, routingStrategy, destination, probability);
    }

    /**
     * Sets a routing strategy taking a single parameter, twin of the MATLAB
     * {@code setRouting(class, strategy, par1)} call.
     *
     * SQ reads the parameter as the number d of sampled destinations; every
     * other strategy takes none, so the parameter is ignored and the strategy is
     * declared as in {@link #setRouting(JobClass, RoutingStrategy)}.
     *
     * @param jobClass        the job class to configure routing for
     * @param routingStrategy the routing strategy to use
     * @param param           the strategy parameter
     */
    public void setRouting(JobClass jobClass, RoutingStrategy routingStrategy, Object param) {
        if (routingStrategy == RoutingStrategy.SQ) {
            if (!(param instanceof Number)) {
                throw new IllegalArgumentException("SQ parameter d must be a positive integer.");
            }
            double d = ((Number) param).doubleValue();
            if (d < 1 || Math.abs(d - Math.round(d)) > jline.GlobalConstants.CoarseTol) {
                throw new IllegalArgumentException("SQ parameter d must be a positive integer.");
            }
            this.setSQRouting(jobClass, (int) Math.round(d));
            return;
        }
        this.setRouting(jobClass, routingStrategy);
    }

    /**
     * Checks whether this node switches the class of the jobs traversing it.
     *
     * @return true when the service section is a class switcher
     */
    public boolean hasClassSwitching() {
        return this.server instanceof jline.lang.sections.ClassSwitcher;
    }

    /**
     * Checks whether this node is a station, i.e. a node holding jobs in a queue
     * or in service.
     *
     * @return true if this node is a Station
     */
    public boolean isStation() {
        return this instanceof Station;
    }

    /**
     * Links this node to another one in the model it belongs to.
     *
     * @param nodeTo the destination node
     * @return this node, so calls can be chained
     */
    public Node link(Node nodeTo) {
        this.model.addLink(this, nodeTo);
        return this;
    }

    /**
     * Prints the one-line node summary, twin of MATLAB {@code Node.summary}.
     */
    public void summary() {
        System.out.format("\nNode: %s\n", this.getName());
    }

    /**
     * Returns a copy of this node, twin of the MATLAB {@code Node.copyElement}:
     * every field is copied shallowly, then the three sections are copied so the
     * clone can be re-wired without disturbing the original. The model handle
     * stays shared, as it does in MATLAB.
     *
     * @return a copy of this node
     */
    protected Node copyElement() {
        Node clone;
        try {
            clone = (Node) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException("Failed to copy node " + this.getName(), e);
        }
        if (this.input != null) {
            clone.input = (InputSection) this.input.copyElement();
        }
        if (this.server != null) {
            clone.server = (ServiceSection) this.server.copyElement();
        }
        if (this.output != null) {
            clone.output = (OutputSection) this.output.copyElement();
        }
        if (this.state != null) {
            clone.state = this.state.copy();
        }
        return clone;
    }

    /**
     * Returns a copy of this node. Overrides the serialization-based
     * {@link jline.lang.Copyable#copy()} because that one would also clone the
     * Network reached through {@link #model}; MATLAB shares the model handle
     * instead, and {@link #copyElement()} reproduces that.
     *
     * @return a copy of this node
     */
    @Override
    @SuppressWarnings("unchecked")
    public <T extends jline.lang.Copyable> T copy() {
        return (T) this.copyElement();
    }

}
