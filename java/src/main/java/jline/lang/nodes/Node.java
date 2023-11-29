package jline.lang.nodes;

import jline.lang.*;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.sections.InputSection;
import jline.lang.sections.OutputSection;
import jline.lang.sections.ServiceSection;
import jline.solvers.ssa.events.ArrivalEvent;
import jline.solvers.ssa.events.NodeArrivalEvent;
import jline.solvers.ssa.events.SynchedEvent;
import jline.util.Pair;

import java.io.Serializable;
import java.util.*;

/**
 * Superclass for a node element within a Network model
 */
public class Node extends NetworkElement implements Serializable {
    public Network model;
    protected InputSection input;
    protected OutputSection output;
    protected ServiceSection server;
    protected DropStrategy dropStrategy;

    protected Map<JobClass, ArrivalEvent> arrivalEvents;

    protected int statefulIdx;
    protected int nodeIndex;
    protected int stationIdx;

    private final NodeAttribute attribute;

    public Node(String nodeName) {
        super(nodeName);
        this.arrivalEvents = new HashMap<JobClass, ArrivalEvent>();

        this.output = new OutputSection("Generic Output");
        this.input = new InputSection("Generic Input");
        this.dropStrategy = DropStrategy.Drop;
        this.statefulIdx = -1;
        this.nodeIndex = -1;
        this.stationIdx = -1;
        this.attribute = new NodeAttribute();
    }

    public void setModel(Network model) {
        this.model = model;
    }

    public Network getModel() {
        return this.model;
    }

    public void setScheduling(JobClass jobClass, SchedStrategy schedStrategy){
        this.input.setInputJobClasses(jobClass.getIndex(), schedStrategy);
    }

    public void setRouting(JobClass jobClass, RoutingStrategy routingStrategy) {
        this.output.setOutputStrategy(jobClass, routingStrategy);
    }

    public void setRouting(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        this.output.setOutputStrategy(jobClass, routingStrategy, destination, probability);
    }

    public void setProbRouting(JobClass jobClass, Node destination, double probability) {
        this.output.setOutputStrategy(jobClass, RoutingStrategy.PROB, destination, probability);
    }

    public void resetRouting() {
        this.output.resetRouting();
    }

    public RoutingStrategy getRoutingStrategy(JobClass jobClass) {
        for (OutputStrategy outputStrategy : this.output.getOutputStrategies()) {
            if (outputStrategy.getDestination() != null) {
                continue;
            }

            if (outputStrategy.getJobClass() == jobClass) {
                return outputStrategy.getRoutingStrategy();
            }
        }

        return RoutingStrategy.RAND;
    }

    public void printSummary() {
        System.out.format("jline.Node: %s\n", this.getName());
        this.output.printSummary();
    }

    public double getClassCap(JobClass jobClass) {
        return Double.POSITIVE_INFINITY;
    }

    public double getCap() { return Double.POSITIVE_INFINITY; }

    public ArrivalEvent getArrivalEvent(JobClass jobClass) {
        if (!this.arrivalEvents.containsKey(jobClass)) {
            this.arrivalEvents.put(jobClass, new NodeArrivalEvent(this, jobClass));
        }
        return this.arrivalEvents.get(jobClass);
    }

    public SynchedEvent getOutputEvent(JobClass jobClass, Random random) {
        return this.output.getOutputEvent(jobClass, random);
    }

    public ArrayList<Pair<SynchedEvent,Double>>  getOutputEvents(JobClass jobClass, Random random) {
        return this.output.getOutputEvents(jobClass, random);
    }


    public List<OutputStrategy> getOutputStrategies() {
        return this.output.getOutputStrategies();
    }

    public DropStrategy getDropStrategy() {
        return this.dropStrategy;
    }

    public boolean isRefstat() { return false; }

    public int getStatefulIdx() {
        this.statefulIdx = this.model.getStatefulNodeIndex(this);

        return this.statefulIdx;
    }

    public boolean isStateful() {
        return this.statefulIdx != -1;
    }

    public void setNodeIdx(int index) {
        this.nodeIndex = index;
    }

    public int getNodeIdx() {
    	if (this.nodeIndex == -1) {
    		this.nodeIndex = this.model.getNodeIndex(this);
    	}

    	return this.nodeIndex;
    }

    public void setStationIdx(int index) {
        this.stationIdx = index;
    }

    public int getStationIdx() {
    	if (this.stationIdx == -1) {
    		this.stationIdx = this.model.getStationIndex(this);
    	}

    	return this.stationIdx;
    }

    public InputSection getInput() {
    	return this.input;
    }

    public OutputSection getOutput() {
    	return this.output;
    }
    
    public ServiceSection getServer() {
    	return this.server;
    }

    public NodeAttribute getAttribute() {
        return attribute;
    }

    public void reset() {
            // Reset internal data structures when the network model is
            // reset
            /* no-op */
    }
}
