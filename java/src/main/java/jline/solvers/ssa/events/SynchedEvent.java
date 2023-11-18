package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.sections.OutputSection;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

public class SynchedEvent extends Event {
    public OutputSection outputSection;
    public Node node;
    public JobClass jobClass;
    public boolean isClassSwitched;
    public int jobClassIdx;
    public boolean isDummy = false;

    public SynchedEvent(OutputSection outputSection, Node node, JobClass jobClass) {
        super();
        this.jobClass = jobClass;
        this.node = node;
        this.isClassSwitched = false;
        this.jobClassIdx = jobClass.getJobClassIdx();
    }

    public SynchedEvent(OutputSection outputSection, Node node, JobClass jobClass, boolean isClassSwitched) {
        this(outputSection, node, jobClass);
        this.isClassSwitched = isClassSwitched;
    }

    //WARNING: THIS IS ONLY USED AS A DUMMY EVENT FOR SELF-LOOPING EVENTS IN CTMC
    public SynchedEvent(int jobClassIdx) {
        this.jobClassIdx = jobClassIdx;
        this.isDummy = true;
    }

    @Override
    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        return this.node.getArrivalEvent(this.jobClass).updateStateSpace(networkState, random, stateSpace, queue, stateSet);
    }

    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet, Pair<SynchedEvent, Double> outputEventDoublePair) {
        return this.node.getArrivalEvent(this.jobClass).updateEventSpace(networkState, random, eventSpace, event, queue, copy, eventSet, outputEventDoublePair) ;
    }


    public boolean isClassSwitched() {
        return this.isClassSwitched;
    }

    public int getClassIdx() {
        return this.jobClassIdx;
    }

    public OutputSection getOutputSection() {
        return this.outputSection;
    }

    public SSAStateMatrix getNextState(SSAStateMatrix startingState, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {

        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateStateSpace(endingState, new Random(), stateSpace,queue, stateSet)){
            return endingState;
        }

        return null;

    }

    public SSAStateMatrix getNextEventState(SSAStateMatrix startingState, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet, Pair<SynchedEvent, Double> outputEventDoublePair) {

        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateEventSpace(endingState, new Random(), eventSpace,event,queue,copy, eventSet, outputEventDoublePair)){
            return endingState;
        }

        return null;
    }

    @Override
    public Node getNode() {
        return node;
    }

    public boolean isDummy() {
        return this.isDummy;
    }
}
