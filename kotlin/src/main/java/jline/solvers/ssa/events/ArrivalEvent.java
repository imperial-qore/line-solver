package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.lang.nodes.StatefulNode;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.solvers.ssa.state.PhaseList;
import jline.lang.distributions.CumulativeDistribution;
import jline.util.Pair;

import java.io.Serializable;
import java.util.*;

public class ArrivalEvent extends Event implements NodeEvent, Serializable {
    public final int statefulIndex;
    public final int classIndex;
    public final boolean useBuffer;
    public final boolean isStateful;
    public JobClass jobClass;
    public Node node;

    public ArrivalEvent(Node node, JobClass jobClass) {
        super();

        this.jobClass = jobClass;
        this.node = node;

        if (node instanceof StatefulNode) {
            this.statefulIndex = ((StatefulNode)node).getStatefulIndex();
            this.isStateful = true;
        } else {
            this.statefulIndex = -1;
            this.isStateful = false;
        }
        this.classIndex = this.node.getModel().getJobClassIndex(this.jobClass);

        this.useBuffer = true;
    }


    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        if (this.isStateful) {
            if(networkState.getState(this.statefulIndex, this.classIndex) >= networkState.getCapacity(this.statefulIndex, this.classIndex)) {
                return false;
            }
            PhaseList phaseList = networkState.getPhaseList(this.statefulIndex);
            CumulativeDistribution<Integer> startingPhaseProbabilities =  phaseList.getStartingPhaseProbabilities(this.classIndex);
            ArrayList<Pair<Double, Integer>> possibleStarts = startingPhaseProbabilities.getPossibleEventProbability();
            for(Pair<Double, Integer> start : possibleStarts) {
                if(start.getLeft() > 0) {
                    SSAStateMatrix copy = new SSAStateMatrix(networkState);
                    copy.stateArrivalAtPosition(this.statefulIndex, this.classIndex, start.getRight());
                    if (!copy.exceedsCutoff() && !stateSet.contains(copy)) {
                        stateSet.add(copy);
                        stateSpace.add(copy);
                        queue.add(copy);
                    }
                }
            }
            return true;
        } else if (!(this.node instanceof Sink)){
            throw new RuntimeException(String.format("ArrivalEvent at %s not supported!", this.node.getName()));
        }

        return true;
    }

    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet, Pair<SynchedEvent, Double> outputEventDoublePair) {
        if (this.isStateful) {
            if(networkState.getState(this.statefulIndex, this.classIndex) >= networkState.getCapacity(this.statefulIndex, this.classIndex)) {
                return false;
            }
            PhaseList phaseList = networkState.getPhaseList(this.statefulIndex);
            CumulativeDistribution<Integer> startingPhaseProbabilities =  phaseList.getStartingPhaseProbabilities(this.classIndex);
            ArrayList<Pair<Double, Integer>> possibleStarts = startingPhaseProbabilities.getPossibleEventProbability();
            for(Pair<Double, Integer> start : possibleStarts) {
                if(start.getLeft() > 0) {
                    SSAStateMatrix copy2 = new SSAStateMatrix(networkState);
                    copy2.stateArrivalAtPosition(this.statefulIndex, this.classIndex, start.getRight());
                    Pair<SynchedEvent, Double> newOutputEventPair = new Pair<>(outputEventDoublePair.getLeft(), outputEventDoublePair.getRight() * start.getLeft());
                    EventData eventData = new EventData(event, newOutputEventPair,copy,copy2);
                    if (!copy2.exceedsCutoff() && !eventSet.contains(eventData)) {
                        eventSet.add(eventData);
                        eventSpace.add(eventData);
                        queue.add(copy2);
                    }
                }
            }
            return true;
        } else if (!(this.node instanceof Sink)){
            throw new RuntimeException(String.format("ArrivalEvent at %s not supported!", this.node.getName()));
        }
        EventData eventData = new EventData(event, outputEventDoublePair,copy, networkState);
        if (!networkState.exceedsCutoff() && !eventSet.contains(eventData)) {
            eventSet.add(eventData);
            eventSpace.add(eventData);
            queue.add(networkState);
        }
        return true;
    }

    @Override
    public void printSummary() {
        System.out.format("Arrival event for %s at %s\n", this.jobClass.getName(), this.node.getName());
    }

    public Node getNode() {
        return this.node;
    }

    public JobClass getJobClass() { return this.jobClass; }

    public int getNodeStatefulIdx() {
        return this.statefulIndex;
    }
    public int getClassIdx() {
        return this.classIndex;
    }


    public boolean isStateful() {
        return this.statefulIndex != -1;
    }
}
