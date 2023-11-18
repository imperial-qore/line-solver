package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.sections.StatelessClassSwitcher;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

import java.util.*;

public class ClassSwitchArrivalEvent extends ArrivalEvent {
    public final Map<JobClass, Double> transitions;
    public ClassSwitchArrivalEvent(Node node, JobClass jobClass, StatelessClassSwitcher statelessClassSwitcher) {
        super(node, jobClass);

        this.transitions = new HashMap<>();
        for (JobClass jc : statelessClassSwitcher.getJobClasses()) {
            double transitionProbability = statelessClassSwitcher.applyCsFun(jobClass.getIndex() - 1, jc.getIndex() - 1);
            this.transitions.put(jc, transitionProbability);
        }
    }

    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        for (JobClass jobClassIter : this.transitions.keySet()) {
            ArrayList<Pair<SynchedEvent,Double>> outputEvents = this.node.getOutputEvents(jobClassIter, random);
            for (Pair<SynchedEvent,Double> outputEventDoublePair: outputEvents) {
                SSAStateMatrix newNetworkState = new SSAStateMatrix(networkState);
                outputEventDoublePair.getLeft().getNextState(newNetworkState, stateSpace,queue, stateSet);
            }        }
        return true;
    }

    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet, Pair<SynchedEvent, Double> outputEventDoublePair) {
        for (JobClass jobClassIter : this.transitions.keySet()) {
            ArrayList<Pair<SynchedEvent,Double>> outputEvents = this.node.getOutputEvents(jobClassIter, random);
            for (Pair<SynchedEvent,Double> outputEventDoublePair2: outputEvents) {
                SSAStateMatrix newNetworkState = new SSAStateMatrix(networkState);
                Pair<SynchedEvent, Double> newOutPutEventPair = new Pair<>(outputEventDoublePair2.getLeft(), outputEventDoublePair.getRight() * outputEventDoublePair2.getRight());
                newOutPutEventPair.getLeft().getNextEventState(newNetworkState, eventSpace, event, queue, copy, eventSet, newOutPutEventPair);
            }
        }
        return true;
    }

}
