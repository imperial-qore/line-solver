package jline.solvers.ssa.events;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

import jline.lang.nodes.Node;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;

public class Event implements Serializable {
    public Node node;
    public Event() {
    }

    public double getRate(SSAStateMatrix networkState) {
        return Double.NaN;
    }

    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {

        return true;
    }

    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {

        return true;
    }


//    public int stateUpdateN(int n, SSAStateMatrix networkState, Random random, SSARunner SSARunner) {
//        /*
//            stateUpdateN -
//                Attempt to apply N repetitions of an event to the stateMatrix
//
//            Returns: (int) - number of repetitions left unapplied
//         */
//        int rem = n;
//        for (int i = 0; i < n; i++) {
//            if (this.stateUpdate(networkState, random, SSARunner)) {
//                rem--;
//            }
//        }
//
//        return rem;
//    }

    public void printSummary() {
        System.out.format("Generic event\n");
    }

    public int getMaxRepetitions(SSAStateMatrix networkState) {
        return Integer.MAX_VALUE;
    }

    public SSAStateMatrix getNextState(SSAStateMatrix startingState, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        //TO_DO why are there overrides if they don't change anything??? remove them
        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateStateSpace(endingState, new Random(), stateSpace,queue, stateSet)){
            return endingState;
        }
        return null;

    }

    public SSAStateMatrix getNextEventState(SSAStateMatrix startingState, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {

        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateEventSpace(endingState, new Random(), eventSpace,event,queue,copy, eventSet)){
            return endingState;
        }

        return null;
    }

    public Node getNode() {
        return node;
    }
}
