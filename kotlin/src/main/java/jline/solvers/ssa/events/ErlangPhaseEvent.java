package jline.solvers.ssa.events;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Erlang;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.lang.nodes.StatefulNode;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.solvers.ssa.state.PhaseList;
import jline.util.Pair;

import java.util.*;

public class ErlangPhaseEvent extends PhaseEvent implements NodeEvent {
    public final int statefulIndex;
    public final int classIndex;
    public final SchedStrategy schedStrategy;
    public final boolean isSource;
    public final Erlang serviceProcess;
    public Node node;
    public final JobClass jobClass;
    public boolean isProcessorSharing;

    public final DepartureEvent departureEvent;

    public ErlangPhaseEvent(Node node, JobClass jobClass, DepartureEvent departureEvent) {
        super();
        this.node = node;
        this.jobClass = jobClass;

        if (node instanceof StatefulNode) {
            this.statefulIndex = ((StatefulNode)this.node).getStatefulIndex();
        } else {
            this.statefulIndex = -1;
        }
        this.classIndex = this.node.getModel().getJobClassIndex(this.jobClass);

        this.isSource = node instanceof Source;
        if (!(node instanceof HasSchedStrategy)) {
            throw new RuntimeException("Scheduling strategy required");
        }

        this.schedStrategy = ((HasSchedStrategy)node).getSchedStrategy();
        Distribution distServiceProcess = ((HasSchedStrategy)node).getServiceProcess(this.jobClass);
        if (!(distServiceProcess instanceof Erlang)) {
            throw new RuntimeException("Erlang distribution required");
        }

        this.serviceProcess = (Erlang)distServiceProcess;
        this.departureEvent = departureEvent;

        this.isProcessorSharing = this.schedStrategy == SchedStrategy.PS;
    }

    @Override
    public long getNPhases() {
        return this.serviceProcess.getNumberOfPhases();
    }

    @Override
    public double getRate(SSAStateMatrix networkState) {
        int activeServers = 1;

        if (this.isProcessorSharing) {
            double serviceRatio = (double) networkState.getState(this.statefulIndex, this.classIndex)/(double) networkState.totalStateAtNode(this.statefulIndex);
            serviceRatio *= networkState.psTotalCapacity(this.statefulIndex);
            return ((Double)this.serviceProcess.getParam(1).getValue())*serviceRatio;        }

        if (this.node instanceof StatefulNode) {
            activeServers = networkState.inProcess(this.statefulIndex, this.classIndex);
            if (this.node instanceof Source) {
                // NOTE: Pay active attention to this part
                activeServers = 1;//stateMatrix.getPhaseListSize(this)+1;
            } else if (activeServers == 0) {
                return Double.NaN;
            }
        }

        return ((Double)this.serviceProcess.getParam(1).getValue())*activeServers;
    }

    public double getDepartureRate(SSAStateMatrix networkState) {
        int activeServers = 1;
        if (this.isProcessorSharing) {
            double serviceRatio = (double) networkState.getState(this.statefulIndex, this.classIndex)/(double) networkState.totalStateAtNode(this.statefulIndex);
            serviceRatio *= networkState.psTotalCapacity(this.statefulIndex);
            return ((Double)this.serviceProcess.getParam(1).getValue())*serviceRatio* networkState.getInPhase(this.statefulIndex, this.classIndex, (int) this.serviceProcess.getNumberOfPhases() - 1);
        }

        if (this.node instanceof StatefulNode) {
            activeServers = networkState.inProcess(this.statefulIndex, this.classIndex);
            if (this.node instanceof Source) {
                // NOTE: Pay active attention to this part
                activeServers = 1;//stateMatrix.getPhaseListSize(this)+1;
            } else if (activeServers == 0) {
                return Double.NaN;
            }
        }

        return ((Double)this.serviceProcess.getParam(1).getValue())* networkState.getInPhase(this.statefulIndex, this.classIndex, (int) this.serviceProcess.getNumberOfPhases() - 1);
    }


    @Override
    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        boolean ok = false;
        PhaseList phaseList = networkState.getPhaseList(statefulIndex);
        for(int i = 0; i < phaseList.getNPhases(classIndex); i++) {
            SSAStateMatrix copy = new SSAStateMatrix(networkState);
            int currPhase = phaseList.getNInPhase(classIndex, i);
            if(currPhase > 0) {
                ok = true;
                if(i == phaseList.getNPhases(classIndex) - 1) {
                    copy.updatePhase(statefulIndex, classIndex, i, -1);
                    if(this.node instanceof Source) {
                        copy.updatePhase(statefulIndex, classIndex, -1, 0);
                    }
                    this.departureEvent.getNextState(copy, stateSpace, queue, stateSet);
                }
                else {
                    copy.updatePhase(statefulIndex, classIndex, i, i + 1);
                    if (!copy.exceedsCutoff() && !stateSet.contains(copy)) {
                        stateSet.add(copy);
                        stateSpace.add(copy);
                        queue.add(copy);
                    }
                }
            }
        }
        if(!ok && this.node instanceof Source) {
            SSAStateMatrix copy = new SSAStateMatrix(networkState);
            copy.updatePhase(statefulIndex, classIndex, -1, 0);
            stateSpace.remove(0);
            stateSet.remove(networkState);
            stateSet.add(copy);
            stateSpace.add(copy);
            queue.add(copy);
        }
        return true;
    }

    @Override
    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {
        PhaseList phaseList = networkState.getPhaseList(statefulIndex);
        boolean ok = false;
        for(int i = 0; i < phaseList.getNPhases(classIndex); i++) {
            SSAStateMatrix copy2 = new SSAStateMatrix(networkState);
            int currPhase = phaseList.getNInPhase(classIndex, i);
            if(currPhase > 0) {
                ok = true;
                if(i == phaseList.getNPhases(classIndex) - 1) {
                    copy2.updatePhase(statefulIndex, classIndex, i, -1);
                    if(this.node instanceof Source) {
                        copy2.updatePhase(statefulIndex, classIndex, -1, 0);
                    }
                    this.departureEvent.getNextEventState(copy2, eventSpace,event, queue,copy, eventSet);
                }
                else {
                    copy2.updatePhase(statefulIndex, classIndex, i, i + 1);
                    int phase = copy.findPhaseChange(copy2, this.statefulIndex, this.classIndex);
                    double rate =  copy.getInPhase(this.statefulIndex, classIndex, phase) * ((Double)this.serviceProcess.getParam(1).getValue());
                    SynchedEvent dummyEvent = new SynchedEvent(classIndex);
                    Pair<SynchedEvent, Double> pair = new Pair<>(dummyEvent, rate);
                    EventData eventData = new EventData(event, pair, copy, copy2);
                    if (!copy2.exceedsCutoff() && !eventSet.contains(eventData)) {
                        eventSet.add(eventData);
                        eventSpace.add(eventData);
                        queue.add(copy2);
                    }
                }
            }
        }
        if(!ok && this.node instanceof Source) {
            SSAStateMatrix copy2 = new SSAStateMatrix(networkState);
            copy2.updatePhase(statefulIndex, classIndex, -1, 0);
            EventData eventData = new EventData(null, null, copy2, null);
            eventSet = new HashSet<>();
            eventSpace.remove(0);
            eventSet.add(eventData);
            eventSpace.add(eventData);
            queue.add(copy2);
        }
        return true;
    }
    @Override
    public SSAStateMatrix getNextState(SSAStateMatrix startingState, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {

        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateStateSpace(endingState, new Random(), stateSpace,queue, stateSet)){
            return endingState;
        }

        return null;

    }

    @Override
    public SSAStateMatrix getNextEventState(SSAStateMatrix startingState, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {

        SSAStateMatrix endingState = new SSAStateMatrix(startingState);

        if(updateEventSpace(endingState, new Random(), eventSpace,event,queue,copy, eventSet)){
            return endingState;
        }

        return null;
    }


    public Node getNode() {
        return this.node;
    }

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
