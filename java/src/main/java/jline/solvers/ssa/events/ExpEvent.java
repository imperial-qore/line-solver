package jline.solvers.ssa.events;

import jline.lang.HasSchedStrategy;
import jline.lang.JobClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.lang.nodes.StatefulNode;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;

import java.util.ArrayList;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

public class ExpEvent extends PhaseEvent implements NodeEvent {
    public final int statefulIndex;
    public final int classIndex;
    public final SchedStrategy schedStrategy;
    public final boolean isSource;
    public Node node;
    public final JobClass jobClass;
    public boolean isProcessorSharing;

    public final Distribution serviceProcess;

    public final DepartureEvent departureEvent;

    public ExpEvent(Node node, JobClass jobClass, DepartureEvent departureEvent) {
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
        if (!(distServiceProcess instanceof Exp)) {
            throw new RuntimeException("Exp distribution required");
        }

        this.departureEvent = departureEvent;
        this.serviceProcess = ((HasSchedStrategy) node).getServiceProcess(jobClass);

        this.isProcessorSharing = this.schedStrategy == SchedStrategy.PS;
    }

    @Override
    public long getNPhases() {
        return 1;
    }

    @Override
    public double getRate(SSAStateMatrix networkState) {
        int activeServers = 1;

        if (this.isProcessorSharing) {
            double serviceRatio = (double) networkState.getState(this.statefulIndex, this.classIndex)/(double) networkState.totalStateAtNode(this.statefulIndex);
            serviceRatio *= networkState.psTotalCapacity(this.statefulIndex);
            return this.serviceProcess.getRate()*serviceRatio;
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

        return this.serviceProcess.getRate()*activeServers;
    }

    @Override
    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        if (this.node instanceof StatefulNode) {
            if (this.node instanceof Source) {
                if (networkState.incrementPhase(this.statefulIndex, this.classIndex)) {
                    this.departureEvent.getNextState(networkState, stateSpace, queue, stateSet);
                }

                return true;
            } else if (networkState.getState(this.statefulIndex, this.classIndex) == 0) {
                return true;
            }
        }

        if (networkState.incrementPhase(this.statefulIndex, this.classIndex)) {
            this.departureEvent.getNextState(networkState, stateSpace, queue, stateSet);
            return true;
        }
        return true;
    }

    @Override
    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {
        if (this.node instanceof StatefulNode) {
            if (this.node instanceof Source) {
                if (networkState.incrementPhase(this.statefulIndex, this.classIndex)) {
                    this.departureEvent.getNextEventState(networkState, eventSpace,event, queue,copy, eventSet);
                }

                return true;
            } else if (networkState.getState(this.statefulIndex, this.classIndex) == 0) {
                return true;
            }
        }

        if (networkState.incrementPhase(this.statefulIndex, this.classIndex)) {
            this.departureEvent.getNextEventState(networkState, eventSpace,event, queue,copy, eventSet);
            return true;
        }
        return true;
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
}
