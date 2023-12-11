package jline.solvers.ssa.events;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.lang.nodes.StatefulNode;
import jline.lang.processes.MAP;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

public class DepartureEvent extends Event implements NodeEvent {
    public int statefulIndex;
    public int classIndex;
    public boolean useBuffer;
    public SchedStrategy schedStrategy;
    public boolean isSource;
    public final Distribution serviceProcess;
    public Node node;
    public JobClass jobClass;
    public PhaseEvent phaseEvent;
    public boolean isMAP;
    public boolean isReference;
    public boolean isProcessorSharing;

    public static Event fromNodeAndClass(Node node, JobClass jobClass) {
        if (node instanceof HasSchedStrategy) {
            Distribution serviceDist = ((HasSchedStrategy)node).getServiceProcess(jobClass);
            if (serviceDist instanceof Erlang) {
                DepartureEvent depEvent = new DepartureEvent(node, jobClass);
                ErlangPhaseEvent ePhase = new ErlangPhaseEvent(node, jobClass, depEvent);
                depEvent.setPhaseEvent(ePhase);
                return ePhase;
            } else if (serviceDist instanceof MAP) {
                MAPPhaseEvent mapPhaseEvent = new MAPPhaseEvent(node, jobClass, (MAP) serviceDist);
                return new DepartureEvent(node, jobClass, mapPhaseEvent);
            } else if (serviceDist instanceof Coxian) {
                DepartureEvent depEvent = new DepartureEvent(node, jobClass);
                CoxianPhaseEvent cxPhaseEvent = new CoxianPhaseEvent(node, jobClass, depEvent);
                depEvent.setPhaseEvent(cxPhaseEvent);
                return cxPhaseEvent;
            } else if (serviceDist instanceof Exp) {
                DepartureEvent depEvent = new DepartureEvent(node, jobClass);
                ExpEvent activeEvent = new ExpEvent(node, jobClass, depEvent);
                depEvent.setPhaseEvent(activeEvent);
                return activeEvent;
            }

        }

        return new DepartureEvent(node, jobClass);
    }

    public DepartureEvent(Node node, JobClass jobClass) {
        super();
        this.node = node;
        this.jobClass = jobClass;


        if (node instanceof StatefulNode) {
            this.statefulIndex = ((StatefulNode)this.node).getStatefulIndex();
        } else {
            this.statefulIndex = -1;
        }

        if (node instanceof Source) {
            this.isReference = true;
        } else if (node instanceof ClassSwitch) {
            this.isReference = true;
        } else if (jobClass instanceof ClosedClass) {
            if (node == ((ClosedClass) jobClass).getRefstat()) {
                this.isReference = true;
            }
        } else {
            this.isReference = false;
        }

        this.classIndex = this.node.getModel().getJobClassIndex(this.jobClass);

        this.isSource = node instanceof Source;
        this.useBuffer = !this.isSource;
        this.schedStrategy = SchedStrategy.FCFS;
        if (node instanceof HasSchedStrategy) {
            this.schedStrategy = ((HasSchedStrategy)node).getSchedStrategy();
            this.serviceProcess = ((HasSchedStrategy)node).getServiceProcess(this.jobClass);
        } else {
            this.serviceProcess = new Immediate();
        }

        this.phaseEvent = null;
        this.isMAP = (this.serviceProcess) instanceof MAP;

        this.isProcessorSharing = this.schedStrategy == SchedStrategy.PS;
    }

    public DepartureEvent(Node node, JobClass jobClass, PhaseEvent phaseEvent) {
        this(node, jobClass);
        this.phaseEvent = phaseEvent;
    }

    public void setPhaseEvent(PhaseEvent phaseEvent) {
        this.phaseEvent = phaseEvent;
    }

    public PhaseEvent getPhaseEvent() {
        return this.phaseEvent;
    }

    @Override
    public double getRate(SSAStateMatrix networkState) {
        if (this.isProcessorSharing) {
            double serviceRatio = (double) networkState.getState(this.statefulIndex, this.classIndex)/(double) networkState.totalStateAtNode(this.statefulIndex);
            serviceRatio *= networkState.psTotalCapacity(this.statefulIndex);
            return this.serviceProcess.getRate()*serviceRatio;
        }


        int activeServers = 1;

        if (this.node instanceof StatefulNode) {
            activeServers = networkState.inProcess(this.statefulIndex, this.classIndex);
            if (this.node instanceof Source) {
                activeServers = 1;
            } else if (activeServers == 0) {
                return Double.NaN;
            }
        }

        if (this.serviceProcess instanceof Immediate) {
            return Double.POSITIVE_INFINITY;
        } else if (this.serviceProcess instanceof Exp) {
            return this.serviceProcess.getRate()*activeServers;
        } else if (this.serviceProcess instanceof Disabled) {
            return Double.NaN;
        } else if (this.serviceProcess instanceof Erlang) {
            // Rate logic should be handled by PhaseEvent
            return Double.NaN;
        } else if (this.serviceProcess instanceof MAP) {
            //System.out.format("Map phase: %f\n", ((MAP)this.serviceProcess).getDepartureRate(stateMatrix.getGlobalPhase(this.statefulIndex, this.classIndex))*activeServers);
            return ((MAP)this.serviceProcess).getDepartureRate(networkState.getGlobalPhase(this.statefulIndex, this.classIndex))*activeServers;
        } else if (this.serviceProcess instanceof PH) {
            return Double.NaN;
        } else if (this.serviceProcess instanceof APH) {
            // Rate logic should be handled by PhaseEvent
            return Double.NaN;
        }

        return Double.NaN;
    }



    @Override
    public boolean updateStateSpace(SSAStateMatrix networkState, Random random, ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        if (this.isMAP) {
            MAP MAP = (MAP)(this.serviceProcess);
            int nextPhase = MAP.getNextPhaseAfterDeparture(networkState.getGlobalPhase(this.statefulIndex, this.classIndex), random);
            networkState.updateGlobalPhase(this.statefulIndex, this.classIndex, nextPhase);
        }

        if (this.node instanceof Source) {
            ArrayList<Pair<SynchedEvent,Double>> eventArrayList = this.node.getOutputEvents(this.jobClass, random);

            for (Pair<SynchedEvent,Double> outputEventDoublePair : eventArrayList) {
                outputEventDoublePair.getLeft().getNextState(networkState, stateSpace,queue, stateSet);
            }
            return true;
        }
        SSAStateMatrix newNetworkState = new SSAStateMatrix(networkState);
        boolean res = newNetworkState.stateDeparture(this.statefulIndex, classIndex);
        if (!res) {
            return false;
        }

        ArrayList<Pair<SynchedEvent,Double>>  eventArrayList = this.node.getOutputEvents(this.jobClass, random);
        for (Pair<SynchedEvent,Double> outputEventDoublePair: eventArrayList) {
            if(outputEventDoublePair.getRight() == 0) {
                continue;
            }
            outputEventDoublePair.getLeft().getNextState(newNetworkState, stateSpace,queue, stateSet);
        }

        return true;
    }

    @Override
    public boolean updateEventSpace(SSAStateMatrix networkState, Random random, ArrayList<EventData> eventSpace, Event event, Queue<SSAStateMatrix> queue, SSAStateMatrix copy, Set<EventData> eventSet) {
        if (this.isMAP) {

            MAP MAP = (MAP)(this.serviceProcess);
            int nextPhase = MAP.getNextPhaseAfterDeparture(networkState.getGlobalPhase(this.statefulIndex, this.classIndex), random);
            networkState.updateGlobalPhase(this.statefulIndex, this.classIndex, nextPhase);
        }

        if (this.node instanceof Source) {
            ArrayList<Pair<SynchedEvent,Double>> eventArrayList = this.node.getOutputEvents(this.jobClass, random);

            for (Pair<SynchedEvent,Double> outputEventDoublePair : eventArrayList) {outputEventDoublePair.getLeft().getNextEventState(networkState, eventSpace,event,queue,copy, eventSet, outputEventDoublePair);
            }

            return true;
        }
        SSAStateMatrix newNetworkState = new SSAStateMatrix(networkState);
        boolean res = newNetworkState.stateDeparture(this.statefulIndex, classIndex);
        networkState = newNetworkState;
        if (!res) {
            return false;
        }
        ArrayList<Pair<SynchedEvent, Double>> eventArrayList = this.node.getOutputEvents(this.jobClass, random);
        for (Pair<SynchedEvent, Double> outputEventDoublePair : eventArrayList) {
            if(outputEventDoublePair.getRight() == 0) {
                continue;
            }
            outputEventDoublePair.getLeft().getNextEventState(networkState, eventSpace,event,queue,copy, eventSet, outputEventDoublePair);
        }

        return true;
    }


    @Override
    public void printSummary() {
        System.out.format("Departure event for %s at %s\n", this.jobClass.getName(), this.node.getName());
    }


    @Override
    public int getMaxRepetitions(SSAStateMatrix networkState) {
        if (this.node instanceof Source) {
            return Integer.MAX_VALUE;
        }

        return networkState.getState(this.statefulIndex, this.classIndex);
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

    public boolean isReference() {
        return this.isReference;
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
