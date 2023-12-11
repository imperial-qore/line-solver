package jline.solvers.ssa;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.StatefulNode;
import jline.lang.processes.MAP;
import jline.lang.sections.OutputSection;
import jline.solvers.ssa.events.*;
import jline.solvers.ssa.metrics.Metric;
import jline.solvers.ssa.metrics.Metrics;
import jline.solvers.ssa.metrics.QueueLengthMetric;
import jline.solvers.ssa.metrics.ResidenceTimeMetric;
import jline.solvers.ssa.metrics.ResponseTimeMetric;
import jline.solvers.ssa.metrics.ThroughputMetric;
import jline.solvers.ssa.metrics.TotalClassMetric;
import jline.solvers.ssa.metrics.UtilizationMetric;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.lang.distributions.CumulativeDistribution;
import jline.util.Matrix;
import jline.util.Pair;

import java.util.*;

public class Timeline {
    /*
        Maintains a list of all events in the simulation, acts as an interface point for Metric objects,
            handles steady-state
     */
    protected List<List<Integer>[]> transientState;
    protected List<Event> eventTimeline;
    protected List<Double> timeList;
    protected int nstateful;
    protected int nstations;
    protected int nclasses;
    protected Matrix nservers;
    protected SchedStrategy[] schedStrategies;
    protected double maxTime;
    protected double timeCache;
    protected Metric[][][] metrics;
    protected TotalClassMetric[] totalClassMetrics; // total counts for each class
    protected SSAStateMatrix networkState;
    protected List<Pair<Event, Integer>> eventCache; // list of unapplied events, for tau leaping
    protected Map<Event, Integer> eventClassMap; // cache for the class index of each event
    protected Map<Event, Integer> eventNodeMap;  // cache for the node of each event
    protected double currentTime;
    protected double nextTime;
    protected boolean useMSER5;
    protected boolean useR5;
    protected int R5value;
    protected long seed;
    protected Random random;
    protected boolean metricRecord;
    protected boolean cacheRecordings;
    protected boolean recordTransientState;
    protected boolean inferTimes;

    public Timeline(NetworkStruct sn, long seed) {
        this.nstateful = sn.nstateful;
        this.nclasses = sn.nclasses;
        this.nservers = sn.nservers;
        this.seed = seed;
        this.random = new Random(seed);
        this.schedStrategies = new SchedStrategy[nstateful];
        for (int i = 0; i < nstateful; i++) {
            this.schedStrategies[i] = sn.sched.get(sn.stations.get(i));
        }
        this.eventTimeline = new ArrayList<Event>();
        this.transientState = new ArrayList<List<Integer>[]>();
        this.timeList = new ArrayList<Double>();

        this.useMSER5 = false;
        this.useR5 = false;
        this.R5value = 19;
        this.metricRecord = true;
        this.recordTransientState = true;

        this.metrics = new Metric[this.nstateful][this.nclasses][5];
        this.totalClassMetrics = new TotalClassMetric[this.nclasses];
        this.timeCache = 0;
        this.cacheRecordings = false;

        this.currentTime = 0;
        this.nextTime = 0;

        this.inferTimes = false;

        for (int i = 0; i < this.nclasses; i++) {
            this.totalClassMetrics[i] = new TotalClassMetric(i);
        }

        // Build all 5 metrics. In the future, it might be desirable to allow configuration of this
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                this.metrics[i][j][0] = new QueueLengthMetric(i, j, (int) this.nservers.get(i), this.metricRecord);
                this.metrics[i][j][1] = new UtilizationMetric(i, j, (int) this.nservers.get(i), this.metricRecord, sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
                if (!this.inferTimes) {
                    this.metrics[i][j][2] = new ResponseTimeMetric(i, j, (int) this.nservers.get(i), schedStrategies[i], this.metricRecord);
                    this.metrics[i][j][3] = new ResidenceTimeMetric(i, j, (int) this.nservers.get(i), schedStrategies[i], this.metricRecord, totalClassMetrics[j]);
                }
                this.metrics[i][j][4] = new ThroughputMetric(i, j, (int) this.nservers.get(i), this.metricRecord);

                for (int k = 0; k < 5; k++) {
                    if (this.metrics[i][j][k] == null) {
                        continue;
                    }
                    this.metrics[i][j][k].setRecord(this.metricRecord);
                    if (this.useMSER5) {
                        this.metrics[i][j][k].configureMSER5();
                    } else if (this.useR5) {
                        this.metrics[i][j][k].configureR5(this.R5value);
                    }
                }
            }
        }
        this.timeList.add(0.0);
        this.maxTime = 0;

        this.eventCache = new ArrayList<Pair<Event, Integer>>(this.nstateful * this.nclasses);
        this.eventClassMap = new HashMap<Event, Integer>();
        this.eventNodeMap = new HashMap<Event, Integer>();
    }

    public void disableResidenceTime() {
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                ((ResidenceTimeMetric) this.metrics[i][j][3]).disable();
            }
        }
    }

    public void cacheRecordings() {
        this.cacheRecordings = true;
    }

    public void useMSER5() {
        this.useMSER5 = true;
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int l = 0; l < 5; l++) {
                    this.metrics[i][j][l].configureMSER5();
                }
            }
        }
    }

    public void useR5(int k) {
        this.useR5 = true;
        this.R5value = k;
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int l = 0; l < 5; l++) {
                    this.metrics[i][j][l].configureR5(this.R5value);
                }
            }
        }
    }

    public void setMetricRecord(boolean record) {
        this.metricRecord = record;
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 5; k++) {
                    this.metrics[i][j][k].setRecord(this.metricRecord);
                }
            }
        }
    }

    public void disableTransientState() {
        this.recordTransientState = false;
    }

    public void setTime(double t) {
        this.currentTime = t;
    }

    public void setNextTime(double t) {
        this.nextTime = t;
    }

    public void record(double t, Event e, SSAStateMatrix networkState) {
        //this.eventTimeline.add(e);
        this.timeList.add(t);
        this.maxTime = t;

        if (e instanceof DepartureEvent) {
            if (((DepartureEvent) e).isReference()) {
                this.totalClassMetrics[((DepartureEvent) e).getClassIdx()].increment();
            }
            NodeEvent ne = (NodeEvent) e;

            if (!this.inferTimes) {
                for (int k = 2; k < 5; k++) {
                    this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][k].fromEvent(t, e);
                }
            }

            return;
        } else if (e instanceof SynchedEvent) {
            if (((SynchedEvent) e).isClassSwitched()) {
                this.totalClassMetrics[((SynchedEvent) e).getClassIdx()].increment();
            }
            return;
        } else if (e instanceof PhaseEvent) {
            if (this.recordTransientState) {
                this.transientState.add(networkState.getStateVectors());
            }
            return;
        }

        if (this.recordTransientState) {
            this.transientState.add(networkState.getStateVectors());
        }
        boolean foundNode = (e instanceof NodeEvent) && ((NodeEvent) e).isStateful();

        if (foundNode) {
            if (this.inferTimes) {
                this.metrics[((NodeEvent) e).getNodeStatefulIdx()][((NodeEvent) e).getClassIdx()][4].fromEvent(t, e);
            } else {
                for (int k = 2; k < 5; k++) {
                    this.metrics[((NodeEvent) e).getNodeStatefulIdx()][((NodeEvent) e).getClassIdx()][k].fromEvent(t, e);
                }
            }
        }

        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 2; k++) {
                    this.metrics[i][j][k].fromStateMatrix(t, networkState);
                }
                if (!foundNode) {
                    if (this.inferTimes) {
                        this.metrics[i][j][4].fromEvent(t, e);
                    } else {
                        for (int k = 2; k < 5; k++) {
                            this.metrics[i][j][k].fromEvent(t, e);
                        }
                    }
                }
            }
        }
    }

    public void record(double t, Event e, SSAStateMatrix networkState, int n) {
        //this.eventTimeline.add(e);
        if (this.recordTransientState) {
            this.transientState.add(networkState.getStateVectors());
        }
        this.timeList.add(t);

        if (e instanceof DepartureEvent) {
            if (((DepartureEvent) e).isReference()) {
                this.totalClassMetrics[((DepartureEvent) e).getClassIdx()].increment(n);
            }
            NodeEvent ne = (NodeEvent) e;

            if (!this.inferTimes) {
                for (int k = 2; k < 5; k++) {
                    this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][k].fromEvent(t, e, n);
                }
            }
            return;
        } else if (e instanceof SynchedEvent) {
            if (((SynchedEvent) e).isClassSwitched()) {
                this.totalClassMetrics[((SynchedEvent) e).getClassIdx()].increment(n);
            }
            return;
        } else if (!(e instanceof ArrivalEvent)) {
            return;
        } else if (e instanceof PhaseEvent) {
            return;
        }

        NodeEvent ne = (NodeEvent) e;

        if (!this.inferTimes) {
            for (int k = 2; k < 5; k++) {
                this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][k].fromEvent(t, e, n);
            }
        } else {
            this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][4].fromEvent(t, e);
        }
    }


    public void record(Event e, SSAStateMatrix networkState) {
        this.record(this.currentTime, e, networkState);
    }

    public void preRecord(double t, Event e, SSAStateMatrix networkState, int n) {
        if (!this.cacheRecordings) {
            this.timeCache = t;
            this.networkState = networkState;
            this.record(t, e, networkState, n);
            return;
        }
        if (n == 0) {
            this.timeCache = t;
            this.networkState = networkState;
            return;
        }
        if ((this.timeCache != t) && (this.networkState != null)) {
            this.recordCache();
        }
        this.networkState = networkState;
        this.timeCache = t;
        this.eventCache.add(new Pair<Event, Integer>(e, n));
    }


    public void preRecord(Event e, SSAStateMatrix networkState, int n) {
        this.preRecord(this.nextTime, e, networkState, n);
    }

    public void clearCache() {
        //this.eventCache = new ArrayList<Pair<Event,Integer>>();
        this.eventCache.clear();
    }

    public void recordCache() {
        this.currentTime = this.nextTime;
        this.maxTime = currentTime;

        if (!this.cacheRecordings) {
            if (this.networkState == null) {
                return;
            }
            for (int i = 0; i < this.nstateful; i++) {
                for (int j = 0; j < this.nclasses; j++) {
                    for (int k = 0; k < 2; k++) {
                        this.metrics[i][j][k].fromStateMatrix(this.currentTime, this.networkState);
                    }
                }
            }
            return;
        } else if (this.eventCache.isEmpty()) {
            return;
        }

        for (Pair<Event, Integer> ePair : this.eventCache) {
            Event e = ePair.getLeft();
            int n = ePair.getRight();
            double t = this.currentTime;
            if (this.recordTransientState) {
                this.transientState.add(networkState.getStateVectors());
            }

            //this.eventTimeline.add(e);
            this.timeList.add(t);

            if (e instanceof DepartureEvent) {
                if (((DepartureEvent) e).isReference()) {
                    this.totalClassMetrics[((DepartureEvent) e).getClassIdx()].increment(n);
                }
            } else if (e instanceof SynchedEvent) {
                if (((SynchedEvent) e).isClassSwitched()) {
                    this.totalClassMetrics[((SynchedEvent) e).getClassIdx()].increment(n);
                }
                continue;
            } else if (!(e instanceof ArrivalEvent)) {
                continue;
            }

            NodeEvent ne = (NodeEvent) e;

            if (!this.inferTimes) {
                for (int k = 2; k < 5; k++) {
                    this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][k].fromEvent(t, e, n);
                }
            } else {
                this.metrics[ne.getNodeStatefulIdx()][ne.getClassIdx()][4].fromEvent(t, e, n);
            }
        }

        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 2; k++) {
                    this.metrics[i][j][k].fromStateMatrix(this.currentTime, this.networkState);
                }
            }
        }

        this.clearCache();
    }


    public Metrics getMetrics(int nodeIdx, int classIdx) {
        Metrics mMetrics = new Metrics();
        for (int k = 0; k < 5; k++) {

            if (this.metrics[nodeIdx][classIdx] == null) {
                continue;
            }
            mMetrics.addMetric(this.metrics[nodeIdx][classIdx][k]);
        }
        return mMetrics;
    }

    public void finalizeMetrics(double t) {
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 5; k++) {
                    if (this.metrics[i][j][k] == null) {
                        continue;
                    }
                    this.metrics[i][j][k].taper(t);
                }
            }
        }
    }

    public List<Double> allQueueLengths() {
        List<Double> outList = new ArrayList<Double>(this.nstateful * this.nclasses);
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 5; k++) {
                    if (this.metrics[i][j][k] instanceof QueueLengthMetric) {
                        double mVal = ((QueueLengthMetric) this.metrics[i][j][k]).getMetric();
                        outList.add(mVal);
                    }
                }
            }
        }
        return outList;
    }

    public void resetHistory() {
        for (int i = 0; i < this.nstateful; i++) {
            for (int j = 0; j < this.nclasses; j++) {
                for (int k = 0; k < 5; k++) {
                    this.metrics[i][j][k].resetHistory();
                }
            }
        }
    }

    public SynchedEvent getLastOutputEvent() {
        for (int i = this.eventTimeline.size() - 1; i >= 0; i--) {
            if (this.eventTimeline.get(i) instanceof SynchedEvent) {
                return (SynchedEvent) this.eventTimeline.get(i);
            }
        }
        return null;
    }

    public boolean afterEvent(Event ev, SSAStateMatrix networkState) {

        if (ev instanceof ClassSwitchArrivalEvent) {
            ClassSwitchArrivalEvent this_ev = (ClassSwitchArrivalEvent) ev;
            CumulativeDistribution<JobClass> transitionCumulativeDistribution = new CumulativeDistribution<JobClass>(random);

            for (JobClass jobClassIter : this_ev.transitions.keySet()) {
                transitionCumulativeDistribution.addElement(jobClassIter, this_ev.transitions.get(jobClassIter));
            }

            this.record(this_ev, networkState);

            JobClass outClass = transitionCumulativeDistribution.sample(random);
            SynchedEvent synchedEvent = this_ev.node.getOutputEvent(outClass, random);
            this.record(synchedEvent, networkState);
            return this.afterEvent(synchedEvent, networkState);
        } else if (ev instanceof DepartureEvent) {
            DepartureEvent this_ev = (DepartureEvent) ev;
            if (this_ev.isMAP) {
                MAP MAP = (MAP) (this_ev.serviceProcess);
                int nextPhase = MAP.getNextPhaseAfterDeparture(networkState.getGlobalPhase(this_ev.statefulIndex, this_ev.classIndex), random);
                networkState.updateGlobalPhase(this_ev.statefulIndex, this_ev.classIndex, nextPhase);
            }

            if (this_ev.node instanceof Source) {
                if (this.afterEvent(this_ev.node.getOutputEvent(this_ev.jobClass, random), networkState)) {
                    this.record(this_ev, networkState);
                    return true;
                }
                return false;
            }

            boolean res = networkState.stateDeparture(this_ev.statefulIndex, this_ev.classIndex);
            if (!res) {
                return false;
            }

            this.afterEvent(this_ev.node.getOutputEvent(this_ev.jobClass, random), networkState);


            this.record(this_ev, networkState);

            return true;
        } else if (ev instanceof ErlangPhaseEvent) {
            ErlangPhaseEvent this_ev = (ErlangPhaseEvent) ev;
            if (this_ev.node instanceof StatefulNode) {
                if (this_ev.node instanceof Source) {
                    if (networkState.incrementPhase(this_ev.statefulIndex, this_ev.classIndex)) {
                        this.afterEvent(this_ev.departureEvent, networkState);
                        this.record(this_ev, networkState);
                    }

                    return true;
                } else if (networkState.getState(this_ev.statefulIndex, this_ev.classIndex) == 0) {
                    return true;
                }
            }

            if (networkState.incrementPhase(this_ev.statefulIndex, this_ev.classIndex)) {
                this.afterEvent(this_ev.departureEvent, networkState);
                this.record(this_ev, networkState);
                return true;
            }

            this.record(this_ev, networkState);
            return true;
        } else if (ev instanceof ExpEvent) {
            ExpEvent this_ev = (ExpEvent) ev;
            if (this_ev.node instanceof StatefulNode) {
                if (this_ev.node instanceof Source) {
                    if (networkState.incrementPhase(this_ev.statefulIndex, this_ev.classIndex)) {
                        this.afterEvent(this_ev.departureEvent, networkState);
                        this.record(this_ev, networkState);
                    }

                    return true;
                } else if (networkState.getState(this_ev.statefulIndex, this_ev.classIndex) == 0) {
                    return true;
                }
            }

            if (networkState.incrementPhase(this_ev.statefulIndex, this_ev.classIndex)) {
                this.afterEvent(this_ev.departureEvent, networkState);
                this.record(this_ev, networkState);
                return true;
            }

            this.record(this_ev, networkState);
//        System.out.println("hello");
            return true;
        } else if (ev instanceof ForkEvent) {
            ForkEvent this_ev = (ForkEvent) ev;
            for (OutputStrategy outputStrategy : this_ev.outputStrategies) {
                afterEvent(outputStrategy.getDestination().getArrivalEvent(this_ev.jobClass), networkState);
            }
            this.record(this_ev, networkState);
            return true;
        } else if (ev instanceof JoinEvent) {
            JoinEvent this_ev = (JoinEvent) ev;
            OutputSection os = this.getLastOutputEvent().getOutputSection();
            this_ev.seenSet.add(os);
            this_ev.waitingJobs.put(os, this_ev.waitingJobs.get(os) + 1);

            if (this_ev.seenSet.size() == this_ev.nSources) {
                for (OutputSection sOs : this_ev.waitingJobs.keySet()) {
                    this_ev.waitingJobs.put(sOs, this_ev.waitingJobs.get(sOs) - 1);
                }

                for (OutputStrategy outputStrategy : this_ev.outputSection.getOutputStrategies()) {
                    afterEvent(outputStrategy.getDestination().getArrivalEvent(this_ev.jobClass), networkState);
                }
            }

            return true;
        } else if (ev instanceof MAPPhaseEvent) {
            MAPPhaseEvent this_ev = (MAPPhaseEvent) ev;
            this.record(this_ev, networkState);

            int nextPhase = this_ev.MAP.getNextPhase(networkState.getGlobalPhase(this_ev.statefulIndex, this_ev.classIndex), random);
            return networkState.updateGlobalPhase(this_ev.statefulIndex, this_ev.classIndex, nextPhase);

        } else if (ev instanceof NodeArrivalEvent) {
            NodeArrivalEvent this_ev = (NodeArrivalEvent) ev;
            SynchedEvent nodeSynchedEvent = this_ev.node.getOutputEvent(this_ev.jobClass, random);
            return this.afterEvent(nodeSynchedEvent, networkState);
        } else if (ev instanceof CoxianPhaseEvent) {
            CoxianPhaseEvent this_ev = (CoxianPhaseEvent) ev;
            int nInPhase = 1;

            if (this_ev.node instanceof StatefulNode) {
                nInPhase = networkState.inProcess(this_ev.statefulIndex, this_ev.classIndex);
                if (this_ev.node instanceof Source) {
                    if (networkState.incrementPhase(this_ev.statefulIndex, this_ev.classIndex)) {
                        this.afterEvent(this_ev.departureEvent, networkState);
                        this.record(this_ev, networkState);
                    }

                    return true;
                } else if (nInPhase == 0) {
                    return true;
                }
            }

            CumulativeDistribution<Integer> startingPhaseCumulativeDistribution = new CumulativeDistribution<Integer>(random);
            int totalInPhase = 0;

            for (int i = 0; i < this_ev.subGenerator.length(); i++) {
                int inPhase = networkState.getInPhase(this_ev.statefulIndex, this_ev.classIndex, i);
                startingPhaseCumulativeDistribution.addElement(i, inPhase);
                totalInPhase += inPhase;
            }
            startingPhaseCumulativeDistribution.normalize(totalInPhase);

            int startingPhase = startingPhaseCumulativeDistribution.sample(random);

            CumulativeDistribution<Integer> endingPhaseCumulativeDistribution = new CumulativeDistribution<Integer>(random);
            double departureRate = -this_ev.subGenerator.get(startingPhase,startingPhase);
            double totalRate = -this_ev.subGenerator.get(startingPhase,startingPhase);

            for (int i = 0; i < this_ev.subGenerator.length(); i++) {
                if (i == startingPhase) {
                    continue;
                }

                departureRate -= this_ev.subGenerator.get(startingPhase,i);
                endingPhaseCumulativeDistribution.addElement(i, this_ev.subGenerator.get(startingPhase,i));
            }
            endingPhaseCumulativeDistribution.addElement(-1, departureRate);
            endingPhaseCumulativeDistribution.normalize(totalRate);

            int endingPhase = endingPhaseCumulativeDistribution.sample(random);

            if (endingPhase == -1) {
                networkState.updatePhase(this_ev.statefulIndex, this_ev.classIndex, startingPhase, -1);
                this.afterEvent(this_ev.departureEvent, networkState);
                this.record(this_ev, networkState);
                return true;
            }

            networkState.updatePhase(this_ev.statefulIndex, this_ev.classIndex, startingPhase, endingPhase);

            this.record(this_ev, networkState);
            return true;
        } else if (ev instanceof SinkArrivalEvent) {
            SinkArrivalEvent this_ev = (SinkArrivalEvent) ev;
            this.record(this_ev, networkState);
            return true;
        } else if (ev instanceof SynchedEvent) {
            SynchedEvent this_ev = (SynchedEvent) ev;
            this.record(this_ev, networkState);
            return this.afterEvent(this_ev.node.getArrivalEvent(this_ev.jobClass), networkState);
        } else if (ev instanceof ArrivalEvent) {
            ArrivalEvent this_ev = (ArrivalEvent) ev;
            if (this_ev.isStateful) {
                if (!networkState.stateArrival(this_ev.statefulIndex, this_ev.classIndex)) {
                    return false;
                }
            } else if (!(this_ev.node instanceof Sink)) {
                throw new RuntimeException(String.format("ArrivalEvent at %s not supported!", this_ev.node.getName()));
            }

            this.record(this_ev, networkState);
            return true;
        }
        return true;
    }

    public int afterEventN(int n, Event ev, SSAStateMatrix networkState) {
        int res = n;

        if (ev instanceof ClassSwitchArrivalEvent) {
            ClassSwitchArrivalEvent this_ev = (ClassSwitchArrivalEvent) ev;
            CumulativeDistribution<JobClass> transitionCumulativeDistribution = new CumulativeDistribution<JobClass>(random);
            Map<JobClass, Integer> transitionCount = new HashMap<JobClass, Integer>();

            for (JobClass jobClassIter : this_ev.transitions.keySet()) {
                transitionCumulativeDistribution.addElement(jobClassIter, this_ev.transitions.get(jobClassIter));
            }

            List<JobClass> jobClasses = this_ev.node.getModel().getClasses();

            for (JobClass jobClass : jobClasses) {
                transitionCount.put(jobClass, 0);
            }

            res = 0;
            for (int i = 0; i < n; i++) {
                JobClass selectedClass = transitionCumulativeDistribution.sample(random);
                transitionCount.put(selectedClass, transitionCount.get(selectedClass) + 1);
            }

            for (JobClass jobClass : jobClasses) {
                SynchedEvent synchedEvent = this_ev.node.getOutputEvent(jobClass, random);
                int nSwitched = this.afterEventN(transitionCount.get(jobClass), synchedEvent, networkState);
                this.preRecord(synchedEvent, networkState, nSwitched);
                res += nSwitched;
            }

            this.preRecord(this_ev, networkState, n - res);

            return res;
        } else if (ev instanceof DepartureEvent) {
            DepartureEvent this_ev = (DepartureEvent) ev;
            res = 0;

            if (this_ev.isMAP) {
                MAP MAP = (MAP) (this_ev.serviceProcess);
                int nextPhase = MAP.getNextPhaseAfterDeparture(networkState.getGlobalPhase(this_ev.statefulIndex, this_ev.classIndex), random);
                networkState.updateGlobalPhase(this_ev.statefulIndex, this_ev.classIndex, nextPhase);
            }

            if (this_ev.node instanceof Source) {
                res = this.afterEventN(n, this_ev.node.getOutputEvent(this_ev.jobClass, random), networkState);
            } else {
                res = networkState.stateDepartureN(n, this_ev.statefulIndex, this_ev.classIndex);
                this.afterEventN(n - res, this_ev.node.getOutputEvent(this_ev.jobClass, random), networkState);
            }

            this.preRecord(this_ev, networkState, n - res);

            return res;
        } else if (ev instanceof ErlangPhaseEvent) {
            ErlangPhaseEvent this_ev = (ErlangPhaseEvent) ev;
            if (this_ev.node instanceof StatefulNode) {
                if (this_ev.node instanceof Source) {
                    int nDepartures = networkState.incrementPhaseN(n, this_ev.statefulIndex, this_ev.classIndex);
                    int nRemDepartures = afterEventN(nDepartures, this_ev.departureEvent, networkState);
                    this.preRecord(this_ev, networkState, nDepartures - nRemDepartures);
                    return 0;
                } else if (networkState.getState(this_ev.statefulIndex, this_ev.classIndex) == 0) {
                    return 0;
                }
            }

            int nDepartures = networkState.incrementPhaseN(n, this_ev.statefulIndex, this_ev.classIndex);
            int nRemDepartures = this.afterEventN(nDepartures, this_ev.departureEvent, networkState);
            this.preRecord(this_ev, networkState, nDepartures - nRemDepartures);
            return nRemDepartures;
        } else if (ev instanceof ExpEvent) {
            ExpEvent this_ev = (ExpEvent) ev;
            if (this_ev.node instanceof StatefulNode) {
                if (this_ev.node instanceof Source) {
                    int nDepartures = networkState.incrementPhaseN(n, this_ev.statefulIndex, this_ev.classIndex);
                    int nRemDepartures = this.afterEventN(nDepartures, this_ev.departureEvent,networkState);
                    this.preRecord(this_ev, networkState, nDepartures - nRemDepartures);
                    return 0;
                } else if (networkState.getState(this_ev.statefulIndex, this_ev.classIndex) == 0) {
                    return 0;
                }
            }

            int nDepartures = networkState.incrementPhaseN(n, this_ev.statefulIndex, this_ev.classIndex);
            int nRemDepartures = afterEventN(nDepartures, this_ev.departureEvent, networkState);
            this.preRecord(this_ev, networkState, nDepartures - nRemDepartures);
            return nRemDepartures;
        } else if (ev instanceof ForkEvent) {
            ForkEvent this_ev = (ForkEvent) ev;
            res = 0;
            for (OutputStrategy outputStrategy : this_ev.outputStrategies) {
                res += this.afterEventN(n, outputStrategy.getDestination().getArrivalEvent(this_ev.jobClass), networkState);
            }
            this.record(this_ev, networkState);
            return res;
        } else if (ev instanceof JoinEvent) {
            JoinEvent this_ev = (JoinEvent) ev;
            int nUnapplied = 0;
            for (int i = 0; i < n; i++) {
                if (this.afterEvent(this_ev, networkState)) {
                    nUnapplied += (n - i);
                }
            }
            return nUnapplied;
        } else if (ev instanceof MAPPhaseEvent) {
            MAPPhaseEvent this_ev = (MAPPhaseEvent) ev;
            res = n;
            for (int i = 0; i < n; i++) {

                int nextPhase = this_ev.MAP.getNextPhase(networkState.getGlobalPhase(this_ev.statefulIndex, this_ev.classIndex), random);
                if (networkState.updateGlobalPhase(this_ev.statefulIndex, this_ev.classIndex, nextPhase)) {
                    res--;
                }
            }

            this.preRecord(this_ev, networkState, n - res);

            return res;
        } else if (ev instanceof NodeArrivalEvent) {
            NodeArrivalEvent this_ev = (NodeArrivalEvent) ev;
            SynchedEvent nodeSynchedEvent = this_ev.node.getOutputEvent(this_ev.jobClass, random);
            return afterEventN(n, nodeSynchedEvent, networkState);
        } else if (ev instanceof CoxianPhaseEvent) {
            CoxianPhaseEvent this_ev = (CoxianPhaseEvent) ev;
            res = n;
            for (int i = 0; i < n; i++) {
                if (this.afterEvent(this_ev, networkState)) {
                    res -= 1;
                } else {
                    return res;
                }
            }
            return res;
        } else if (ev instanceof SinkArrivalEvent) {
            SinkArrivalEvent this_ev = (SinkArrivalEvent) ev;
            this.preRecord(this_ev, networkState, n);
            return 0;
        } else if (ev instanceof SynchedEvent) {
            SynchedEvent this_ev = (SynchedEvent) ev;
            this.record(n, this_ev, networkState);
            return this.afterEventN(n, this_ev.node.getArrivalEvent(this_ev.jobClass), networkState);
        } else if (ev instanceof ArrivalEvent) {
            ArrivalEvent this_ev = (ArrivalEvent) ev;
            res = 0;

            if (this_ev.isStateful) {
                res = networkState.stateArrivalN(n, this_ev.statefulIndex, this_ev.classIndex);
            } else if (!(this_ev.node instanceof Sink)) {
                throw new RuntimeException(String.format("ArrivalEvent at %s not supported!", this_ev.node.getName()));
            }

            this.preRecord(this_ev, networkState, n); // NOT n-res to control the buffer.

            return res;
        }
        return res;
    }
}

