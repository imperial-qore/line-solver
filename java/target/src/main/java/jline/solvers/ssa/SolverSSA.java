package jline.solvers.ssa;

import jline.lang.*;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.ssa.events.DepartureEvent;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.metrics.Metrics;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.solvers.ssa.strategies.TauLeapingStateStrategy;
import jline.solvers.ssa.strategies.TauLeapingType;
import jline.util.Matrix;
//import jline.util.JLineAPI;

import java.util.*;

public class SolverSSA extends NetworkSolver {

    public double timeout;

    // transient filtering
    public boolean useMSER5;
    public boolean useR5;
    public int R5value;

    // metrics configurations
    public boolean recordMetricTimeline;
    public boolean disableResTime;
    public boolean disableTransientState;
    public EventStack eventStack;

    public SolverSSA(Network model) {
        this(model, new SolverOptions(SolverType.SSA));
    }

    public SolverSSA(Network model, SolverOptions options) {
        super(model, "SolverSSA", options);

        this.model = model;
        this.sn = model.getStruct(true);

        this.disableResTime = false;
        this.timeout = Double.POSITIVE_INFINITY;
        this.useMSER5 = false;
        this.useR5 = false;
        this.R5value = 19;
        this.recordMetricTimeline = true;
        this.disableTransientState = false;
    }

    private void initEventStack() {
        // loop through each node and add active events to the eventStack
        ListIterator<Node> nodeIter = this.model.getNodes().listIterator();
        int nodeIdx = -1;
        while (nodeIter.hasNext()) {
            Node node = nodeIter.next();
            if (!(node instanceof StatefulNode)) {
                continue;
            }

            nodeIdx++;
            Iterator<JobClass> jobClassIter = this.model.getClasses().listIterator();

            while (jobClassIter.hasNext()) {
                JobClass jobClass = jobClassIter.next();
                int jobClassIdx = jobClass.getJobClassIdx();
//                if (network.getClassLinks(node, jobClass) == 0) {
//                    this.simStruct.classcap[nodeIdx][jobClassIdx] = 0;
//                } else {
//                    double jobCap = jobClass.getNumberOfJobs();
//                    jobCap = Math.min(jobCap, node.getClassCap(jobClass));
//                    if ((jobCap == Double.POSITIVE_INFINITY) || (node.getDropStrategy() == DropStrategy.WaitingQueue)) {
//                        this.simStruct.classcap[nodeIdx][jobClassIdx] = Integer.MAX_VALUE;
//                    } else {
//                        this.simStruct.classcap[nodeIdx][jobClassIdx] = (int) jobCap;
//                    }
//                }
                Event dEvent = DepartureEvent.fromNodeAndClass(node, jobClass);
                this.eventStack.addEvent(dEvent);
                if (dEvent instanceof DepartureEvent) {
                    if (((DepartureEvent) dEvent).getPhaseEvent() != null) {
                        this.eventStack.addEvent(((DepartureEvent) dEvent).getPhaseEvent());
                    }
                }
            }

//            double nodeCap = node.getCap();
//            if (nodeCap == Double.POSITIVE_INFINITY) {
//                this.simStruct.cap[nodeIdx] = Integer.MAX_VALUE;
//            } else {
//                this.simStruct.cap[nodeIdx] = (int) nodeCap;
//            }
        }
    }

    public void runAnalyzer() throws IllegalAccessException {
        if (this.sn == null) {
            this.sn = this.model.getStruct(true);
        }

        this.random = new Random(this.options.seed);
        int samplesCollected = 1;
        int maxSamples = options.samples;
        double curTime = options.timespan[0];
        double maxTime = options.timespan[1];

        // Add ClosedClass instances to the reference station
        SSAStateMatrix networkState = new SSAStateMatrix(this.sn, this.random);
        for (JobClass jobClass : this.model.getClasses()) {
            if (jobClass instanceof ClosedClass) {
                int classIdx = this.model.getJobClassIndex(jobClass);
                ClosedClass cClass = (ClosedClass) jobClass;
                int stationIdx = this.model.getStatefulNodeIndex(cClass.getRefstat());
                networkState.setState(stationIdx, classIdx, (int) cClass.getPopulation());
                for (int i = 0; i < cClass.getPopulation(); i++) {
                    networkState.addToBuffer(stationIdx, classIdx);
                }
            }
        }

        if (this.options.config.tau_leaping == null) {
            this.eventStack = new EventStack();
        } else {
            this.eventStack = new TauLeapEventStack();

        }

        initEventStack();

        Timeline ssarunner = new Timeline(this.sn, this.options.seed);

        if (this.disableResTime) {
            ssarunner.disableResidenceTime();
        }

        if (this.disableTransientState) {
            ssarunner.disableTransientState();
        }

        if (this.useMSER5) {
            ssarunner.useMSER5();
        } else if (this.useR5) {
            ssarunner.useR5(this.R5value);
        }

        if (!this.recordMetricTimeline) {
            ssarunner.setMetricRecord(false);
        }

        if (this.options.config.tau_leaping != null) {
            TauLeapingType tauLeapingType = new TauLeapingType(this.options.config.tau_leaping.var_type,
                    this.options.config.tau_leaping.order_strategy,
                    this.options.config.tau_leaping.state_strategy,
                    this.options.config.tau_leaping.tau);

            ((TauLeapEventStack) this.eventStack).configureTauLeap(tauLeapingType);
            if ((tauLeapingType.getStateStrategy() == TauLeapingStateStrategy.TimeWarp) ||
                    (tauLeapingType.getStateStrategy() == TauLeapingStateStrategy.TauTimeWarp)) {
                ssarunner.cacheRecordings();
            }
        }

        double sysTime = 0;
        double startTime = System.currentTimeMillis();

        boolean isSteadyState;

        // collect samples and update states
        while ((samplesCollected < maxSamples) && (curTime < maxTime) && (sysTime < this.timeout)) {
            isSteadyState = curTime < this.options.timespan[1];

            if (this.options.config.tau_leaping != null) {
                curTime = ((TauLeapEventStack) this.eventStack).tauLeapUpdate(networkState, ssarunner, curTime, random);
            } else {
                curTime = this.eventStack.updateState(networkState, ssarunner, curTime, random);
            }

            if (isSteadyState && (curTime > this.options.timespan[1])) {
                ssarunner.resetHistory();
            }
            samplesCollected++;
            sysTime = (System.currentTimeMillis() - startTime) / 1000.0;

        }

        ssarunner.finalizeMetrics(curTime);

        this.result = new SolverResult();

        int M = this.sn.nstateful;
        int K = this.sn.nclasses;

        this.result.UN = new Matrix(M, K);
        this.result.QN = new Matrix(M, K);
        this.result.RN = new Matrix(M, K);
        this.result.TN = new Matrix(M, K);
        this.result.XN = new Matrix(1, K);
        this.result.CN = new Matrix(1, K);

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                Metrics metrics = ssarunner.getMetrics(i, r);
                this.result.QN.set(i, r, metrics.getMetricValueByName("Queue Length"));
                this.result.UN.set(i, r, metrics.getMetricValueByName("Utilization"));
                this.result.TN.set(i, r, metrics.getMetricValueByName("Throughput"));
                this.result.RN.set(i, r, metrics.getMetricValueByName("Response Time"));
            }
        }
    }

    public void enableMSER5() {
        this.useMSER5 = true;
        this.useR5 = false;
    }

    public void enableR5(int k) {
        this.useR5 = true;
        this.useMSER5 = false;
        this.R5value = k;
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.SSA);
    }
}
