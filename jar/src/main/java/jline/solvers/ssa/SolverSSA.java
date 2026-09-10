/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.GlobalConstants;
import jline.lang.Event;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SolverType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.State;
import jline.lang.state.EventCache;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import jline.streaming.Collector;
import jline.streaming.StreamingOptions;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.io.InputOutput.line_debug;
import static jline.solvers.ssa.analyzers.Solver_ssa_analyzer.solver_ssa_analyzer;

public class SolverSSA extends NetworkSolver {

    private final int DEFAULT_THREADS = (int) FastMath.ceil(Runtime.getRuntime().availableProcessors() / 2.0);
    public ExecutorService threadPool;
    public int numThreads = DEFAULT_THREADS;
    public EventCache eventCache;
    private Collector streamingCollector;

    public SolverSSA(Network model) {
        // If no options provided, use default options
        this(model, new SSAOptions());
        this.result = new SSAResult();
    }

    public SolverSSA(Network model, Object... args) {
        super(model, "SolverSSA");
        this.setOptions(Solver.parseOptions(new SolverOptions(SolverType.SSA), args));
        this.result = new SSAResult();
    }

    /**
     * Creates a new SolverSSA that warm-starts the simulated trajectory from
     * the steady-state solution of an auxiliary solver (see
     * {@link jline.solvers.NetworkSolver#initFromSolver}).
     *
     * @param model The network model to analyze
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @param args Variable arguments for solver options
     */
    public SolverSSA(Network model, NetworkSolver initSolver, Object... args) {
        this(model, args);
        this.initFromSolver(initSolver);
    }

    public SolverSSA(Network model, String method) {
        super(model, "SolverSSA", SolverSSA.defaultOptions().method(method));
        this.result = new SSAResult();
    }

    public SolverSSA(Network model, SolverOptions options) {
        super(model, "SolverSSA", options);
        this.result = new SSAResult();
    }

    /**
     * Returns the feature set supported by the SSA solver
     *
     * @return - the feature set supported by the SSA solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source", "Router",
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "Cache", "CacheClassSwitcher", "CacheRetrieval",
                // see _kb/06-solver-catalog.md for rationale
                "Place", "Transition", "Linkage", "Enabling", "Inhibiting", "Timing", "Firing", "Storage",
                "MAP", "MMPP2", "MMAP", "APH", "PH", "Replayer",
                "Coxian", "Erlang", "Exp", "HyperExp",
                "Det", "Gamma", "Lognormal", "Pareto", "Uniform", "Weibull",  // converted via snNonmarkovToPh
                "StatelessClassSwitcher", "InfiniteServer",
                "SharedServer", "Buffer", "Dispatcher",
                // Finite capacity regions: the NRM sample path honours them.
                "Region",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS",
                "SchedStrategy_GPS", "SchedStrategy_LPS", "SchedStrategy_SIRO",
                "SchedStrategy_HOL", "SchedStrategy_LCFS",
                "SchedStrategy_SEPT", "SchedStrategy_LEPT",
                "SchedStrategy_LCFSPR", "SchedStrategy_PSPRIO",
                "SchedStrategy_DPSPRIO", "SchedStrategy_GPSPRIO",
                "SchedStrategy_LCFSPRPRIO", "SchedStrategy_FCFSPRPRIO",
                "SchedStrategy_PAS", "SchedStrategy_OI", "SchedStrategy_POLLING",
                "RoutingStrategy_RROBIN",
                "RoutingStrategy_WRROBIN",
                "RoutingStrategy_JSQ",
                "RoutingStrategy_SQ",
                "RoutingStrategy_SDR",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_SFIFO", "ReplacementStrategy_LRU",
                "ReplacementStrategy_HLRU", "ReplacementStrategy_CLIMB", "ReplacementStrategy_QLRU",
                "SchedStrategy_EXT", "ClosedClass", "SelfLoopingClass", "OpenClass",
                "OpenSignal", "ClosedSignal",
                "SignalType_NEGATIVE", "SignalType_CATASTROPHE",
                "SignalBatchRemoval", "SignalRemovalPolicy",
                "Fork", "Join", "Forker", "Joiner",
                "Balking", "Reneging", "Retrial",
                "LoadDependence", "ClassDependence", "JointDependence", "GlobalDependence",
                // c-server stations and binding buffers: the serial engine walks
                // the same State arms SolverCTMC declares and the NRM honours both
                "MultiServer", "FiniteCapacity"
        });
        return featSupported;
    }

    public EventCache getEventCache() {
        return eventCache;
    }

    public void setEventCache(EventCache eventCache) {
        this.eventCache = eventCache;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public void setNumThreads(int numThreads) {
        this.numThreads = numThreads;
    }

    public NetworkStruct getStruct() {
        return this.model.getStruct(true);
    }

    public ExecutorService getThreadPool() {
        return threadPool;
    }

    public void setThreadPool(ExecutorService threadPool) {
        this.threadPool = threadPool;
    }

    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    public List<String> listValidMethods(Network model) {

        // SolverSSA.m:55 verbatim. 'ssa' is the reference's ALIAS for the serial
        // engine and 'para' the short spelling of 'parallel'; both used to be
        // wrong here -- 'ssa' and 'ssa.parallel' were advertised while
        // Solver_ssa_analyzer dispatched neither (they fell to its "Unknown
        // analysis method" arm), and 'para', which it does dispatch, was not
        // advertised at all.
        return Arrays.asList("default", "ssa", "serial", "para", "parallel", "nrm");
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        threadPool = Executors.newFixedThreadPool(this.numThreads);

        String origmethod = options.method;
        long T0 = java.lang.System.nanoTime();
        if (this.options == null) {
            this.options = new SolverOptions(SolverType.SSA);
        }
        // options.events (DES event budget) overrides options.samples when set;
        // samples remains accepted as a deprecated alias for the event budget.
        if (this.options.events > 0) {
            this.options.samples = this.options.events;
        }
        // Propagate solver verbose level to global
        GlobalConstants.Verbose = options.verbose;
        if (this.enableChecks && !supports(this.model)) {
            throw new RuntimeException("This model is not supported by the SSA solver.");
        }
        // see _kb/06-solver-catalog.md for rationale
        this.resetRandomGeneratorSeed(options.seed);

        // Native fork-join support: simulate the tag-augmented copy and fold
        // the auxiliary sibling classes back into the original classes
        boolean isFJ = false;
        for (jline.lang.constant.NodeType nt : this.sn.nodetype) {
            if (nt == jline.lang.constant.NodeType.Fork || nt == jline.lang.constant.NodeType.Join) {
                isFJ = true;
                break;
            }
        }
        jline.solvers.tr.FJTagTransform.Context fjctx = null;
        if (isFJ) {
            if ("parallel".equals(options.method) || "ssa.parallel".equals(options.method) || "nrm".equals(options.method)) {
                jline.io.InputOutput.line_warning(jline.io.InputOutput.mfilename(new Object(){}), "The " + options.method + " method does not support fork-join models, switching to the serial method.");
            }
            options.method = "serial";
            fjctx = jline.solvers.tr.FJTagTransform.expand(this.model, this.sn, options.verbose, "SSA");
            this.sn = fjctx.fjsn;
        }

        String method = options.method;
        line_debug(options.verbose, String.format("SSA solver starting: method=%s, samples=%d, seed=%d",
            method, options.samples, options.seed));
        line_debug(options.verbose, "Running SSA simulation, calling solver_ssa_analyzer");
        SSAResult result;
        try {
            result = solver_ssa_analyzer(this.sn, this.options, this);
        } catch (RuntimeException e) {
            // carry the reason in the message: the getAvg boundary reports only
            // getMessage(), so a bare wrapper hides why the model was refused
            throw new RuntimeException("SSA simulation failed: " + e.getMessage(), e);
        }
        
        // Validate that the analyzer returned valid results
        Matrix QN = result.QN;
        Matrix UN = result.UN;
        Matrix RN = result.RN;
        Matrix TN = result.TN;
        Matrix CN = result.CN;
        Matrix XN = result.XN;
        
        // Validate that essential result matrices are not null or empty
        Map<Integer, Matrix> tranSysState = result.tranSysState;
        Matrix tranSync = result.tranSync;
        NetworkStruct sn = result.sn;

        // see _kb/06-solver-catalog.md for rationale
        Matrix rtOrig = sn.rt.copy();
        if (!isFJ) {
            // (skipped on fork-join models: the analyzed struct is the
            // tag-augmented copy, whose stateful set does not fit the original model)
            for (int isf = 0; isf < sn.nstateful; isf++) {
                StatefulNode statefulNode = this.model.getStatefulNodes().get(isf);
                if (statefulNode instanceof Cache) {
                    Cache cache = (Cache) statefulNode;
                    CacheNodeParam _cnp = (CacheNodeParam) sn.nodeparam.get(statefulNode);
                    cache.setResultHitProb(_cnp.actualhitprob);
                    cache.setResultMissProb(_cnp.actualmissprob);
                    if (_cnp.actualdelayedhitprob != null) {
                        cache.setResultDelayedHitProb(_cnp.actualdelayedhitprob);
                    }
                    cache.setResultResidT(_cnp.actualresidt);
                    this.model.refreshChains(true);
                }
            }
        }
        double runtime = result.runtime;
        int M = sn.nstations;
        int R = sn.nclasses;
        AvgHandle T = getAvgTputHandles();
        Matrix AN;
        if (isFJ) {
            // fold the auxiliary sibling classes back into the original
            // classes and report the Join per-sibling waiting time
            jline.solvers.tr.FJTagTransform.Lifted lifted =
                    jline.solvers.tr.FJTagTransform.lift(fjctx, QN, UN, RN, TN, CN, XN, T);
            QN = lifted.QN;
            UN = lifted.UN;
            RN = lifted.RN;
            TN = lifted.TN;
            CN = lifted.CN;
            XN = lifted.XN;
            AN = lifted.AN;
            // restore the original struct: downstream consumers (tables,
            // handles) index the folded matrices by the original classes
            this.sn = fjctx.snOrig;
        } else {
            Matrix rtRefreshed = sn.rt;
            sn.rt = rtOrig;
            jline.api.sn.SnPnAvgRates.snPnAvgRates(sn, QN, TN, null, RN);
            AN = snGetArvRFromTput(sn, TN, T);
            if (rtRefreshed != null) {
                sn.rt = rtRefreshed;
            }
        }
        Matrix WN = new Matrix(0, 0);

        // The analyzer already handles the "default/" convention, so just use the result
        this.result.method = result.method;
        // see _kb/06-solver-catalog.md for rationale
        ((SSAResult) this.result).tranSysState = result.tranSysState;
        ((SSAResult) this.result).tranSync = result.tranSync;
        ((SSAResult) this.result).sn = result.sn;
        this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, this.result.method, options.samples);

        // Transfer confidence interval data if available
        if (result instanceof SSAResult) {
            SSAResult ssaResult = (SSAResult) result;
            ((SSAResult) this.result).QNCI = ssaResult.QNCI;
            ((SSAResult) this.result).UNCI = ssaResult.UNCI;
            ((SSAResult) this.result).RNCI = ssaResult.RNCI;
            ((SSAResult) this.result).TNCI = ssaResult.TNCI;
            ((SSAResult) this.result).ANCI = ssaResult.ANCI;
            ((SSAResult) this.result).WNCI = ssaResult.WNCI;
            // Derived START/PREEMPT rates and the per-step tags. Kept in their
            // own fields: they are annotations on transitions the engine already
            // fires, not metrics, so they add no getAvgTable column.
            ((SSAResult) this.result).startRate = ssaResult.startRate;
            ((SSAResult) this.result).preemptRate = ssaResult.preemptRate;
            ((SSAResult) this.result).tranTags = ssaResult.tranTags;
            // HOW LONG THE RUN SHOULD HAVE BEEN, when the caller asked. The
            // batch-means half-width at the run's confidence level pins the
            // ASYMPTOTIC variance, which is the quantity a run length is planned
            // from -- not the stationary variance, which on M/M/1 differs from
            // it by a factor blowing up like (1-rho)^-2.
            this.result.runLengthPlan = planRunLength(options, this.result.QN,
                    ssaResult.QNCI, options.samples);
        }
    }

    /**
     * runAnalyzer with its checked exceptions wrapped: the two accessors below
     * are plain getters and a caller of getStartRate has no separate recovery
     * for a failure of the analyzer itself.
     */
    private void runAnalyzerChecked() {
        try {
            this.runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * (stations x classes) rate at which a class-r job BEGINS or RESUMES
     * holding a server at station i, estimated over the simulated path exactly
     * as the throughput is.
     *
     * <p>At a lossless station with no in-service abandonment
     * getStartRate == getAvgTput + getPreemptRate up to simulation error; the
     * CTMC accessor of the same name reports the exact value.</p>
     */
    public Matrix getStartRate() {
        if (this.result == null || ((SSAResult) this.result).startRate == null) {
            runAnalyzerChecked();
        }
        Matrix startRate = ((SSAResult) this.result).startRate;
        if (startRate == null) {
            throw new RuntimeException("This solver run produced no START rates.");
        }
        return startRate;
    }

    /**
     * (stations x classes) rate at which a class-r job HOLDING A SERVER at
     * station i is pushed back into the buffer. Zero at a non-preemptive station.
     */
    public Matrix getPreemptRate() {
        if (this.result == null || ((SSAResult) this.result).preemptRate == null) {
            runAnalyzerChecked();
        }
        Matrix preemptRate = ((SSAResult) this.result).preemptRate;
        if (preemptRate == null) {
            throw new RuntimeException("This solver run produced no PREEMPT rates.");
        }
        return preemptRate;
    }

    // set the number of threads for para SolverSSA
    public void setParallelism(int numThreads) {
        this.numThreads = numThreads;
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverSSA.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * All SSA methods are stochastic simulation.
     *
     * @param method the method name to classify
     * @return true always
     */
    @Override
    public boolean isStochasticMethod(String method) {
        return true;
    }

    /**
     * The fork-join model class, which EVERY SSA method has to clear.
     *
     * <p>runAnalyzer tag-augments a fork-join model through
     * {@code ModelAdapter.fjtag}, whose first act is {@code fjValidate}, so a
     * model that validator refuses is refused whichever method was asked for.
     * The feature set cannot state it -- Fork and Join are declared, and the
     * rules are about how they are WIRED (the pairing, the join strategy, the
     * tasks per link, whether an open class is routed through the fork) -- so it
     * is structural, and it is the validator's own body of rules rather than a
     * copy of them.</p>
     *
     * <p>Without it the report offered every ssa.* row on a fork-join model
     * whose Join names no fork, and each one then threw; SolverCTMC gates on the
     * same predicate for the same reason.</p>
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        String reason = super.supportsModelMethod(method);
        if (!reason.isEmpty() || this.model == null) {
            return reason;
        }
        String fj = jline.lang.ModelAdapter.fjSupportsReason(this.model.getStruct(false));
        if (!fj.isEmpty()) {
            return fj;
        }
        // 'nrm' is the one SSA method with a model class of its own, and the test
        // for it already existed: the dispatch consulted it to PREFER the NRM
        // while nothing consulted it to decide whether the name could be OFFERED,
        // so the report listed ssa.nrm on every model and an explicit request then
        // raised from Solver_ssa_analyzer_nrm.
        if ("nrm".equalsIgnoreCase(method) || "ssa.nrm".equalsIgnoreCase(method)) {
            return jline.solvers.ssa.analyzers.Solver_ssa_analyzer.nrmMethodRefusal(
                    this.model.getStruct(false));
        }
        return "";
    }

    /**
     * Sample node state evolution using SSA simulation
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState sample(jline.lang.nodes.Node node, Integer numEvents, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numEvents != null) {
                this.options.samples = numEvents;
            }
            
            // Only support serial method for sampling (as in MATLAB)
            String originalMethod = this.options.method;
            if (!this.options.method.equals("serial")) {
                this.options.method = "serial";
            }
            
            // Force reanalysis 
            this.options.force = true;
            
            // Run the analyzer to get transient state information
            SSAResult result = solver_ssa_analyzer(this.sn, this.options, this);
            
            if (result.tranSysState == null || result.tranSync == null) {
                throw new RuntimeException("Transient state data not available from SSA analyzer");
            }
            
            NetworkStruct sn = this.getStruct();
            int nodeIndex = node.getNodeIndex();
            
            // Get the stateful index for this node (using nodeToStateful like MATLAB)
            int statefulIndex = (int) sn.nodeToStateful.get(nodeIndex);
            
            // Create the sample result
            SampleNodeState sampleResult = new SampleNodeState();
            sampleResult.handle = node;
            sampleResult.isaggregate = false;
            
            // Get time points (index 0 in tranSysState)
            sampleResult.t = result.tranSysState.get(0);
            
            // Get state data for this node (index 1 + statefulIndex)
            sampleResult.state = result.tranSysState.get(1 + statefulIndex);
            
            // Process events from sync data - following MATLAB implementation
            List<jline.lang.Event> events = new ArrayList<>();
            Matrix tranSync = result.tranSync;
            
            for (int e = 0; e < tranSync.length(); e++) {
                int syncIndex = (int) tranSync.get(e) - 1; // tranSync stores 1-based indices (MATLAB convention)
                if (syncIndex >= 0 && syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
                    // Add active events
                    if (sn.sync.get(syncIndex).active != null) {
                        for (jline.lang.Event activeEvent : sn.sync.get(syncIndex).active.values()) {
                            // Create a copy of the event with the timestamp
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                activeEvent.getEvent(), 
                                activeEvent.getNode(), 
                                activeEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(activeEvent.getProb());
                            eventCopy.setState(activeEvent.getState());
                            eventCopy.setJob(activeEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                    
                    // Add passive events  
                    if (sn.sync.get(syncIndex).passive != null) {
                        for (jline.lang.Event passiveEvent : sn.sync.get(syncIndex).passive.values()) {
                            // Create a copy of the event with the timestamp
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                passiveEvent.getEvent(), 
                                passiveEvent.getNode(), 
                                passiveEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(passiveEvent.getProb());
                            eventCopy.setState(passiveEvent.getState());
                            eventCopy.setJob(passiveEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                }
            }
            
            sampleResult.event = events;
            
            // Handle active/passive marking if requested
            if (markActivePassive) {
                // Create array for active/passive event categorization
                List<jline.lang.Event[]> categorizedEvents = new ArrayList<>();
                int numTimePoints = sampleResult.t.length() - 1;
                
                for (int ti = 0; ti < numTimePoints; ti++) {
                    categorizedEvents.add(new jline.lang.Event[]{null, null}); // [active, passive]
                }
                
                // Categorize events by time and type
                for (jline.lang.Event event : events) {
                    // Find time index for this event
                    for (int ti = 0; ti < numTimePoints; ti++) {
                        if (Math.abs(event.getT() - sampleResult.t.get(ti)) < 1e-10) {
                            jline.lang.Event[] timeEvents = categorizedEvents.get(ti);
                            if (event.getEvent() == jline.lang.constant.EventType.ARV) {
                                timeEvents[1] = event; // passive
                            } else {
                                timeEvents[0] = event; // active  
                            }
                            break;
                        }
                    }
                }
                
                // Convert to simple event list (flatten active/passive structure)
                List<jline.lang.Event> flatEvents = new ArrayList<>();
                for (jline.lang.Event[] timeEvents : categorizedEvents) {
                    if (timeEvents[0] != null) flatEvents.add(timeEvents[0]);
                    if (timeEvents[1] != null) flatEvents.add(timeEvents[1]);
                }
                sampleResult.event = flatEvents;
            }
            
            return sampleResult;
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Sample node state evolution using SSA simulation
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate
     * @return SampleNodeState containing the sampling results  
     * @throws Exception if sampling fails
     */
    public SampleNodeState sample(jline.lang.nodes.Node node, int numEvents) throws Exception {
        return sample(node, numEvents, false);
    }

    /**
     * Sample node state evolution using SSA simulation using default sample count
     *
     * @param node The node to sample from
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails  
     */
    public SampleNodeState sample(jline.lang.nodes.Node node) throws Exception {
        return sample(node, null, false);
    }

    /**
     * Sample aggregated node state evolution using SSA simulation
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState sampleAggr(jline.lang.nodes.Node node, Integer numEvents, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numEvents != null) {
                this.options.samples = numEvents;
            }
            
            // Only support serial method for sampling (as in MATLAB)
            if (!this.options.method.equals("serial")) {
                this.options.method = "serial";
            }
            
            // Force reanalysis 
            this.options.force = true;
            
            // Run the analyzer to get transient state information
            SSAResult result = solver_ssa_analyzer(this.sn, this.options, this);
            
            if (result.tranSysState == null || result.tranSync == null) {
                throw new RuntimeException("Transient state data not available from SSA analyzer");
            }
            
            NetworkStruct sn = this.getStruct();
            int nodeIndex = node.getNodeIndex();
            
            // Get the stateful index for this node (using nodeToStateful like MATLAB)
            int statefulIndex = (int) sn.nodeToStateful.get(nodeIndex);
            
            // Create the sample result
            SampleNodeState sampleResult = new SampleNodeState();
            sampleResult.handle = node;
            sampleResult.isaggregate = true;
            
            // Get time points (index 0 in tranSysState)
            sampleResult.t = result.tranSysState.get(0);
            
            // Get state data for this node (index 1 + statefulIndex) and apply marginal aggregation
            Matrix nodeState = result.tranSysState.get(1 + statefulIndex);
            
            // Apply marginal aggregation as in MATLAB sampleAggr
            jline.lang.state.State.StateMarginalStatistics marginal = 
                jline.lang.state.ToMarginal.toMarginal(sn, nodeIndex, nodeState, null, null, null, null, null);
            
            sampleResult.state = marginal.nir;
            
            // Process events from sync data - following MATLAB implementation
            List<jline.lang.Event> events = new ArrayList<>();
            Matrix tranSync = result.tranSync;
            
            for (int e = 0; e < tranSync.length(); e++) {
                int syncIndex = (int) tranSync.get(e) - 1; // tranSync stores 1-based indices (MATLAB convention)
                if (syncIndex >= 0 && syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
                    // Add active events
                    if (sn.sync.get(syncIndex).active != null) {
                        for (jline.lang.Event activeEvent : sn.sync.get(syncIndex).active.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                activeEvent.getEvent(), 
                                activeEvent.getNode(), 
                                activeEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(activeEvent.getProb());
                            eventCopy.setState(activeEvent.getState());
                            eventCopy.setJob(activeEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                    
                    // Add passive events  
                    if (sn.sync.get(syncIndex).passive != null) {
                        for (jline.lang.Event passiveEvent : sn.sync.get(syncIndex).passive.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                passiveEvent.getEvent(), 
                                passiveEvent.getNode(), 
                                passiveEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(passiveEvent.getProb());
                            eventCopy.setState(passiveEvent.getState());
                            eventCopy.setJob(passiveEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                }
            }
            
            sampleResult.event = events;
            
            // Handle active/passive marking if requested
            if (markActivePassive) {
                List<jline.lang.Event[]> categorizedEvents = new ArrayList<>();
                int numTimePoints = sampleResult.t.length() - 1;
                
                for (int ti = 0; ti < numTimePoints; ti++) {
                    categorizedEvents.add(new jline.lang.Event[]{null, null}); // [active, passive]
                }
                
                // Categorize events by time and type
                for (jline.lang.Event event : events) {
                    for (int ti = 0; ti < numTimePoints; ti++) {
                        if (Math.abs(event.getT() - sampleResult.t.get(ti)) < 1e-10) {
                            jline.lang.Event[] timeEvents = categorizedEvents.get(ti);
                            if (event.getEvent() == jline.lang.constant.EventType.ARV) {
                                timeEvents[1] = event; // passive
                            } else {
                                timeEvents[0] = event; // active  
                            }
                            break;
                        }
                    }
                }
                
                // Convert to simple event list (flatten active/passive structure)
                List<jline.lang.Event> flatEvents = new ArrayList<>();
                for (jline.lang.Event[] timeEvents : categorizedEvents) {
                    if (timeEvents[0] != null) flatEvents.add(timeEvents[0]);
                    if (timeEvents[1] != null) flatEvents.add(timeEvents[1]);
                }
                sampleResult.event = flatEvents;
            }
            
            return sampleResult;
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Sample aggregated node state evolution using SSA simulation
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate
     * @return SampleNodeState containing the aggregated sampling results  
     * @throws Exception if sampling fails
     */
    public SampleNodeState sampleAggr(jline.lang.nodes.Node node, int numEvents) throws Exception {
        return sampleAggr(node, numEvents, false);
    }

    /**
     * Sample aggregated node state evolution using SSA simulation using default sample count
     *
     * @param node The node to sample from
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails  
     */
    public SampleNodeState sampleAggr(jline.lang.nodes.Node node) throws Exception {
        return sampleAggr(node, null, false);
    }

    /**
     * Sample node state with real-time streaming to OTLP receiver.
     * Streams phase-detailed state information during simulation.
     *
     * @param node The node to sample from
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState stream(jline.lang.nodes.Node node, StreamingOptions streamingOptions) throws Exception {
        return stream(node, null, streamingOptions);
    }

    /**
     * Sample node state with real-time streaming to OTLP receiver.
     * Streams phase-detailed state information during simulation.
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState stream(jline.lang.nodes.Node node, Integer numEvents, StreamingOptions streamingOptions) throws Exception {
        // Try to create collector, but proceed without streaming if it fails (e.g., gRPC issues)
        try {
            this.streamingCollector = new Collector(streamingOptions, this.sn);
        } catch (Throwable e) {
            // Collector creation failed (likely gRPC/OTLP issue) - proceed without streaming
            java.util.logging.Logger.getLogger(SolverSSA.class.getName())
                .log(java.util.logging.Level.WARNING,
                    "Streaming disabled due to initialization error: " + e.getMessage());
            this.streamingCollector = null;
        }

        try {
            SampleNodeState result = sample(node, numEvents, false);
            if (this.streamingCollector != null && result.t != null && result.t.length() > 0) {
                double finalTime = result.t.get(result.t.length() - 1);
                this.streamingCollector.flush(finalTime);
            }
            return result;
        } finally {
            if (this.streamingCollector != null) {
                this.streamingCollector.shutdown();
                this.streamingCollector = null;
            }
        }
    }

    /**
     * Sample aggregated node state with real-time streaming to OTLP receiver.
     * Streams state information aggregated by job class during simulation.
     *
     * @param node The node to sample from
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState streamAggr(jline.lang.nodes.Node node, StreamingOptions streamingOptions) throws Exception {
        return streamAggr(node, null, streamingOptions);
    }

    /**
     * Sample aggregated node state with real-time streaming to OTLP receiver.
     * Streams state information aggregated by job class during simulation.
     *
     * @param node The node to sample from
     * @param numEvents Number of samples to generate
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState streamAggr(jline.lang.nodes.Node node, Integer numEvents, StreamingOptions streamingOptions) throws Exception {
        // Try to create collector, but proceed without streaming if it fails (e.g., gRPC issues)
        try {
            this.streamingCollector = new Collector(streamingOptions, this.sn);
        } catch (Throwable e) {
            // Collector creation failed (likely gRPC/OTLP issue) - proceed without streaming
            java.util.logging.Logger.getLogger(SolverSSA.class.getName())
                .log(java.util.logging.Level.WARNING,
                    "Streaming disabled due to initialization error: " + e.getMessage());
            this.streamingCollector = null;
        }

        try {
            SampleNodeState result = sampleAggr(node, numEvents, false);
            if (this.streamingCollector != null && result.t != null && result.t.length() > 0) {
                double finalTime = result.t.get(result.t.length() - 1);
                this.streamingCollector.flush(finalTime);
            }
            return result;
        } finally {
            if (this.streamingCollector != null) {
                this.streamingCollector.shutdown();
                this.streamingCollector = null;
            }
        }
    }

    /**
     * Get the active streaming collector (for use by Solver_ssa.kt).
     * Returns null if streaming is not active.
     *
     * @return The active Collector, or null
     */
    public Collector getStreamingCollector() {
        return this.streamingCollector;
    }

    /**
     * Sample system-wide state evolution using SSA simulation
     *
     * @param numEvents Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleSysState containing the system-wide sampling results
     * @throws Exception if sampling fails
     */
    private SampleSysState sampleSysInternal(Integer numEvents, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numEvents != null) {
                this.options.samples = numEvents;
            } else {
                numEvents = this.options.samples;
            }
            
            // Only support serial method for sampling (as in MATLAB)
            if (!this.options.method.equals("serial")) {
                this.options.method = "serial";
            }
            
            // Force reanalysis 
            this.options.force = true;
            
            // Run the analyzer to get transient state information
            SSAResult result = solver_ssa_analyzer(this.sn, this.options, this);
            
            if (result.tranSysState == null || result.tranSync == null) {
                throw new RuntimeException("Transient state data not available from SSA analyzer");
            }
            
            NetworkStruct sn = this.getStruct();
            
            // Create the sample result
            SampleSysState sampleResult = new SampleSysState();
            sampleResult.handle = new ArrayList<>(this.model.getStatefulNodes());
            sampleResult.isaggregate = false;
            
            // Get time points (index 0 in tranSysState)
            sampleResult.t = result.tranSysState.get(0);
            
            // Get state data for all stateful nodes (indices 1 through nstateful+1)
            List<Matrix> stateList = new ArrayList<>();
            for (int i = 1; i <= sn.nstateful; i++) {
                Matrix nodeState = result.tranSysState.get(i);
                
                // Truncate if needed to match numEvents
                if (nodeState.getNumRows() > numEvents) {
                    Matrix truncatedTime = new Matrix(numEvents, 1);
                    Matrix truncatedState = new Matrix(numEvents, nodeState.getNumCols());
                    
                    for (int j = 0; j < numEvents; j++) {
                        truncatedTime.set(j, 0, sampleResult.t.get(j, 0));
                        for (int k = 0; k < nodeState.getNumCols(); k++) {
                            truncatedState.set(j, k, nodeState.get(j, k));
                        }
                    }
                    sampleResult.t = truncatedTime;
                    nodeState = truncatedState;
                }
                
                stateList.add(nodeState);
            }
            sampleResult.state = stateList;
            
            // Process events from sync data - following MATLAB implementation
            List<jline.lang.Event> events = new ArrayList<>();
            Matrix tranSync = result.tranSync;
            
            for (int e = 0; e < tranSync.length(); e++) {
                int syncIndex = (int) tranSync.get(e) - 1; // tranSync stores 1-based indices (MATLAB convention)
                if (syncIndex >= 0 && syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
                    // Add active events
                    if (sn.sync.get(syncIndex).active != null) {
                        for (jline.lang.Event activeEvent : sn.sync.get(syncIndex).active.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                activeEvent.getEvent(), 
                                activeEvent.getNode(), 
                                activeEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(activeEvent.getProb());
                            eventCopy.setState(activeEvent.getState());
                            eventCopy.setJob(activeEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                    
                    // Add passive events  
                    if (sn.sync.get(syncIndex).passive != null) {
                        for (jline.lang.Event passiveEvent : sn.sync.get(syncIndex).passive.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                passiveEvent.getEvent(), 
                                passiveEvent.getNode(), 
                                passiveEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(passiveEvent.getProb());
                            eventCopy.setState(passiveEvent.getState());
                            eventCopy.setJob(passiveEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                }
            }
            
            sampleResult.event = events;
            
            // Handle active/passive marking if requested
            if (markActivePassive) {
                List<jline.lang.Event[]> categorizedEvents = new ArrayList<>();
                int numTimePoints = sampleResult.t.length() - 1;
                
                for (int ti = 0; ti < numTimePoints; ti++) {
                    categorizedEvents.add(new jline.lang.Event[]{null, null}); // [active, passive]
                }
                
                // Categorize events by time and type
                for (jline.lang.Event event : events) {
                    for (int ti = 0; ti < numTimePoints; ti++) {
                        if (Math.abs(event.getT() - sampleResult.t.get(ti)) < 1e-10) {
                            jline.lang.Event[] timeEvents = categorizedEvents.get(ti);
                            if (event.getEvent() == jline.lang.constant.EventType.ARV) {
                                timeEvents[1] = event; // passive
                            } else {
                                timeEvents[0] = event; // active  
                            }
                            break;
                        }
                    }
                }
                
                // Convert to simple event list (flatten active/passive structure)
                List<jline.lang.Event> flatEvents = new ArrayList<>();
                for (jline.lang.Event[] timeEvents : categorizedEvents) {
                    if (timeEvents[0] != null) flatEvents.add(timeEvents[0]);
                    if (timeEvents[1] != null) flatEvents.add(timeEvents[1]);
                }
                sampleResult.event = flatEvents;
            }
            
            return sampleResult;
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Sample system-wide state evolution using SSA simulation
     *
     * @param numEvents Number of samples to generate
     * @return SampleSysState containing the system-wide sampling results  
     * @throws Exception if sampling fails
     */
    public SampleResult sampleSys(int numEvents) {
        try {
            SampleSysState sysState = sampleSysInternal(numEvents, false);
            return convertToSampleResult(sysState, numEvents);
        } catch (Exception e) {
            return new SampleResult();
        }
    }
    
    private SampleSysState sampleSysInternal(int numEvents, boolean markActivePassive) throws Exception {
        return sampleSysInternal((Integer) numEvents, markActivePassive);
    }

    /**
     * Sample system-wide state evolution using SSA simulation using default sample count
     *
     * @return SampleSysState containing the system-wide sampling results
     * @throws Exception if sampling fails  
     */
    public SampleSysState sampleSys() throws Exception {
        return sampleSysInternal(null, false);
    }

    /**
     * Sample aggregated system-wide state evolution using SSA simulation
     *
     * @param numEvents Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleSysState containing the aggregated system-wide sampling results
     * @throws Exception if sampling fails
     */
    private SampleSysState sampleSysAggrInternal(Integer numEvents, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            if (numEvents == null) {
                numEvents = this.options.samples;
            } else {
                this.options.samples = numEvents;
            }
            
            // Only support serial method for sampling (as in MATLAB)
            if (!this.options.method.equals("serial")) {
                this.options.method = "serial";
            }
            
            // Force reanalysis 
            this.options.force = true;
            
            NetworkStruct sn = this.getStruct();
            
            // Run the analyzer to get transient state information
            SSAResult result = solver_ssa_analyzer(this.sn, this.options, this);
            
            if (result.tranSysState == null || result.tranSync == null) {
                throw new RuntimeException("Transient state data not available from SSA analyzer");
            }
            
            // Create the sample result
            SampleSysState sampleResult = new SampleSysState();
            sampleResult.handle = new ArrayList<>(this.model.getStatefulNodes());
            sampleResult.isaggregate = true;
            
            // Get time points (index 0 in tranSysState)
            sampleResult.t = result.tranSysState.get(0);
            
            // Get and aggregate state data for all stateful nodes
            List<Matrix> stateList = new ArrayList<>();
            for (int isf = 0; isf < sn.nstateful; isf++) {
                Matrix nodeState = result.tranSysState.get(1 + isf);
                
                // Truncate if needed to match numEvents
                if (nodeState.getNumRows() > numEvents) {
                    Matrix truncatedTime = new Matrix(numEvents, 1);
                    Matrix truncatedState = new Matrix(numEvents, nodeState.getNumCols());
                    
                    for (int j = 0; j < numEvents; j++) {
                        truncatedTime.set(j, 0, sampleResult.t.get(j, 0));
                        for (int k = 0; k < nodeState.getNumCols(); k++) {
                            truncatedState.set(j, k, nodeState.get(j, k));
                        }
                    }
                    sampleResult.t = truncatedTime;
                    nodeState = truncatedState;
                }
                
                // Apply marginal aggregation as in MATLAB
                int nodeIndex = ((Double) sn.statefulToNode.get(isf)).intValue();
                Matrix aggregatedState;
                try {
                    jline.lang.state.State.StateMarginalStatistics marginal =
                        jline.lang.state.ToMarginal.toMarginal(sn, nodeIndex, nodeState, null, null, null, null, null);
                    aggregatedState = (marginal != null && marginal.nir != null) ? marginal.nir : nodeState;
                } catch (Exception e) {
                    // If toMarginal fails, use the original nodeState
                    aggregatedState = nodeState;
                }
                stateList.add(aggregatedState);
            }
            sampleResult.state = stateList;
            
            // Process events from sync data - same as sampleSys
            List<jline.lang.Event> events = new ArrayList<>();
            Matrix tranSync = result.tranSync;
            
            for (int e = 0; e < tranSync.length(); e++) {
                int syncIndex = (int) tranSync.get(e) - 1; // tranSync stores 1-based indices (MATLAB convention)
                if (syncIndex >= 0 && syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
                    // Add active events
                    if (sn.sync.get(syncIndex).active != null) {
                        for (jline.lang.Event activeEvent : sn.sync.get(syncIndex).active.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                activeEvent.getEvent(), 
                                activeEvent.getNode(), 
                                activeEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(activeEvent.getProb());
                            eventCopy.setState(activeEvent.getState());
                            eventCopy.setJob(activeEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                    
                    // Add passive events  
                    if (sn.sync.get(syncIndex).passive != null) {
                        for (jline.lang.Event passiveEvent : sn.sync.get(syncIndex).passive.values()) {
                            jline.lang.Event eventCopy = new jline.lang.Event(
                                passiveEvent.getEvent(), 
                                passiveEvent.getNode(), 
                                passiveEvent.getJobClass()
                            );
                            eventCopy.setT(sampleResult.t.get(e));
                            eventCopy.setProb(passiveEvent.getProb());
                            eventCopy.setState(passiveEvent.getState());
                            eventCopy.setJob(passiveEvent.getJob());
                            events.add(eventCopy);
                        }
                    }
                }
            }
            
            sampleResult.event = events;
            
            // Handle active/passive marking if requested
            if (markActivePassive) {
                List<jline.lang.Event[]> categorizedEvents = new ArrayList<>();
                int numTimePoints = sampleResult.t.length() - 1;
                
                for (int ti = 0; ti < numTimePoints; ti++) {
                    categorizedEvents.add(new jline.lang.Event[]{null, null}); // [active, passive]
                }
                
                // Categorize events by time and type
                for (jline.lang.Event event : events) {
                    for (int ti = 0; ti < numTimePoints; ti++) {
                        if (Math.abs(event.getT() - sampleResult.t.get(ti)) < 1e-10) {
                            jline.lang.Event[] timeEvents = categorizedEvents.get(ti);
                            if (event.getEvent() == jline.lang.constant.EventType.ARV) {
                                timeEvents[1] = event; // passive
                            } else {
                                timeEvents[0] = event; // active  
                            }
                            break;
                        }
                    }
                }
                
                // Convert to simple event list (flatten active/passive structure)
                List<jline.lang.Event> flatEvents = new ArrayList<>();
                for (jline.lang.Event[] timeEvents : categorizedEvents) {
                    if (timeEvents[0] != null) flatEvents.add(timeEvents[0]);
                    if (timeEvents[1] != null) flatEvents.add(timeEvents[1]);
                }
                sampleResult.event = flatEvents;
            }
            
            return sampleResult;
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Sample aggregated system-wide state evolution using SSA simulation
     *
     * @param numEvents Number of samples to generate
     * @return SampleSysState containing the aggregated system-wide sampling results  
     * @throws Exception if sampling fails
     */
    public SampleResult sampleSysAggr(int numEvents) {
        try {
            SampleSysState sysState = sampleSysAggrInternal(numEvents, false);
            return convertToSampleResult(sysState, numEvents);
        } catch (Exception e) {
            return new SampleResult();
        }
    }

    /**
     * Sample aggregated system-wide state evolution using SSA simulation using default sample count
     *
     * @return SampleSysState containing the aggregated system-wide sampling results
     * @throws Exception if sampling fails  
     */
    public SampleSysState sampleSysAggr() throws Exception {
        return sampleSysAggrInternal(null, false);
    }

    /**
     * Get probability for a specific node state
     *
     * @param node The node to get probability for
     * @param state The state vector (optional - uses node's default state if null)
     * @return Probability of being in the specified state
     * @throws Exception if probability calculation fails
     */
    public double getProb(jline.lang.nodes.Node node, Matrix state) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Follow MATLAB behavior - switch to serial for default/nrm methods
            String originalMethod = this.options.method;
            if (this.options.method.equals("default") || this.options.method.equals("nrm")) {
                this.options.method = "serial";
            }
            
            // Force reanalysis
            this.options.force = true;
            
            // Run the analyzer to get transient state information
            SSAResult result = solver_ssa_analyzer(this.sn, this.options, this);
            
            if (result.tranSysState == null) {
                throw new RuntimeException("Transient state data not available from SSA analyzer");
            }
            
            NetworkStruct sn = this.getStruct();
            int nodeIndex = node.getNodeIndex();
            
            // Get the stateful index for this node (using nodeToStateful like MATLAB)
            int statefulIndex = (int) sn.nodeToStateful.get(nodeIndex);
            
            // Get time points and state data
            Matrix timePoints = result.tranSysState.get(0);
            Matrix nodeState = result.tranSysState.get(1 + statefulIndex);
            
            // Create TSS matrix: [time_diffs, state_data]
            Matrix TSS = new Matrix(timePoints.getNumRows(), 1 + nodeState.getNumCols());
            
            // Calculate time differences (MATLAB: diff operation)
            TSS.set(0, 0, timePoints.get(0));
            for (int i = 1; i < timePoints.getNumRows(); i++) {
                TSS.set(i, 0, timePoints.get(i) - timePoints.get(i - 1));
            }
            
            // Copy state data
            for (int i = 0; i < nodeState.getNumRows(); i++) {
                for (int j = 0; j < nodeState.getNumCols(); j++) {
                    TSS.set(i, j + 1, nodeState.get(i, j));
                }
            }
            
            // Determine the target state
            Matrix targetState;
            if (state == null) {
                // Use default state from sn.state
                StatefulNode statefulNode = this.model.getStatefulNodes().get(statefulIndex);
                Matrix nodeDefaultState = sn.state.get(statefulNode);
                if (nodeDefaultState.getNumRows() > 1) {
                    throw new RuntimeException("There are multiple station states, choose an initial state as a parameter to getProb.");
                }
                targetState = nodeDefaultState;
            } else {
                targetState = state;
            }
            
            // Add padding of zeros for FCFS stations (MATLAB behavior)
            int stateColDiff = TSS.getNumCols() - 1 - targetState.getNumCols();
            if (stateColDiff > 0) {
                Matrix paddedState = new Matrix(1, targetState.getNumCols() + stateColDiff);
                // Fill with zeros first (padding)
                for (int j = 0; j < stateColDiff; j++) {
                    paddedState.set(0, j, 0.0);
                }
                // Copy original state
                for (int j = 0; j < targetState.getNumCols(); j++) {
                    paddedState.set(0, j + stateColDiff, targetState.get(0, j));
                }
                targetState = paddedState;
            }
            
            // Extract state portion of TSS (columns 2 to end)
            Matrix stateData = new Matrix(TSS.getNumRows(), TSS.getNumCols() - 1);
            for (int i = 0; i < stateData.getNumRows(); i++) {
                for (int j = 0; j < stateData.getNumCols(); j++) {
                    stateData.set(i, j, TSS.get(i, j + 1));
                }
            }
            
            // Find matching rows
            List<Integer> matchingRows = Matrix.findRows(stateData, targetState);
            
            if (!matchingRows.isEmpty()) {
                // Calculate probability: sum(time_diffs_of_matching_rows) / sum(all_time_diffs)
                double numerator = 0.0;
                double denominator = 0.0;
                
                for (int i = 0; i < TSS.getNumRows(); i++) {
                    double timeDiff = TSS.get(i, 0);
                    denominator += timeDiff;
                    if (matchingRows.contains(i)) {
                        numerator += timeDiff;
                    }
                }
                
                return numerator / denominator;
            } else {
                // State was not seen during simulation
                line_debug(options.verbose, "The state was not seen during the simulation.");
                return 0.0;
            }
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Get probability for a specific node state using the node's default state
     *
     * @param node The node to get probability for
     * @return Probability of being in the node's default state
     * @throws Exception if probability calculation fails
     */
    public double getProb(jline.lang.nodes.Node node) throws Exception {
        return getProb(node, null);
    }

    /**
     * Get marginal probability for a specific node state (by node index).
     * This overrides the base NetworkSolver method to use SSA sampling.
     *
     * @param node The node index to get probability for
     * @param state The state vector (optional - uses node's default state if null)
     * @return Probability result for being in the specified state
     */
    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        try {
            // Convert node index to Node object
            if (node < 0 || node >= this.model.getNumberOfNodes()) {
                return new ProbabilityResult(Double.NaN);
            }
            jline.lang.nodes.Node nodeObj = this.model.getNodes().get(node);
            double prob = getProb(nodeObj, state);
            return new ProbabilityResult(prob);
        } catch (Exception e) {
            return new ProbabilityResult(Double.NaN);
        }
    }

    /**
     * Get aggregated probability for a specific node state
     *
     * @param node The node to get probability for
     * @param state The state vector (optional - uses node's default state if null)
     * @return Aggregated probability of being in the specified state
     * @throws Exception if probability calculation fails
     */
    public double getProbAggr(jline.lang.nodes.Node node, Matrix state) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();

        try {
            // Force reanalysis
            this.options.force = true;

            // Get aggregated system state sample
            SampleSysState tranSysStateAggr;
            try {
                tranSysStateAggr = this.sampleSysAggr();
            } catch (Exception e) {
                line_debug(options.verbose, "SSA getProbAggr: sampleSysAggr failed (" + e.getMessage() + "), returning 0.");
                return 0.0;
            }

            // Validate sample result
            if (tranSysStateAggr == null || tranSysStateAggr.t == null ||
                tranSysStateAggr.state == null || tranSysStateAggr.state.isEmpty()) {
                line_debug(options.verbose, "SSA getProbAggr: no sample data available, returning 0.");
                return 0.0;
            }

            NetworkStruct sn = this.getStruct();
            int nodeIndex = node.getNodeIndex();

            // Get the stateful index for this node (using nodeToStateful like MATLAB)
            int statefulIndex = (int) sn.nodeToStateful.get(nodeIndex);

            // Validate stateful index
            if (statefulIndex < 0 || statefulIndex >= tranSysStateAggr.state.size()) {
                line_debug(options.verbose, "SSA getProbAggr: invalid stateful index, returning 0.");
                return 0.0;
            }

            // Get time points and aggregated state data for this node
            Matrix timePoints = tranSysStateAggr.t;
            Matrix nodeState = tranSysStateAggr.state.get(statefulIndex);

            // Validate time points and node state
            if (timePoints == null || timePoints.getNumRows() == 0 ||
                nodeState == null || nodeState.getNumRows() == 0) {
                line_debug(options.verbose, "SSA getProbAggr: no state samples available for node, returning 0.");
                return 0.0;
            }

            // Create TSS matrix: [time_diffs, state_data]
            Matrix TSS = new Matrix(timePoints.getNumRows(), 1 + nodeState.getNumCols());
            
            // Calculate time differences (MATLAB: diff operation)
            TSS.set(0, 0, timePoints.get(0));
            for (int i = 1; i < timePoints.getNumRows(); i++) {
                TSS.set(i, 0, timePoints.get(i) - timePoints.get(i - 1));
            }
            
            // Copy state data
            for (int i = 0; i < nodeState.getNumRows(); i++) {
                for (int j = 0; j < nodeState.getNumCols(); j++) {
                    TSS.set(i, j + 1, nodeState.get(i, j));
                }
            }
            
            // Extract state portion of TSS (columns 2 to end)
            Matrix stateData = new Matrix(TSS.getNumRows(), TSS.getNumCols() - 1);
            for (int i = 0; i < stateData.getNumRows(); i++) {
                for (int j = 0; j < stateData.getNumCols(); j++) {
                    stateData.set(i, j, TSS.get(i, j + 1));
                }
            }

            // Determine the target state (matching MATLAB: state = sn.state{isf})
            Matrix targetState;
            if (state == null) {
                StatefulNode statefulNode = this.model.getStatefulNodes().get(statefulIndex);
                targetState = sn.state.get(statefulNode);
                if (targetState == null || targetState.isEmpty()) {
                    line_debug(options.verbose, "SSA getProbAggr: state not set for node, returning 0.");
                    return 0.0;
                }
            } else {
                targetState = state;
            }
            // Convert to marginal per-class job counts if dimensions don't match
            // This handles the case where state is in detailed format but stateData is aggregated
            if (targetState.getNumCols() != stateData.getNumCols()) {
                try {
                    State.StateMarginalStatistics margStats = jline.lang.state.ToMarginal.toMarginal(sn, nodeIndex, targetState, null, null, null, null, null);
                    if (margStats != null && margStats.nir != null && !margStats.nir.isEmpty()) {
                        targetState = margStats.nir;
                    }
                } catch (Exception e) {
                    line_debug(options.verbose, "SSA getProbAggr: toMarginal conversion failed, returning 0.");
                    return 0.0;
                }
            }

            // Verify dimensions match after conversion, and targetState is valid
            if (targetState.getNumRows() == 0 || targetState.getNumCols() != stateData.getNumCols()) {
                line_debug(options.verbose, "SSA getProbAggr: state dimensions don't match stateData after conversion, returning 0.");
                return 0.0;
            }

            // Find matching rows - ensure both matrices have valid dimensions
            if (stateData.getNumRows() == 0) {
                line_debug(options.verbose, "SSA getProbAggr: no sample state data available, returning 0.");
                return 0.0;
            }

            List<Integer> matchingRows = Matrix.findRows(stateData, targetState);
            
            if (!matchingRows.isEmpty()) {
                // Calculate probability: sum(time_diffs_of_matching_rows) / sum(all_time_diffs)
                double numerator = 0.0;
                double denominator = 0.0;
                
                for (int i = 0; i < TSS.getNumRows(); i++) {
                    double timeDiff = TSS.get(i, 0);
                    denominator += timeDiff;
                    if (matchingRows.contains(i)) {
                        numerator += timeDiff;
                    }
                }
                
                return numerator / denominator;
            } else {
                // State was not seen during simulation
                line_debug(options.verbose, "The state was not seen during the simulation.");
                return 0.0;
            }
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Get aggregated probability for a specific node state using the node's default state
     *
     * @param node The node to get probability for
     * @return Aggregated probability of being in the node's default state
     * @throws Exception if probability calculation fails
     */
    public double getProbAggr(jline.lang.nodes.Node node) throws Exception {
        return getProbAggr(node, null);
    }

    /**
     * Aggregated probability for a node state, addressed by NODE INDEX.
     *
     * <p>The Node-typed pair above is the implementation; this overload is what
     * {@link jline.solvers.NetworkSolver#getProbAggr(int, Matrix)} declares and
     * what every index-addressed caller reaches, LineCLI's {@code -a prob-aggr}
     * included. Without it the base class answered "getProbAggr is not supported
     * by SolverSSA" for a solver that implements it, and a delegating client read
     * that refusal as a missing feature.</p>
     *
     * @param node  the node index
     * @param state per-class job counts, or null for the node's own state
     * @return scalar probability in [0,1]
     */
    @Override
    public ProbabilityResult getProbAggr(int node, Matrix state) {
        try {
            if (node < 0 || node >= this.model.getNumberOfNodes()) {
                return new ProbabilityResult(Double.NaN);
            }
            return new ProbabilityResult(getProbAggr(this.model.getNodes().get(node), state));
        } catch (Exception e) {
            return new ProbabilityResult(Double.NaN);
        }
    }

    /**
     * Get system-wide probability for the current system state
     *
     * @return Probability of being in the system state
     * @throws Exception if probability calculation fails
     */
    public ProbabilityResult getProbSys() {
        try {
            double prob = getProbSysInternal();
            return new ProbabilityResult(prob);
        } catch (Exception e) {
            return new ProbabilityResult(Double.NaN);
        }
    }
    
    private double getProbSysInternal() throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Force reanalysis 
            this.options.force = true;
            
            // Get system state sample
            SampleSysState tranSysState = this.sampleSys();
            
            // Combine time and state data
            List<Matrix> allStateData = new ArrayList<>();
            allStateData.add(tranSysState.t);
            allStateData.addAll(tranSysState.state);
            
            // Create combined matrix with all state information
            int totalCols = 1; // Start with 1 for time column
            for (Matrix stateMatrix : tranSysState.state) {
                totalCols += stateMatrix.getNumCols();
            }
            
            Matrix TSS = new Matrix(tranSysState.t.getNumRows(), totalCols);
            
            // Calculate time differences and copy to first column
            TSS.set(0, 0, tranSysState.t.get(0));
            for (int i = 1; i < tranSysState.t.getNumRows(); i++) {
                TSS.set(i, 0, tranSysState.t.get(i) - tranSysState.t.get(i - 1));
            }
            
            // Copy all state data to remaining columns
            int colOffset = 1;
            for (Matrix stateMatrix : tranSysState.state) {
                for (int i = 0; i < stateMatrix.getNumRows(); i++) {
                    for (int j = 0; j < stateMatrix.getNumCols(); j++) {
                        TSS.set(i, colOffset + j, stateMatrix.get(i, j));
                    }
                }
                colOffset += stateMatrix.getNumCols();
            }
            
            NetworkStruct sn = this.getStruct();
            
            // Build target system state vector from sn.state - follow MATLAB logic
            List<Double> targetStateList = new ArrayList<>();
            
            for (int isf = 0; isf < sn.nstateful; isf++) {
                StatefulNode statefulNode = this.model.getStatefulNodes().get(isf);
                Matrix nodeState = sn.state.get(statefulNode);
                
                if (nodeState.getNumRows() > 1) {
                    throw new RuntimeException("There are multiple station states, choose an initial state as a parameter to getProb.");
                }
                
                Matrix stateMatrix = tranSysState.state.get(isf);
                
                // Add padding of zeros for FCFS stations (MATLAB behavior)
                int paddingSize = stateMatrix.getNumCols() - nodeState.getNumCols();
                for (int p = 0; p < paddingSize; p++) {
                    targetStateList.add(0.0);
                }
                
                // Add the actual state values
                for (int j = 0; j < nodeState.getNumCols(); j++) {
                    targetStateList.add(nodeState.get(0, j));
                }
            }
            
            // Convert target state list to matrix
            Matrix targetState = new Matrix(1, targetStateList.size());
            for (int i = 0; i < targetStateList.size(); i++) {
                targetState.set(0, i, targetStateList.get(i));
            }
            
            // Extract state portion of TSS (columns 2 to end)
            Matrix stateData = new Matrix(TSS.getNumRows(), TSS.getNumCols() - 1);
            for (int i = 0; i < stateData.getNumRows(); i++) {
                for (int j = 0; j < stateData.getNumCols(); j++) {
                    stateData.set(i, j, TSS.get(i, j + 1));
                }
            }
            
            // Find matching rows
            List<Integer> matchingRows = Matrix.findRows(stateData, targetState);
            
            if (!matchingRows.isEmpty()) {
                // Calculate probability: sum(time_diffs_of_matching_rows) / sum(all_time_diffs)
                double numerator = 0.0;
                double denominator = 0.0;
                
                for (int i = 0; i < TSS.getNumRows(); i++) {
                    double timeDiff = TSS.get(i, 0);
                    denominator += timeDiff;
                    if (matchingRows.contains(i)) {
                        numerator += timeDiff;
                    }
                }
                
                return numerator / denominator;
            } else {
                // State was not seen during simulation
                line_debug(options.verbose, "The state was not seen during the simulation.");
                return 0.0;
            }
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }

    /**
     * Get aggregated system-wide probability for the current system state
     *
     * @return Aggregated probability of being in the system state
     * @throws Exception if probability calculation fails
     */
    public ProbabilityResult getProbSysAggr() {
        try {
            double prob = getProbSysAggrInternal();
            return new ProbabilityResult(prob);
        } catch (Exception e) {
            return new ProbabilityResult(Double.NaN);
        }
    }
    
    private double getProbSysAggrInternal() throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Force reanalysis 
            this.options.force = true;
            
            // Get aggregated system state sample
            SampleSysState tranSysStateAggr = this.sampleSysAggr();
            
            // Combine time and state data
            List<Matrix> allStateData = new ArrayList<>();
            allStateData.add(tranSysStateAggr.t);
            allStateData.addAll(tranSysStateAggr.state);
            
            // Create combined matrix with all state information
            int totalCols = 1; // Start with 1 for time column
            for (Matrix stateMatrix : tranSysStateAggr.state) {
                totalCols += stateMatrix.getNumCols();
            }
            
            Matrix TSS = new Matrix(tranSysStateAggr.t.getNumRows(), totalCols);
            
            // Calculate time differences and copy to first column
            TSS.set(0, 0, tranSysStateAggr.t.get(0));
            for (int i = 1; i < tranSysStateAggr.t.getNumRows(); i++) {
                TSS.set(i, 0, tranSysStateAggr.t.get(i) - tranSysStateAggr.t.get(i - 1));
            }
            
            // Copy all state data to remaining columns
            int colOffset = 1;
            for (Matrix stateMatrix : tranSysStateAggr.state) {
                for (int i = 0; i < stateMatrix.getNumRows(); i++) {
                    for (int j = 0; j < stateMatrix.getNumCols(); j++) {
                        TSS.set(i, colOffset + j, stateMatrix.get(i, j));
                    }
                }
                colOffset += stateMatrix.getNumCols();
            }
            
            NetworkStruct sn = this.getStruct();
            
            // Build aggregated target state using State.toMarginal (following MATLAB logic)
            Matrix nir = new Matrix(sn.nstateful, sn.nclasses);
            
            for (int isf = 0; isf < sn.nstateful; isf++) {
                int ind = ((Double) sn.statefulToNode.get(isf)).intValue();
                StatefulNode statefulNode = this.model.getStatefulNodes().get(isf);
                Matrix nodeState = sn.state.get(statefulNode);
                
                if (nodeState.getNumRows() > 1) {
                    line_debug(options.verbose, String.format("Some states at node %d will be ignored. Please assign the node with a specific state.", ind));
                }
                
                // Use toMarginal to get marginal statistics - get first row if multiple rows
                Matrix stateRow = nodeState.getNumRows() > 1 ? Matrix.extractRows(nodeState, 0, 1) : nodeState;
                jline.lang.state.State.StateMarginalStatistics marginal = 
                    jline.lang.state.ToMarginal.toMarginal(sn, ind, stateRow, null, null, null, null, null);
                
                // Extract nir values (marginal state). toMarginal returns nir as a
                // 1xK ROW, so (r,0) reads past the end for every class after the
                // first and the exception surfaced as a NaN probability.
                for (int r = 0; r < sn.nclasses; r++) {
                    nir.set(isf, r, marginal.nir.get(0, r));
                }
            }
            
            // Transpose nir and flatten to row vector (following MATLAB: nir = nir'; nir(:)').
            // Column-major over the TRANSPOSE walks the classes of one stateful node
            // before moving to the next, which is also the order the sampled columns
            // are concatenated in above. Walking stateful-first instead compares the
            // target against a permutation of itself, so no row ever matched.
            Matrix targetState = new Matrix(1, sn.nstateful * sn.nclasses);
            int idx = 0;
            for (int isf = 0; isf < sn.nstateful; isf++) {
                for (int r = 0; r < sn.nclasses; r++) {
                    targetState.set(0, idx++, nir.get(isf, r));
                }
            }
            
            // Extract state portion of TSS (columns 2 to end)
            Matrix stateData = new Matrix(TSS.getNumRows(), TSS.getNumCols() - 1);
            for (int i = 0; i < stateData.getNumRows(); i++) {
                for (int j = 0; j < stateData.getNumCols(); j++) {
                    stateData.set(i, j, TSS.get(i, j + 1));
                }
            }
            
            // Find matching rows
            List<Integer> matchingRows = Matrix.findRows(stateData, targetState);
            
            if (!matchingRows.isEmpty()) {
                // Calculate probability: sum(time_diffs_of_matching_rows) / sum(all_time_diffs)
                double numerator = 0.0;
                double denominator = 0.0;
                
                for (int i = 0; i < TSS.getNumRows(); i++) {
                    double timeDiff = TSS.get(i, 0);
                    denominator += timeDiff;
                    if (matchingRows.contains(i)) {
                        numerator += timeDiff;
                    }
                }
                
                return numerator / denominator;
            } else {
                // State was not seen during simulation
                line_debug(options.verbose, "The state was not seen during the simulation.");
                return 0.0;
            }
            
        } finally {
            // Restore original options
            this.options = originalOptions;
        }
    }
    
    /**
     * Convert SampleSysState to SampleResult for API compatibility
     */
    private SampleResult convertToSampleResult(SampleSysState sysState, int numEvents) {
        if (sysState == null) {
            return new SampleResult();
        }
        
        // Convert event list to Matrix 
        Matrix eventMatrix = new Matrix(0, 0);
        if (sysState.event != null && !sysState.event.isEmpty()) {
            eventMatrix = new Matrix(sysState.event.size(), 3); // time, node, class
            for (int i = 0; i < sysState.event.size(); i++) {
                Event evt = sysState.event.get(i);
                eventMatrix.set(i, 0, evt.getT());
                eventMatrix.set(i, 1, evt.getNode());
                eventMatrix.set(i, 2, evt.getJobClass());
            }
        }
        
        return new SampleResult("ssa", sysState.t != null ? sysState.t : new Matrix(0, 0), 
                               sysState.state != null ? sysState.state : new ArrayList<>(), 
                               eventMatrix, sysState.isaggregate, numEvents);
    }
    
    /**
     * Rename existing methods to avoid conflicts
     */
    private SampleSysState sampleSysAggrInternal(int numEvents, boolean markActivePassive) throws Exception {
        return sampleSysAggrInternal((Integer) numEvents, markActivePassive);
    }

    /**
     * Returns the default solver options for the SSA solver.
     *
     * @return Default solver options with SolverType.SSA
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.SSA);
    }


    /**
     * The run length the caller would need for the precision they asked for.
     *
     * <p>Null unless {@code options.config.runLengthPlan} is set; it is either a
     * Double target RELATIVE precision or a map with keys {@code relprecision}
     * and {@code confidence}.
     *
     * @param options     the solver options
     * @param means       the measured means
     * @param ciHalfWidth their confidence-interval half-widths
     * @param samplesUsed the run length those half-widths came from
     * @return the plan, or null
     * @see jline.api.sim.SimRunlength#sim_runlength_plan
     */
    public static java.util.Map<String, Object> planRunLength(SolverOptions options, Matrix means,
                                                              Matrix ciHalfWidth,
                                                              double samplesUsed) {
        if (options == null || options.config == null || options.config.runLengthPlan == null
                || means == null || ciHalfWidth == null || ciHalfWidth.isEmpty()
                || samplesUsed <= 0) {
            return null;
        }
        double relPrecision = 0.05;
        double confidence = 0.95;
        Object spec = options.config.runLengthPlan;
        if (spec instanceof java.util.Map) {
            java.util.Map<?, ?> m = (java.util.Map<?, ?>) spec;
            if (m.get("relprecision") instanceof Number) {
                relPrecision = ((Number) m.get("relprecision")).doubleValue();
            }
            if (m.get("confidence") instanceof Number) {
                confidence = ((Number) m.get("confidence")).doubleValue();
            }
        } else if (spec instanceof Number && ((Number) spec).doubleValue() > 0) {
            relPrecision = ((Number) spec).doubleValue();
        }
        return jline.api.sim.SimRunlength.sim_runlength_plan(means, ciHalfWidth, samplesUsed,
                relPrecision, confidence);
    }

    /**
     * Not available: SolverSSA does not record per-job response times.
     *
     * <p>A simulator must report what it measured. The inherited NetworkSolver
     * implementation fabricates an exponential law with the right mean, which
     * carries no information about the tail and would be indistinguishable, to
     * the caller, from a measured distribution. SSA samples state trajectories,
     * not per-job sojourn times, so there is nothing to build an empirical CDF
     * from -- the reference {@code @SolverSSA/getCdfRespT.m} refuses by name
     * and so does this port.</p>
     */
    @Override
    public jline.io.Ret.DistributionResult getCdfRespT() {
        return getCdfRespT((AvgHandle) null);
    }

    @Override
    public jline.io.Ret.DistributionResult getCdfRespT(AvgHandle R) {
        throw new RuntimeException("SolverSSA does not record per-job response times, so it cannot "
                + "return an empirical response time CDF. Use SolverJMT for a measured CDF, or "
                + "getPerctRespT(...,'forktail') for the analytical fork-join tail.");
    }
}
