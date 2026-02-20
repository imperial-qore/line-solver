/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

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

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.io.InputOutputKt.line_debug;
import static jline.solvers.ssa.analyzers.Solver_ssa_analyzerKt.solver_ssa_analyzer;

public class SolverSSA extends NetworkSolver {

    private final int DEFAULT_THREADS = (int) FastMath.ceil(Runtime.getRuntime().availableProcessors() / 2.0);
    private ExecutorService threadPool;
    private int numThreads = DEFAULT_THREADS;
    private EventCache eventCache;
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
                "Cache", "CacheClassSwitcher",
                "Place", "Transition",  // SPN support
                "MAP", "MMPP2", "APH", "PH", "BMAP",
                "Coxian", "Erlang", "Exp", "HyperExp", "Pareto",
                "StatelessClassSwitcher", "InfiniteServer",
                "SharedServer", "Buffer", "Dispatcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS",
                "SchedStrategy_GPS", "SchedStrategy_SIRO",
                "SchedStrategy_HOL", "SchedStrategy_LCFS",
                "SchedStrategy_SEPT", "SchedStrategy_LEPT",
                "SchedStrategy_LCFSPR", "SchedStrategy_PSPRIO",
                "SchedStrategy_DPSPRIO", "SchedStrategy_GPSPRIO",
                "RoutingStrategy_RROBIN",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_SFIFO", "ReplacementStrategy_LRU",
                "SchedStrategy_EXT", "ClosedClass", "SelfLoopingClass", "OpenClass"
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

        // Implementation of listValidMethods
        return Arrays.asList("default", "ssa", "serial", "ssa.parallel", "parallel", "nrm");
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        threadPool = Executors.newFixedThreadPool(this.numThreads);

        String origmethod = options.method;
        long T0 = java.lang.System.nanoTime();
        if (this.options == null) {
            this.options = new SolverOptions(SolverType.SSA);
        }
        if (this.enableChecks && !supports(this.model)) {
            throw new RuntimeException("This model is not supported by the SSA solver.");
        }
        this.resetRandomGeneratorSeed(options.seed);
        String method = options.method;
        line_debug(options.verbose, String.format("SSA solver starting: method=%s, samples=%d, seed=%d", 
            method, options.samples, options.seed));
        line_debug(options.verbose, "Running SSA simulation, calling solver_ssa_analyzer");
        SSAResult result;
        try {
            result = solver_ssa_analyzer(this.sn, this.options, this);
        } catch (RuntimeException e) {
            e.printStackTrace();
            throw new RuntimeException("SSA simulation failed.", e);
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

        // Save original rt before refreshChains modifies it in-place.
        // In MATLAB, sn is a value-copy struct unaffected by refreshChains;
        // in Java, sn is a reference, so we must save/restore manually.
        Matrix rtOrig = sn.rt.copy();
        for (int isf = 0; isf < sn.nstateful; isf++) {
            StatefulNode statefulNode = this.model.getStatefulNodes().get(isf);
            if (statefulNode instanceof Cache) {
                Cache cache = (Cache) statefulNode;
                cache.setResultHitProb(((CacheNodeParam) sn.nodeparam.get(statefulNode)).actualhitprob);
                cache.setResultMissProb(((CacheNodeParam) sn.nodeparam.get(statefulNode)).actualmissprob);
                this.model.refreshChains(true);
            }
        }
        double runtime = result.runtime;
        int M = sn.nstations;
        int R = sn.nclasses;
        AvgHandle T = getAvgTputHandles();
        sn.rt = rtOrig;
        Matrix AN = snGetArvRFromTput(sn, TN, T);
        Matrix WN = new Matrix(0, 0);

        // The analyzer already handles the "default/" convention, so just use the result
        this.result.method = result.method;
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
        }
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
     * Sample node state evolution using SSA simulation
     *
     * @param node The node to sample from
     * @param numSamples Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState sample(jline.lang.nodes.Node node, Integer numSamples, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numSamples != null) {
                this.options.samples = numSamples;
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
                int syncIndex = (int) tranSync.get(e);
                if (syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
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
     * @param numSamples Number of samples to generate
     * @return SampleNodeState containing the sampling results  
     * @throws Exception if sampling fails
     */
    public SampleNodeState sample(jline.lang.nodes.Node node, int numSamples) throws Exception {
        return sample(node, numSamples, false);
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
     * @param numSamples Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState sampleAggr(jline.lang.nodes.Node node, Integer numSamples, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numSamples != null) {
                this.options.samples = numSamples;
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
                int syncIndex = (int) tranSync.get(e);
                if (syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
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
     * @param numSamples Number of samples to generate
     * @return SampleNodeState containing the aggregated sampling results  
     * @throws Exception if sampling fails
     */
    public SampleNodeState sampleAggr(jline.lang.nodes.Node node, int numSamples) throws Exception {
        return sampleAggr(node, numSamples, false);
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
     * @param numSamples Number of samples to generate
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState stream(jline.lang.nodes.Node node, Integer numSamples, StreamingOptions streamingOptions) throws Exception {
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
            SampleNodeState result = sample(node, numSamples, false);
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
     * @param numSamples Number of samples to generate
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleNodeState containing the aggregated sampling results
     * @throws Exception if sampling fails
     */
    public SampleNodeState streamAggr(jline.lang.nodes.Node node, Integer numSamples, StreamingOptions streamingOptions) throws Exception {
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
            SampleNodeState result = sampleAggr(node, numSamples, false);
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
     * @param numSamples Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleSysState containing the system-wide sampling results
     * @throws Exception if sampling fails
     */
    private SampleSysState sampleSysInternal(Integer numSamples, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            // Set number of samples if provided
            if (numSamples != null) {
                this.options.samples = numSamples;
            } else {
                numSamples = this.options.samples;
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
                
                // Truncate if needed to match numSamples
                if (nodeState.getNumRows() > numSamples) {
                    Matrix truncatedTime = new Matrix(numSamples, 1);
                    Matrix truncatedState = new Matrix(numSamples, nodeState.getNumCols());
                    
                    for (int j = 0; j < numSamples; j++) {
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
                int syncIndex = (int) tranSync.get(e);
                if (syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
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
     * @param numSamples Number of samples to generate
     * @return SampleSysState containing the system-wide sampling results  
     * @throws Exception if sampling fails
     */
    public SampleResult sampleSys(int numSamples) {
        try {
            SampleSysState sysState = sampleSysInternal(numSamples, false);
            return convertToSampleResult(sysState, numSamples);
        } catch (Exception e) {
            return new SampleResult();
        }
    }
    
    private SampleSysState sampleSysInternal(int numSamples, boolean markActivePassive) throws Exception {
        return sampleSysInternal((Integer) numSamples, markActivePassive);
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
     * @param numSamples Number of samples to generate (overrides solver options if provided)
     * @param markActivePassive Whether to mark events as active/passive
     * @return SampleSysState containing the aggregated system-wide sampling results
     * @throws Exception if sampling fails
     */
    private SampleSysState sampleSysAggrInternal(Integer numSamples, boolean markActivePassive) throws Exception {
        SolverOptions originalOptions = (SolverOptions) this.options.copy();
        
        try {
            if (numSamples == null) {
                numSamples = this.options.samples;
            } else {
                this.options.samples = numSamples;
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
                
                // Truncate if needed to match numSamples
                if (nodeState.getNumRows() > numSamples) {
                    Matrix truncatedTime = new Matrix(numSamples, 1);
                    Matrix truncatedState = new Matrix(numSamples, nodeState.getNumCols());
                    
                    for (int j = 0; j < numSamples; j++) {
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
                int syncIndex = (int) tranSync.get(e);
                if (syncIndex < sn.sync.size() && sn.sync.get(syncIndex) != null) {
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
     * @param numSamples Number of samples to generate
     * @return SampleSysState containing the aggregated system-wide sampling results  
     * @throws Exception if sampling fails
     */
    public SampleResult sampleSysAggr(int numSamples) {
        try {
            SampleSysState sysState = sampleSysAggrInternal(numSamples, false);
            return convertToSampleResult(sysState, numSamples);
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
                
                // Extract nir values (marginal state)
                for (int r = 0; r < sn.nclasses; r++) {
                    nir.set(isf, r, marginal.nir.get(r, 0));
                }
            }
            
            // Transpose nir and flatten to row vector (following MATLAB: nir = nir'; nir(:)')
            Matrix targetState = new Matrix(1, sn.nstateful * sn.nclasses);
            int idx = 0;
            for (int r = 0; r < sn.nclasses; r++) {
                for (int isf = 0; isf < sn.nstateful; isf++) {
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
    private SampleResult convertToSampleResult(SampleSysState sysState, int numSamples) {
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
                               eventMatrix, sysState.isaggregate, numSamples);
    }
    
    /**
     * Rename existing methods to avoid conflicts
     */
    private SampleSysState sampleSysAggrInternal(int numSamples, boolean markActivePassive) throws Exception {
        return sampleSysAggrInternal((Integer) numSamples, markActivePassive);
    }

    /**
     * Returns the default solver options for the SSA solver.
     *
     * @return Default solver options with SolverType.SSA
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.SSA);
    }

}
