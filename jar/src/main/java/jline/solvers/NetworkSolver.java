/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.JobClass;
import jline.lang.processes.Distribution;
import jline.api.fj.FJ_tail_forktail;
import jline.io.LineCitations;
import jline.lang.processes.DistributionScaling;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.MatrixCell;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import static jline.api.pfqn.sens.Pfqn_sens.pfqn_sens;
import static jline.api.pfqn.sens.Pfqn_sens_linearizer.pfqn_sens_linearizer;
import jline.api.pfqn.sens.Pfqn_sens_mom;
import static jline.api.pfqn.sens.Pfqn_sens_mom.pfqn_sens_mom;
import static jline.api.pfqn.sens.Pfqn_sens_mva.pfqn_sens_mva;
import static jline.api.pfqn.sens.Pfqn_sens_mvaldmx.pfqn_sens_mvaldmx;
import static jline.api.pfqn.sens.Pfqn_sens_respt.pfqn_sens_respt;
import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.api.sn.SnGetProductFormParams.snGetProductFormParams;
import static jline.api.sn.SnHasProductForm.snHasProductForm;
import java.lang.reflect.Constructor;
import static jline.api.sn.SnGetDemandsChain.snGetDemandsChain;
import static jline.api.sn.SnGetNodeArvRFromTput.snGetNodeArvRFromTput;
import static jline.api.sn.SnGetNodeTputFromTput.snGetNodeTputFromTput;
import static jline.api.sn.SnGetResidTFromRespT.snGetResidTFromRespT;
import static jline.io.InputOutput.*;
import static jline.util.Utils.isInf;

/**
 * Abstract base class for solvers applicable to queueing network models.
 * <p>
 * This class provides the core functionality for analyzing queueing networks using various
 * solution algorithms. It manages performance metrics computation including queue lengths,
 * utilizations, response times, throughputs, and arrival rates at both steady-state and
 * transient conditions.
 * <p>
 * The solver operates on {@link Network} models and produces results through various
 * table formats for different levels of aggregation (station-level, node-level, chain-level).
 *
 * @see Network
 * @see SolverOptions
 * @see SolverResult
 */
public abstract class NetworkSolver extends Solver {

    // ========== Class Fields ==========
    
    /**
     * The queueing network model to be solved
     */
    public Network model;

    /**
     * Internal data structure describing the network model
     */
    public NetworkStruct sn;

    /**
     * Handles for steady-state average performance metrics
     */
    public SolverAvgHandles avgHandles;

    /**
     * Handles for transient performance metrics
     */
    public SolverTranHandles tranHandles;

    /**
     * Constructs a NetworkSolver with the specified model, name, and options.
     *
     * @param model   the queueing network model to solve
     * @param name    the name identifier for this solver instance
     * @param options configuration options for the solver
     * @throws RuntimeException if the model is empty (has no nodes)
     */
    protected NetworkSolver(Network model, String name, SolverOptions options) {
        super(name, options);
        this.model = model;
        // Allow null model for LayeredNetwork-based solvers (e.g., SolverLDES with LQN)
        if (model != null) {
            if (model.getNumberOfNodes() == 0) {
                throw new RuntimeException("The model supplied in input is empty.");
            }
            this.avgHandles = model.getAvgHandles();
            this.tranHandles = model.getTranHandles();
            this.sn = model.getStruct(true); // Force model to refresh
        }
    }

    /**
     * Constructs a NetworkSolver with the specified model and name using default options.
     *
     * @param model the queueing network model to solve
     * @param name  the name identifier for this solver instance
     */
    protected NetworkSolver(Network model, String name) {
        this(model, name, defaultOptions());
    }

    /**
     * Returns a list containing instances of all available network solvers for the given model.
     *
     * @param model the queueing network model
     * @return list of all available solver implementations
     */
    public static List<NetworkSolver> getAllSolvers(Network model) {
        SolverOptions options = new SolverOptions();
        List<NetworkSolver> solvers = new ArrayList<>();
        solvers.add(new SolverCTMC(model, options));
        solvers.add(new SolverLDES(model, options));
        solvers.add(new SolverFluid(model, options));
        solvers.add(new SolverJMT(model, options));
        solvers.add(new SolverMAM(model, options));
        solvers.add(new SolverMVA(model, options));
        solvers.add(new SolverNC(model, options));
        solvers.add(new SolverSSA(model, options));
        return solvers;
    }

    // Basic solver results
    public SolverResult avg() {
        return getAvg();
    }

    public SolverResult avg(SolverAvgHandles avgHandles) {
        return getAvg(avgHandles);
    }

    public SolverResult avg(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvg(Q, U, R, W, T, A);
    }

    // Arrival rates
    public Matrix avgArvR() {
        return getAvgArvR();
    }

    public Matrix avgArvRChain() {
        return getAvgArvRChain();
    }

    public AvgHandle avgArvRHandles() {
        return getAvgArvRHandles();
    }

    // Chain-level results
    public SolverResult avgChain() {
        return getAvgChain();
    }

    public NetworkAvgChainTable avgChainTable() {
        return getAvgChainTable();
    }

    public NetworkAvgChainTable avgChainTable(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    public NetworkAvgChainTable avgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgChainTable avgChainTable(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    public NetworkAvgChainTable avgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgChainTable avgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgChainTable -> avgChainT aliases
    public NetworkAvgChainTable avgChainT() {
        return avgChainTable();
    }

    public NetworkAvgChainTable avgChainT(SolverAvgHandles avgHandles) {
        return avgChainTable(avgHandles);
    }

    public NetworkAvgChainTable avgChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgChainTable avgChainT(boolean keepDisabled) {
        return avgChainTable(keepDisabled);
    }

    public NetworkAvgChainTable avgChainT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgChainTable avgChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Handles
    public SolverAvgHandles avgHandles() {
        return getAvgHandles();
    }

    // Node-level results
    public SolverResult avgNode() {
        return getAvgNode();
    }

    public Matrix avgNodeArvRChain() {
        return getAvgNodeArvRChain();
    }

    public SolverResult avgNodeChain() {
        return getAvgNodeChain();
    }

    public NetworkAvgNodeChainTable avgNodeChainTable() {
        return getAvgNodeChainTable();
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgNodeChainTable -> avgNodeChainT aliases
    public NetworkAvgNodeChainTable avgNodeChainT() {
        return avgNodeChainTable();
    }

    public NetworkAvgNodeChainTable avgNodeChainT(SolverAvgHandles avgHandles) {
        return avgNodeChainTable(avgHandles);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgNodeChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(boolean keepDisabled) {
        return avgNodeChainTable(keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgNodeChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    public Matrix avgNodeQLenChain() {
        return getAvgNodeQLenChain();
    }

    public Matrix avgNodeResidTChain() {
        return getAvgNodeResidTChain();
    }

    public Matrix avgNodeRespTChain() {
        return getAvgNodeRespTChain();
    }

    public NetworkAvgNodeTable avgNodeTable() {
        return getAvgNodeTable();
    }

    public NetworkAvgNodeTable avgNodeTable(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    public NetworkAvgNodeTable avgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeTable avgNodeTable(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgNodeTable -> avgNodeT aliases
    public NetworkAvgNodeTable avgNodeT() {
        return avgNodeTable();
    }

    public NetworkAvgNodeTable avgNodeT(SolverAvgHandles avgHandles) {
        return avgNodeTable(avgHandles);
    }

    public NetworkAvgNodeTable avgNodeT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgNodeTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeTable avgNodeT(boolean keepDisabled) {
        return avgNodeTable(keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgNodeTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    public Matrix avgNodeTputChain() {
        return getAvgNodeTputChain();
    }

    public Matrix avgNodeUtilChain() {
        return getAvgNodeUtilChain();
    }

    // Queue lengths
    public Matrix avgQLen() {
        return getAvgQLen();
    }

    public Matrix avgQLenChain() {
        return getAvgQLenChain();
    }

    public AvgHandle avgQLenHandles() {
        return getAvgQLenHandles();
    }

    // Residence times
    public Matrix avgResidT() {
        return getAvgResidT();
    }

    public Matrix avgResidTChain() {
        return getAvgResidTChain();
    }

    public AvgHandle avgResidTHandles() {
        return getAvgResidTHandles();
    }

    // Response times
    public Matrix avgRespT() {
        return getAvgRespT();
    }

    public Matrix avgRespTChain() {
        return getAvgRespTChain();
    }

    public AvgHandle avgRespTHandles() {
        return getAvgRespTHandles();
    }

    // System-level results
    public void avgSys() {
        getAvgSys();
    }

    public Matrix avgSysRespT() {
        return getAvgSysRespT();
    }

    public NetworkAvgSysTable avgSysTable() {
        return getAvgSysTable();
    }

    public NetworkAvgSysTable avgSysTable(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    public NetworkAvgSysTable avgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    // avgSysTable(boolean keepDisabled) - method signature not available
    // avgSysTable with handles and boolean - method signature not available
    // avgSysTable with 6 handles and boolean - method signature not available

    // avgSysTable -> avgSysT aliases
    public NetworkAvgSysTable avgSysT() {
        return avgSysTable();
    }

    public NetworkAvgSysTable avgSysT(SolverAvgHandles avgHandles) {
        return avgSysTable(avgHandles);
    }

    public NetworkAvgSysTable avgSysT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgSysTable(Q, U, R, W, T, A);
    }

    public Matrix avgSysTput() {
        return getAvgSysTput();
    }

    /**
     * Returns a fluent configurator for this solver's options, allowing chained
     * calls such as {@code solver.options().method("parallel").samples(20000)
     * .seed(1).build()}. The terminal {@code build()} returns this solver.
     *
     * @return a {@link SolverConfigurator} wrapping this solver
     */
    public SolverConfigurator options() {
        return new SolverConfigurator(this);
    }

    /**
     * Fluent configurator returned by {@link NetworkSolver#options()}.
     */
    public static class SolverConfigurator {
        private final NetworkSolver solver;

        public SolverConfigurator(NetworkSolver solver) {
            this.solver = solver;
        }

        public SolverConfigurator method(String method) {
            this.solver.options.method(method);
            return this;
        }

        public SolverConfigurator samples(int samples) {
            this.solver.options.samples(samples);
            return this;
        }

        public SolverConfigurator seed(int seed) {
            this.solver.options.seed(seed);
            return this;
        }

        public NetworkSolver build() {
            return this.solver;
        }
    }

    // Tables
    public NetworkAvgTable avgTable() {
        return getAvgTable();
    }

    public NetworkAvgTable avgTable(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    public NetworkAvgTable avgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable avgTable(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    public NetworkAvgTable avgTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable avgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Table -> T aliases
    public NetworkAvgTable avgT() {
        return avgTable();
    }

    public NetworkAvgTable avgT(SolverAvgHandles avgHandles) {
        return avgTable(avgHandles);
    }

    public NetworkAvgTable avgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable avgT(boolean keepDisabled) {
        return avgTable(keepDisabled);
    }

    public NetworkAvgTable avgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable avgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aT aliases for getAvgTable
    public NetworkAvgTable aT() {
        return getAvgTable();
    }

    public NetworkAvgTable aT(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    public NetworkAvgTable aT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable aT(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    public NetworkAvgTable aT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable aT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Throughputs
    public Matrix avgTput() {
        return getAvgTput();
    }

    public Matrix avgTputChain() {
        return getAvgTputChain();
    }

    public AvgHandle avgTputHandles() {
        return getAvgTputHandles();
    }

    // Utilizations
    public Matrix avgUtil() {
        return getAvgUtil();
    }

    public Matrix avgUtilChain() {
        return getAvgUtilChain();
    }

    public AvgHandle avgUtilHandles() {
        return getAvgUtilHandles();
    }

    // Wait times
    public Matrix avgWaitT() {
        return getAvgWaitT();
    }

    // Distribution functions
    public DistributionResult cdfPassT() {
        return getCdfPassT();
    }

    public DistributionResult cdfPassT(AvgHandle R) {
        return getCdfPassT(R);
    }

    public DistributionResult cdfRespT() {
        return getCdfRespT();
    }

    public DistributionResult cdfRespT(AvgHandle R) {
        return getCdfRespT(R);
    }

    /**
     * Filters and processes performance metrics based on enabled handles and masks.
     * This internal method handles disabling metrics, applying zero masks, and
     * cleaning up numerical artifacts.
     *
     * @param avgHandle metric handle configuration
     * @param metric    raw metric matrix to filter
     * @param zeroMask  optional mask to zero out specific entries
     * @return filtered and processed metric matrix
     */
    private Matrix filterMetric(AvgHandle avgHandle, Matrix metric, Matrix zeroMask) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix outData = new Matrix(M, K);
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (!avgHandle.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                        && metric != null && !metric.isEmpty()) {
                    double value = metric.get(i, k);
                    if (value < GlobalConstants.FineTol || Double.isNaN(value)) {
                        value = 0; // Round to zero or disable the metric
                    }
                    outData.set(i, k, value);
                } else {
                    outData.set(i, k, Double.NaN); // Indicates that a metric is disabled
                }
            }
        }

        // NaN values indicate that a metric is disabled
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (Double.isNaN(outData.get(i, j))) {
                    outData.set(i, j, 0);
                }
            }
        }

        // Set to zero entries in the mask
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (zeroMask != null && !zeroMask.isEmpty() && (zeroMask.get(i, j) == 1)) {
                    outData.set(i, j, 0);
                }
            }
        }

        // Round to zero numerical perturbations
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (outData.get(i, j) < GlobalConstants.FineTol) {
                    outData.set(i, j, 0);
                }
            }
        }

        // set to zero metrics for classes that are unreachable;
        // skip for fork-join models where the analytical visits calculation
        // doesn't capture the parent class visiting Join and downstream nodes
        // (post-fork class-switch aux classes get measured by JMT but have
        // visits[c][i, aux] == 0 in the analytical chain matrix);
        // skip when no chains are defined (routing not specified);
        // skip for SPN models where Places don't have traditional visits.
        // Mirrors matlab/src/solvers/@NetworkSolver/getAvg.m lines 184-203.
        boolean hasSPN = sn.nodetype.contains(NodeType.Place) || sn.nodetype.contains(NodeType.Transition);
        boolean hasForkJoin = sn.nodetype.contains(NodeType.Fork) && sn.nodetype.contains(NodeType.Join);
        // A Cache node switches classes by item state, which the routing-based
        // visit equations cannot represent: visits downstream of the cache
        // solve to garbage, so trust the simulation there, like fork-join.
        boolean hasCache = sn.nodetype.contains(NodeType.Cache);

        // A chain fed by spawn-on-completion (sn.classspawn, LQN phase-2) gets
        // its jobs by injection at a station rather than by routing, so the
        // routing-based visit equations solve to zero for its classes even
        // though they are served. Trust the simulation there, like fork-join.
        boolean[] spawnFedChain = new boolean[Math.max(sn.nchains, 0)];
        if (sn.classspawn != null && !sn.classspawn.isEmpty()) {
            for (int k2 = 0; k2 < K && k2 < sn.classspawn.getNumRows(); k2++) {
                int sc2 = (int) sn.classspawn.get(k2, 0);
                if (sc2 >= 0 && sc2 < K) {
                    for (int c2 = 0; c2 < sn.nchains; c2++) {
                        if (sn.chains.get(c2, sc2) > 0) {
                            spawnFedChain[c2] = true;
                        }
                    }
                }
            }
        }

        if (sn.nchains > 0 && !hasSPN) {
            for (int k = 0; k < K; k++) {
                Matrix chainCol = sn.chains.getColumn(k);
                if (chainCol.elementMax() > 0) { // Only check if class k belongs to a chain
                    int c = (int) chainCol.find().value();
                    for (int i = 0; i < M; i++) {
                        if (sn.visits.get(c).get(i, k) == 0) {
                            // Fork-join exception: trust simulation when the
                            // metric value is significant. The analytical
                            // visit matrix doesn't capture the post-fork
                            // class-switch dynamics that JMT measures.
                            if ((hasForkJoin || hasCache || spawnFedChain[c]) && metric != null && !metric.isEmpty()
                                    && metric.get(i, k) > GlobalConstants.FineTol) {
                                continue;
                            }
                            outData.set(i, k, 0);
                        }
                    }
                }
            }
        }

        return outData;
    }

    /**
     * Computes and returns average station metrics at steady-state.
     * This is the main method for obtaining steady-state performance metrics.
     *
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    /**
     * True if station index {@code ist} uses a batch Markovian service process
     * (a BMAP assigned as service, i.e. bulk/BMSP service). For such stations
     * the scalar {@code sn.rates} is the per-event service rate, not the per-job
     * capacity {@code rate*E[batch]}, so throughput-based stability heuristics
     * that divide by {@code sn.rates} do not apply.
     */
    private static boolean isBatchServiceStation(NetworkStruct sn, int ist) {
        if (sn == null || sn.procid == null || sn.stations == null
                || ist < 0 || ist >= sn.stations.size()) {
            return false;
        }
        java.util.Map<jline.lang.JobClass, jline.lang.constant.ProcessType> pmap =
                sn.procid.get(sn.stations.get(ist));
        if (pmap == null) {
            return false;
        }
        for (jline.lang.constant.ProcessType pt : pmap.values()) {
            if (pt == jline.lang.constant.ProcessType.BMAP) {
                return true;
            }
        }
        return false;
    }

    public SolverResult getAvg() {

        if (this.avgHandles == null || this.avgHandles.Q == null || this.avgHandles.U == null || this.avgHandles.R == null ||
                this.avgHandles.T == null || this.avgHandles.A == null) {
            reset();
        }
        this.avgHandles = model.getAvgHandles();

        if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
            throw new RuntimeException(
                    "The getAvg method does not support the timespan option, use the getTranAvg method instead.");
        } else {
            this.options.timespan[0] = Inf;
            this.options.timespan[1] = Inf;
        }

        if (!this.hasAvgResults() || !this.options.cache) {
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                line_error(mfilename(new Object() {
                }), "IllegalAccessException upon running runAnalyzer(): " + e.getMessage());
            } catch (ParserConfigurationException e) {
                line_error(mfilename(new Object() {
                }), "ParserConfigurationException upon running runAnalyzer(): " + e.getMessage());
            } catch (IOException e) {
                line_error(mfilename(new Object() {
                }), "IOException upon running runAnalyzer(): " + e.getMessage());
            } catch (RuntimeException e) {
                line_error(mfilename(new Object() {
                }), "RuntimeException upon running runAnalyzer(): " + e.getMessage());
            }
            if (!this.hasAvgResults()) {
                line_error(mfilename(new Object() {
                }), "Line is unable to return results for this model. " +
                        "The solver " + this.getName() + " may not support this model type.");
                return null;
            }
        } // else return cached value

        int M = sn.nstations;
        int K = sn.nclasses;

        int Vrows = sn.visits.get(0).getNumRows();
        int Vcols = sn.visits.get(0).getNumCols();
        int Vcells = sn.visits.size();
        Matrix V = new Matrix(Vrows, Vcols);
        for (int i = 0; i < Vrows; i++) {
            for (int j = 0; j < Vcols; j++) {
                double tmpSum = 0;
                for (int k = 0; k < Vcells; k++) {
                    tmpSum += sn.visits.get(k).get(i, j);
                }
                V.set(i, j, tmpSum);
            }
        }

        Matrix QNclass = new Matrix(0, 0);
        Matrix UNclass = new Matrix(0, 0);
        Matrix RNclass = new Matrix(0, 0);
        Matrix TNclass = new Matrix(0, 0);
        Matrix ANclass = new Matrix(0, 0);
        Matrix WNclass = new Matrix(0, 0);

        if (!this.avgHandles.R.isEmpty() && this.result.RN != null && !this.result.RN.isEmpty()) {
            RNclass = filterMetric(this.avgHandles.R, this.result.RN, null);
        }

        Matrix zeroMask = new Matrix(M, K);
        for (int i = 0; i < M; i++)
            for (int r = 0; r < K; r++)
                if (RNclass.get(i, r) < 10 * GlobalConstants.FineTol)
                    zeroMask.set(i, r, 1);

        // On an SPN, do not zero Q or U from a near-zero response time: a Place holds
        // tokens rather than serving jobs, so it has no response time, yet it still has
        // a token count and (scheduling as INF) a utilization equal to that count. A
        // resident token whose class is never consumed has RespT = 0 with QLen = 1 and
        // Util = 1. Exempting only Q left such a class reporting QLen = 1 with Util = 0.
        // Mirrors matlab/src/solvers/@NetworkSolver/getAvg.m, which disables both masks
        // when hasSPN = any(nodetype == Place) || any(nodetype == Transition).
        boolean hasSPN = false;
        for (int ind = 0; ind < sn.nnodes && !hasSPN; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt == NodeType.Place || nt == NodeType.Transition) {
                hasSPN = true;
            }
        }

        if (!this.avgHandles.Q.isEmpty() && this.result.QN != null && !this.result.QN.isEmpty()) {
            QNclass = filterMetric(this.avgHandles.Q, this.result.QN, hasSPN ? null : zeroMask);
        }

        if (!this.avgHandles.U.isEmpty() && this.result.UN != null && !this.result.UN.isEmpty()) {
            UNclass = filterMetric(this.avgHandles.U, this.result.UN, hasSPN ? null : zeroMask);
        }

        if (!this.avgHandles.T.isEmpty() && this.result.TN != null && !this.result.TN.isEmpty()) {
            TNclass = filterMetric(this.avgHandles.T, this.result.TN, null);
        }

        if (!this.avgHandles.A.isEmpty()) {
            // If result.AN is available, use it; otherwise compute from throughputs
            if (this.result.AN != null && !this.result.AN.isEmpty()) {
                Matrix zeroMaskSource = Matrix.createLike(zeroMask);
                for (int ist = 0; ist < M; ist++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            zeroMaskSource.set(ist, r, 1);
                        }
                    }
                }
                ANclass = filterMetric(this.avgHandles.A, this.result.AN, zeroMaskSource);
            } else if (!TNclass.isEmpty()) {
                // Compute arrival rates from throughputs when result.AN is not available
                Matrix computedAN = snGetArvRFromTput(sn, TNclass, this.avgHandles.T);
                Matrix zeroMaskSource = Matrix.createLike(zeroMask);
                for (int ist = 0; ist < M; ist++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            zeroMaskSource.set(ist, r, 1);
                        }
                    }
                }
                ANclass = filterMetric(this.avgHandles.A, computedAN, zeroMaskSource);
            }
        }

        if (!this.avgHandles.W.isEmpty()) {
            WNclass = snGetResidTFromRespT(sn, RNclass, this.avgHandles.W);
        }

        if (!UNclass.isEmpty()) {
            // A finite-server open queueing station is unstable when its offered
            // load rho = sum_r T[i,r]/(nservers[i]*rate[i,r]) >= 1. Such a station
            // is fully saturated, so its utilization is reported capped at 1.0
            // (split across classes by offered load) with a single instability
            // warning. rho is recomputed from throughput and service rate so the
            // cap is independent of whatever value the algorithm left in UN (some
            // leave 0). Source and infinite-server/delay stations are excluded.
            boolean infJobsFlag = false;
            for (int i = 0; i < sn.njobs.length(); i++) {
                if (isInf(sn.njobs.get(0, i))) {
                    infJobsFlag = true;
                    break;
                }
            }

            boolean anyUnstable = false;
            if (infJobsFlag && !TNclass.isEmpty() && sn.rates != null && !sn.rates.isEmpty()) {
                for (int i = 0; i < M; i++) {
                    double c = sn.nservers.get(i, 0);
                    if (Double.isInfinite(c) || c <= 0) {
                        continue; // infinite-server / delay station: never capped
                    }
                    int nodeIdx = (int) sn.stationToNode.get(i);
                    if (nodeIdx >= 0 && nodeIdx < sn.nodetype.size()
                            && sn.nodetype.get(nodeIdx) == NodeType.Source) {
                        continue; // source station
                    }
                    // Batch (bulk) service station: sn.rates holds the per-event
                    // service rate, not the per-job capacity rate*E[batch], so the
                    // rho = T/(c*rate) test spuriously reads >= 1. The measured
                    // queue length/response time are trusted instead of capping.
                    if (isBatchServiceStation(sn, i)) {
                        continue;
                    }
                    double[] rho = new double[K];
                    double rhoOpen = 0.0;
                    double rhoTot = 0.0;
                    boolean hasOpen = false;
                    for (int r = 0; r < K; r++) {
                        double rate = sn.rates.get(i, r);
                        double t = TNclass.get(i, r);
                        if (rate > 0 && t > 0) {
                            rho[r] = t / (c * rate);
                            rhoTot += rho[r];
                            if (isInf(sn.njobs.get(0, r))) {
                                rhoOpen += rho[r];
                                hasOpen = true;
                            }
                        }
                    }
                    if (hasOpen && rhoOpen >= 1.0 && rhoTot > 0) {
                        anyUnstable = true;
                        for (int r = 0; r < K; r++) {
                            UNclass.set(i, r, rho[r] / rhoTot); // station total capped to 1.0
                            // A saturated open station has unbounded backlog:
                            // report queue length and response time as Inf for
                            // its open classes (matching native-Python), so SLA
                            // constraints treat it as infeasible rather than
                            // reading a spuriously finite response time.
                            if (isInf(sn.njobs.get(0, r))) {
                                if (RNclass != null && !RNclass.isEmpty()) {
                                    RNclass.set(i, r, Double.POSITIVE_INFINITY);
                                }
                                if (QNclass != null && !QNclass.isEmpty()) {
                                    QNclass.set(i, r, Double.POSITIVE_INFINITY);
                                }
                            }
                        }
                    }
                }
            }

            if (anyUnstable) {
                line_warning(mfilename(new Object() {
                        }),
                        "The model has unstable queues (utilization >= 1); station utilization is reported capped at 1.0, queue length and response time as Inf.");
            }
        }
        double runtime = this.result.runtime;
        // Note: this.result.reset() is intentionally not called here as it would erase
        // transient measures (QNt, UNt, TNt, t, pi_t) that may be needed for analysis
        this.result.QN = QNclass;
        this.result.UN = UNclass;
        this.result.RN = RNclass;
        this.result.AN = ANclass;
        this.result.TN = TNclass;
        this.result.WN = WNclass;
        this.result.runtime = runtime;
        return this.result;
    }

    /**
     * Computes and returns average station metrics at steady-state using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    public SolverResult getAvg(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvg();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Computes and returns average station metrics at steady-state using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    public SolverResult getAvg(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvg(customHandles);
    }

    /**
     * Computes and returns average arrival rates at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average arrival rates [stations x classes]
     */
    public Matrix getAvgArvR() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.AN;
    }

    /**
     * Returns average arrival rates aggregated by job chains.
     *
     * @return matrix of arrival rates [stations x chains]
     */
    public Matrix getAvgArvRChain() {
        int C = sn.nchains;
        Matrix ANclass = getAvgArvR();
        Matrix AN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    AN.set(i, c, AN.get(i, c) + ANclass.get(i, k));
                }
            }
        }
        return AN;
    }

    /**
     * Returns the average arrival rate metric handles.
     *
     * @return handles for arrival rate metrics
     */
    public AvgHandle getAvgArvRHandles() {
        return this.avgHandles.A;
    }

    /**
     * Returns average station metrics aggregated by job chains.
     *
     * @return solver result with metrics aggregated by chains
     */
    public SolverResult getAvgChain() {
        SolverResult res = new SolverResult();
        res.QN = this.getAvgQLenChain();
        res.UN = this.getAvgUtilChain();
        res.RN = this.getAvgRespTChain();
        res.TN = this.getAvgTputChain();
        res.AN = this.getAvgArvRChain();
        res.WN = this.getAvgResidTChain();
        // Preserve the method field from the original result
        if (this.result != null) {
            res.method = this.result.method;
            res.solver = this.result.solver;
            res.runtime = this.result.runtime;
            res.iter = this.result.iter;
        }
        // Note: Do NOT overwrite this.result here as it would replace class-level
        // results (M x K) with chain-level results (M x C), breaking subsequent
        // calls to getAvgSys() which expect class-level dimensions.
        return res;
    }

    /**
     * Returns a table of average station metrics aggregated by job chains.
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int C = sn.nchains;

        SolverResult chainResult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                chainResult = getAvgChain();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgChainTable: " + e.getMessage());
            return null;
        }

        // Use the chain-aggregated results, not this.result which has class-level dimensions
        Matrix QN = chainResult != null ? chainResult.QN : null;
        Matrix UN = chainResult != null ? chainResult.UN : null;
        Matrix RN = chainResult != null ? chainResult.RN : null;
        Matrix TN = chainResult != null ? chainResult.TN : null;
        Matrix AN = chainResult != null ? chainResult.AN : null;
        Matrix WN = chainResult != null ? chainResult.WN : null;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgChainTable.");
            return null;
        }

        if (!keepDisabled) {

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // if (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) > 0) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
            NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgChainTable.setOptions(this.options);
            avgChainTable.setStationNames(stationName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgChainTable.setChainNames(chainNames);
            avgChainTable.setInChainNames(inChainNames);

            return avgChainTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
            NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgChainTable.setOptions(this.options);
            avgChainTable.setStationNames(stationName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgChainTable.setChainNames(chainNames);
            avgChainTable.setInChainNames(inChainNames);

            return avgChainTable;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgChainTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgChainTable(customHandles);
    }

    /**
     * Returns a table of average station metrics aggregated by job chains with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(boolean keepDisabled) {
        // Modify the implementation to support keepDisabled
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int C = sn.nchains;

        SolverResult chainResult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                chainResult = getAvgChain();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgChainTable: " + e.getMessage());
            return null;
        }

        // Use the chain-aggregated results, not this.result which has class-level dimensions
        Matrix QN = chainResult != null ? chainResult.QN : null;
        Matrix UN = chainResult != null ? chainResult.UN : null;
        Matrix RN = chainResult != null ? chainResult.RN : null;
        Matrix TN = chainResult != null ? chainResult.TN : null;
        Matrix AN = chainResult != null ? chainResult.AN : null;
        Matrix WN = chainResult != null ? chainResult.WN : null;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgChainTable.");
            return null;
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) + AN.get(i, c) + WN.get(i, c) > 0)) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
        }
        NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgChainTable.setOptions(this.options);
        avgChainTable.setStationNames(stationName);

        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < this.sn.nchains; c++) {
                chainNames.add("Chain" + (c + 1));
                Matrix inchain = sn.inchain.get(c);
                String chainMembers = "(";
                for (int j = 0; j < inchain.length(); j++) {
                    int r = (int) inchain.get(j);
                    if (j == 0) {
                        chainMembers = chainMembers.concat(sn.classnames.get(r));
                    } else {
                        chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                    }
                }
                inChainNames.add(chainMembers + ")");
            }
        }
        avgChainTable.setChainNames(chainNames);
        avgChainTable.setInChainNames(inChainNames);

        return avgChainTable;
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgChainTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgChainTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aCT aliases for getAvgChainTable
    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns the steady-state average performance metric handles.
     *
     * @return the average handles object
     */
    public SolverAvgHandles getAvgHandles() {
        return this.avgHandles;
    }

    /**
     * Sets the steady-state average performance metric handles.
     *
     * @param handles the average handles to set
     */
    public void setAvgHandles(SolverAvgHandles handles) {
        this.avgHandles = handles;
    }


    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) updateModel() - model is public and therefore no need for setter
    // b) all methods that are "not supported by this solver" - lack of existing method will suffice
    // rather than dedicated warning that method is not applicable

    /**
     * Computes average performance metrics at steady-state for all nodes.
     * This method aggregates station-level metrics to node-level metrics.
     *
     * @return solver result containing node-level average metrics
     */
    public SolverResult getAvgNode() {
        if (this.avgHandles == null || this.avgHandles.Q == null || this.avgHandles.U == null || this.avgHandles.R == null ||
                this.avgHandles.T == null || this.avgHandles.A == null) {
            reset();
        }
        this.avgHandles = model.getAvgHandles();
        SolverResult result = getAvg();

        SolverResult noderesult = new SolverResult();
        int I = sn.nnodes;  // Physical nodes only
        int M = sn.nstations;
        int R = sn.nclasses;
        int C = sn.nchains;
        int F = sn.nregions;  // Number of FCR virtual nodes

        // Total nodes includes physical nodes + FCR virtual nodes
        int totalNodes = I + F;

        Matrix QNn = Matrix.zeros(totalNodes, R);
        Matrix UNn = Matrix.zeros(totalNodes, R);
        Matrix RNn = Matrix.zeros(totalNodes, R);
        Matrix WNn = Matrix.zeros(totalNodes, R);
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int r = 0; r < R; r++) {
                QNn.set(ind, r, result.QN.get(ist, r));
                UNn.set(ind, r, result.UN.get(ist, r));
                RNn.set(ind, r, result.RN.get(ist, r));
                WNn.set(ind, r, result.WN.get(ist, r));
            }
        }

        AvgHandle T = getAvgTputHandles();
        Matrix ANn = snGetNodeArvRFromTput(sn, result.TN, T, result.AN);
        Matrix TNn = snGetNodeTputFromTput(sn, result.TN, T, ANn);

        // Fix arrival rates for ClassSwitch and Sink nodes for cache hit/miss
        // classes: the arrival rate at these nodes for a hit/miss class equals the
        // Cache throughput of that class.
        for (int cacheInd = 0; cacheInd < I; cacheInd++) {
            if (sn.nodetype.get(cacheInd) == NodeType.Cache) {
                CacheNodeParam cnp = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(cacheInd));
                for (int ind = 0; ind < I; ind++) {
                    if (sn.nodetype.get(ind) == NodeType.ClassSwitch || sn.nodetype.get(ind) == NodeType.Sink) {
                        for (int classIdx = 0; classIdx < R; classIdx++) {
                            // A class is a hit/miss class if some jobin class maps to it.
                            // Only positive entries are real assignments: the unset
                            // sentinel is 0 for hitclass and -1 for missclass, and the
                            // read class (index 0) is never a hit/miss class.
                            boolean isHit = false, isMiss = false;
                            if (cnp.hitclass != null) {
                                for (int oc = 0; oc < cnp.hitclass.length(); oc++) {
                                    if (cnp.hitclass.get(oc) > 0 && cnp.hitclass.get(oc) == classIdx) { isHit = true; break; }
                                }
                            }
                            if (cnp.missclass != null) {
                                for (int oc = 0; oc < cnp.missclass.length(); oc++) {
                                    if (cnp.missclass.get(oc) > 0 && cnp.missclass.get(oc) == classIdx) { isMiss = true; break; }
                                }
                            }
                            if (isHit || isMiss) {
                                ANn.set(ind, classIdx, TNn.get(cacheInd, classIdx));
                            }
                        }
                    }
                }
            }
        }

        // Extend ANn and TNn to include FCR rows
        if (F > 0) {
            Matrix ANnExtended = Matrix.zeros(totalNodes, R);
            Matrix TNnExtended = Matrix.zeros(totalNodes, R);
            for (int i = 0; i < I; i++) {
                for (int r = 0; r < R; r++) {
                    ANnExtended.set(i, r, ANn.get(i, r));
                    TNnExtended.set(i, r, TNn.get(i, r));
                }
            }
            ANn = ANnExtended;
            TNn = TNnExtended;
        }

        // Merge FCR metrics into node matrices at FCR node indices
        if (result.QNfcr != null && F > 0) {
            for (int f = 0; f < F; f++) {
                int fcrNodeIdx = I + f;  // FCR indices start after physical nodes
                for (int r = 0; r < R; r++) {
                    QNn.set(fcrNodeIdx, r, result.QNfcr.get(f, r));
                    UNn.set(fcrNodeIdx, r, result.UNfcr.get(f, r));
                    RNn.set(fcrNodeIdx, r, result.RNfcr.get(f, r));
                    WNn.set(fcrNodeIdx, r, result.WNfcr.get(f, r));
                    TNn.set(fcrNodeIdx, r, result.TNfcr.get(f, r));
                    ANn.set(fcrNodeIdx, r, result.ANfcr.get(f, r));
                }
            }
        }

        noderesult.QN = QNn;
        noderesult.UN = UNn;
        noderesult.RN = RNn;
        noderesult.TN = TNn;
        noderesult.AN = ANn;
        noderesult.WN = WNn;
        return noderesult;
    }

    // ========== State Probability Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns average node arrival rates aggregated by job chains.
     *
     * @return matrix of node arrival rates [nodes x chains]
     */
    public Matrix getAvgNodeArvRChain() {
        int C = sn.nchains;
        Matrix ANclass = getAvgNode().AN;
        Matrix AN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    AN.set(i, c, AN.get(i, c) + ANclass.get(i, k));
                }
            }
        }
        return AN;
    }

    /**
     * Returns average node metrics aggregated by job chains.
     *
     * @return solver result with node metrics aggregated by chains
     */
    public SolverResult getAvgNodeChain() {
        SolverResult res = new SolverResult();
        res.QN = this.getAvgNodeQLenChain();
        res.UN = this.getAvgNodeUtilChain();
        res.RN = this.getAvgNodeRespTChain();
        res.TN = this.getAvgNodeTputChain();
        res.AN = this.getAvgNodeArvRChain();
        res.WN = this.getAvgNodeResidTChain();
        this.result = res;
        return res;
    }

    /**
     * Returns a table of average node metrics aggregated by chains.
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvgNodeChain();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgNodeChainTable: " + e.getMessage());
            return null;
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgNodeChainTable.");
            return null;
        }

        if (!keepDisabled) {

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();

            for (int i = 0; i < I; i++) {
                for (int c = 0; c < C; c++) {
                    //if (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) > 0) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                    //}
                }
            }
            NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgNodeChainTable.setOptions(this.options);
            avgNodeChainTable.setNodeNames(nodeName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgNodeChainTable.setChainNames(chainNames);
            avgNodeChainTable.setInChainNames(inChainNames);

            return avgNodeChainTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                }
            }
            NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgNodeChainTable.setOptions(this.options);
            avgNodeChainTable.setNodeNames(nodeName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgNodeChainTable.setChainNames(chainNames);
            avgNodeChainTable.setInChainNames(inChainNames);

            return avgNodeChainTable;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeChainTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeChainTable(customHandles);
    }

    /**
     * Returns a table of average node metrics aggregated by chains with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nnodes;
        int C = sn.nchains;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvgNodeChain();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgNodeChainTable: " + e.getMessage());
            return null;
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgNodeChainTable.");
            return null;
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> nodeName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) + AN.get(i, c) + WN.get(i, c) > 0)) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                }
            }
        }
        NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgNodeChainTable.setOptions(this.options);
        avgNodeChainTable.setNodeNames(nodeName);

        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < this.sn.nchains; c++) {
                chainNames.add("Chain" + (c + 1));
                Matrix inchain = sn.inchain.get(c);
                String chainMembers = "(";
                for (int j = 0; j < inchain.length(); j++) {
                    int r = (int) inchain.get(j);
                    if (j == 0) {
                        chainMembers = chainMembers.concat(sn.classnames.get(r));
                    } else {
                        chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                    }
                }
                inChainNames.add(chainMembers + ")");
            }
        }
        avgNodeChainTable.setChainNames(chainNames);
        avgNodeChainTable.setInChainNames(inChainNames);

        return avgNodeChainTable;
    }

    /**
     * Returns a table of average node metrics aggregated by chains using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeChainTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeChainTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT() {
        return getAvgNodeChainTable();
    }

    /**
     * Alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT() {
        return getAvgNodeChainTable();
    }

    /**
     * Alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT() {
        return getAvgNodeChainTable();
    }

    /**
     * Short alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Short alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Short alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns average node queue lengths aggregated by job chains.
     *
     * @return matrix of node queue lengths [nodes x chains]
     */
    public Matrix getAvgNodeQLenChain() {
        int C = sn.nchains;
        Matrix QNclass = getAvgNode().QN;
        Matrix QN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    QN.set(i, c, QN.get(i, c) + QNclass.get(i, k));
                }
            }
        }
        return QN;
    }

    // ========== Sampling Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns average node residence times aggregated by job chains.
     *
     * @return matrix of node residence times [nodes x chains]
     */
    public Matrix getAvgNodeResidTChain() {
        int C = sn.nchains;
        Matrix WNclass = getAvgNode().WN;
        Matrix WN = new Matrix(sn.nnodes, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    if (sn.isstation.get(i) == 1)
                        WN.set(i, c, WN.get(i, c) + WNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return WN;
    }

    /**
     * Returns average node response times aggregated by job chains.
     *
     * @return matrix of node response times [nodes x chains]
     */
    public Matrix getAvgNodeRespTChain() {
        int C = sn.nchains;
        Matrix RNclass = getAvgNode().RN;
        Matrix RN = new Matrix(sn.nnodes, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    if (sn.isstation.get(i) == 1)
                        RN.set(i, c, RN.get(i, c) + RNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return RN;
    }

    /**
     * Returns a table of average node metrics organized by job classes.
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        // Total nodes includes physical nodes + FCR virtual nodes
        int I = sn.nnodes + sn.nregions;
        int K = sn.nclasses;

        // Auxiliary retrieval classes (created internally by Cache.setRetrievalSystem
        // to model the delayed-hit retrieval system) are hidden from the node table.
        boolean[] hiddenClass = new boolean[K];
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                CacheNodeParam cnp = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(i));
                if (cnp != null && cnp.retrievalClassIndices != null) {
                    for (Integer rci : cnp.retrievalClassIndices) {
                        if (rci != null && rci >= 0 && rci < K) hiddenClass[rci] = true;
                    }
                }
            }
        }

        SolverResult noderesult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                noderesult = getAvgNode();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgNodeTable: " + e.getMessage());
            return null;
        }

        // Check if noderesult is null before accessing its fields
        if (noderesult == null) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results - solver execution failed or did not complete.");
            return null;
        }

        Matrix QN = noderesult.QN;
        Matrix UN = noderesult.UN;
        Matrix RN = noderesult.RN;
        Matrix WN = noderesult.WN;
        Matrix TN = noderesult.TN;
        Matrix AN = noderesult.AN;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgNodeTable.");
            return null;
        }

        if (!keepDisabled) {
            Matrix V = new Matrix(sn.nnodes, K);
            for (int i = 0; i < sn.nodevisits.size(); i++) {
                V = V.add(1, sn.nodevisits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.nodevisits.keySet()) {
                    Matrix chainVisits = sn.nodevisits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < K; k++) {
                    if (hiddenClass[k]) {
                        continue; // auxiliary retrieval class - omit from node table
                    }
                    int c = -1;
                    for (int row = 0; row < sn.chains.getNumRows(); row++) {
                        if (sn.chains.get(row, k) > 0) {
                            c = row;
                            break;
                        }
                    }
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setNodeNames(nodeName);

            return avgTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            Matrix V = new Matrix(sn.nnodes, K);
            for (int i = 0; i < sn.nodevisits.size(); i++) {
                V = V.add(1, sn.nodevisits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.nodevisits.keySet()) {
                    Matrix chainVisits = sn.nodevisits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < K; k++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setNodeNames(nodeName);

            return avgTable;
        }
    }

    /**
     * Returns a table of average node metrics organized by job classes using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    // ========== Distribution Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns a table of average node metrics organized by job classes using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeTable(customHandles);
    }

    /**
     * Returns a table of average node metrics organized by job classes with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        // Total nodes includes physical nodes + FCR virtual nodes
        int I = sn.nnodes + sn.nregions;
        int K = sn.nclasses;

        SolverResult noderesult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                noderesult = getAvgNode();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results for AvgNodeTable: " + e.getMessage());
            return null;
        }

        if (noderesult == null) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results - solver execution failed or did not complete.");
            return null;
        }

        Matrix QN = noderesult.QN;
        Matrix UN = noderesult.UN;
        Matrix RN = noderesult.RN;
        Matrix WN = noderesult.WN;
        Matrix TN = noderesult.TN;
        Matrix AN = noderesult.AN;

        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgNodeTable.");
            return null;
        }

        Matrix V = new Matrix(sn.nnodes, K);
        for (int i = 0; i < sn.nodevisits.size(); i++) {
            V = V.add(1, sn.nodevisits.get(i));
        }
        if (V.isEmpty()) { // SSA
            // Implementation for SSA: sum all chain visits across all chains
            for (Integer chainIndex : sn.nodevisits.keySet()) {
                Matrix chainVisits = sn.nodevisits.get(chainIndex);
                if (chainVisits != null) {
                    V = V.add(1.0, chainVisits);
                }
            }
        }

        // Auxiliary retrieval classes (created internally by Cache.setRetrievalSystem
        // to model the delayed-hit retrieval system) are hidden from the node table:
        // they are plumbing, not user-facing classes.
        boolean[] hiddenClass = new boolean[K];
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                CacheNodeParam cnp = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(i));
                if (cnp != null && cnp.retrievalClassIndices != null) {
                    for (Integer rci : cnp.retrievalClassIndices) {
                        if (rci != null && rci >= 0 && rci < K) hiddenClass[rci] = true;
                    }
                }
            }
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> nodeName = new ArrayList<>();
        for (int i = 0; i < I; i++) {
            for (int k = 0; k < K; k++) {
                if (hiddenClass[k]) {
                    continue; // auxiliary retrieval class - omit from node table
                }
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0)) {
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
        }
        NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgTable.setOptions(this.options);
        avgTable.setClassNames(className);
        avgTable.setNodeNames(nodeName);

        return avgTable;
    }

    /**
     * Returns a table of detailed per-class cache performance metrics for every
     * Cache node: a total row (List=0) per cache+read class plus, where the solver
     * reports per-list hit probabilities and the cache has more than one list, one
     * row per cache list. Port of matlab getAvgCacheTable.m.
     *
     * @return cache performance table
     */
    public NetworkAvgCacheTable getAvgCacheTable() {
        this.sn = model.getStruct(true);
        int K = sn.nclasses;

        List<Double> List_ = new ArrayList<>(), ListCap = new ArrayList<>(), Items = new ArrayList<>();
        List<Double> HitProb = new ArrayList<>(), DelayedHitProb = new ArrayList<>(), MissProb = new ArrayList<>();
        List<Double> HitRate = new ArrayList<>(), DelayedHitRate = new ArrayList<>(), MissRate = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>(), Latency = new ArrayList<>();
        List<String> nodeName = new ArrayList<>(), className = new ArrayList<>();

        boolean anyCache = false;
        for (int i = 0; i < sn.nnodes; i++) if (sn.nodetype.get(i) == NodeType.Cache) { anyCache = true; break; }
        if (anyCache) {
            SolverResult nr = getAvgNode();
            Matrix TNnode = nr.TN;
            // Source node index: the cache read-class arrival equals that class's
            // source throughput (every read request enters the cache).
            int srcNode = -1;
            for (int i = 0; i < sn.nnodes; i++) if (sn.nodetype.get(i) == NodeType.Source) { srcNode = i; break; }
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.nodetype.get(ind) != NodeType.Cache) continue;
                CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                Cache cache = (Cache) model.getNodes().get(ind);
                Matrix hitclass = np.hitclass;
                Matrix missclass = np.missclass;
                Matrix itemcap = np.itemcap;
                int h = (itemcap == null) ? 0 : itemcap.length();
                double totcap = 0; for (int l = 0; l < h; l++) totcap += itemcap.get(l);
                int nitems = np.nitems;
                Matrix hitp = cache.getHitRatio(), missp = cache.getMissRatio();
                Matrix dhitp = cache.getDelayedHitRatio(), lat = cache.getResidT();
                Matrix hitplist = cache.getHitRatioByList();
                for (int r = 0; r < K; r++) {
                    if (hitclass == null || r >= hitclass.length() || hitclass.get(r) <= 0) continue;
                    double ph = nanGetAt(hitp, r), pm = nanGetAt(missp, r), pd = nanGetAt(dhitp, r);
                    if (Double.isNaN(ph) && Double.isNaN(pm) && Double.isNaN(pd)) continue;
                    if (Double.isNaN(ph)) ph = 0; if (Double.isNaN(pm)) pm = 0; if (Double.isNaN(pd)) pd = 0;
                    // Read-class arrival into the cache equals that class's source
                    // throughput (every read request enters the cache). Robust across
                    // solvers, including simulators where delayed hits are not folded
                    // into the cache hit/miss throughput.
                    double arvr = 0;
                    if (TNnode != null && srcNode >= 0 && srcNode < TNnode.getNumRows() && r < TNnode.getNumCols()) {
                        arvr = TNnode.get(srcNode, r);
                    }
                    double latr = nanGetAt(lat, r);
                    // total row (List = 0). ArvR is the retrieval-system throughput
                    // arvr*(missprob+delayedprob) = MissRate + DelayedHitRate, i.e. the rate
                    // of requests that enter the retrieval system (a miss or a delayed hit).
                    // This is the arrival rate that is Little-consistent with ResidT = Z (the
                    // delayed-hit expected latency, eq:latency tot / retrieval_fpi_latency):
                    // ArvR*ResidT = sum_i(phi_i+d_i), the mean number of requests in the
                    // retrieval system (fetch job included).
                    double arvrRetr = arvr * (pm + pd);
                    nodeName.add(sn.nodenames.get(ind)); className.add(model.getClasses().get(r).getName());
                    List_.add(0.0); ListCap.add(totcap); Items.add((double) nitems);
                    HitProb.add(ph); DelayedHitProb.add(pd); MissProb.add(pm);
                    HitRate.add(arvr * ph); DelayedHitRate.add(arvr * pd); MissRate.add(arvr * pm);
                    ArvR.add(arvrRetr); Latency.add(latr);
                    // per-list rows (only when a multi-list breakdown is available)
                    boolean haveList = h > 1 && hitplist != null && r < hitplist.getNumRows();
                    if (haveList) {
                        boolean allNaN = true;
                        for (int l = 0; l < hitplist.getNumCols(); l++) if (!Double.isNaN(hitplist.get(r, l))) { allNaN = false; break; }
                        if (!allNaN) {
                            for (int l = 0; l < h; l++) {
                                double phl = (l < hitplist.getNumCols()) ? hitplist.get(r, l) : Double.NaN;
                                if (Double.isNaN(phl)) phl = 0;
                                double capl = (l < itemcap.length()) ? itemcap.get(l) : Double.NaN;
                                nodeName.add(sn.nodenames.get(ind)); className.add(model.getClasses().get(r).getName());
                                List_.add((double) (l + 1)); ListCap.add(capl); Items.add((double) nitems);
                                HitProb.add(phl); DelayedHitProb.add(Double.NaN); MissProb.add(Double.NaN);
                                HitRate.add(arvr * phl); DelayedHitRate.add(Double.NaN); MissRate.add(Double.NaN);
                                ArvR.add(arvr); Latency.add(Double.NaN);
                            }
                        }
                    }
                }
            }
        }

        NetworkAvgCacheTable t = new NetworkAvgCacheTable(List_, ListCap, Items, HitProb, DelayedHitProb,
                MissProb, HitRate, DelayedHitRate, MissRate, ArvR, Latency);
        t.setOptions(this.options);
        t.setNodeNames(nodeName);
        t.setClassNames(className);
        return t;
    }

    private static double nanGetAt(Matrix v, int r) {
        if (v == null || v.isEmpty() || r >= v.length()) return Double.NaN;
        return v.get(r);
    }

    /**
     * Returns a table of item-level cache occupancy: one row per Cache node,
     * item and cache list (level), with the steady-state probability the item
     * resides in that list. Populated only where the solver computes a per-item
     * distribution (exact cache algorithms and the delayed-hit retrieval algorithms).
     * Port of matlab getAvgItemTable.m.
     *
     * @return item-level cache occupancy table
     */
    public NetworkAvgItemTable getAvgItemTable() {
        this.sn = model.getStruct(true);

        List<Double> Item = new ArrayList<>(), List_ = new ArrayList<>();
        List<Double> ListCap = new ArrayList<>(), Prob = new ArrayList<>();
        List<String> nodeName = new ArrayList<>();

        boolean anyCache = false;
        for (int i = 0; i < sn.nnodes; i++) if (sn.nodetype.get(i) == NodeType.Cache) { anyCache = true; break; }
        if (anyCache) {
            getAvgNode(); // ensure solved so cache item probabilities are filled
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.nodetype.get(ind) != NodeType.Cache) continue;
                CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                Cache cache = (Cache) model.getNodes().get(ind);
                Matrix itemcap = np.itemcap;
                int h = (itemcap == null) ? 0 : itemcap.length();
                Matrix itemprob = cache.getItemProb();
                if (itemprob == null || itemprob.isEmpty() || h == 0) continue;
                int n = itemprob.getNumRows();
                for (int i = 0; i < n; i++) {
                    for (int l = 0; l < h; l++) {
                        double p = (l + 1 < itemprob.getNumCols()) ? itemprob.get(i, l + 1) : Double.NaN;
                        nodeName.add(sn.nodenames.get(ind));
                        Item.add((double) (i + 1));
                        List_.add((double) (l + 1));
                        ListCap.add(itemcap.get(l));
                        Prob.add(p);
                    }
                }
            }
        }

        NetworkAvgItemTable t = new NetworkAvgItemTable(Item, List_, ListCap, Prob);
        t.setOptions(this.options);
        t.setNodeNames(nodeName);
        return t;
    }

    /**
     * Returns a table of average node metrics organized by job classes using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics organized by job classes using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aNT aliases for getAvgNodeTable
    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns average node throughputs aggregated by job chains.
     *
     * @return matrix of node throughputs [nodes x chains]
     */
    public Matrix getAvgNodeTputChain() {
        int C = sn.nchains;
        Matrix TNclass = getAvgNode().TN;
        Matrix TN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    TN.set(i, c, TN.get(i, c) + TNclass.get(i, k));
                }
            }
        }
        return TN;
    }

    /**
     * Returns average node utilizations aggregated by job chains.
     *
     * @return matrix of node utilizations [nodes x chains]
     */
    public Matrix getAvgNodeUtilChain() {
        int C = sn.nchains;
        Matrix UNclass = getAvgNode().UN;
        Matrix UN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    UN.set(i, c, UN.get(i, c) + UNclass.get(i, k));
                }
            }
        }
        return UN;
    }

    /**
     * Computes and returns average queue lengths at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average queue lengths [stations x classes]
     */
    public Matrix getAvgQLen() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.QN;
    }

    /**
     * Mean number of jobs waiting in the ORBIT of each retrial station, as an
     * (nstations x nclasses) matrix. Stations that are not retrial queues report 0.
     *
     * <p>A retrial station has no waiting room: a job that finds every server busy
     * joins the orbit instead of queueing, so its station population splits into the
     * jobs currently in service and the jobs orbiting. getAvgQLen reports the whole
     * station population, which is why the orbit had to be recovered by hand as
     * QLen - Util. This method reports it directly.</p>
     *
     * <p>The in-service population is obtained from the station throughput by
     * Little's law applied to the servers alone, E[in service] = X * E[S], which
     * holds for any service distribution and any number of servers, so the orbit
     * length is exact whenever QLen and Tput are.</p>
     *
     * @return the mean orbit length per station and class
     */
    public Matrix getAvgOrbit() {
        Matrix QN = getAvgQLen();
        Matrix TN = getAvgTput();
        NetworkStruct snl = this.model.getStruct(false);
        Matrix ON = new Matrix(QN.getNumRows(), QN.getNumCols());
        ON.zero();
        if (snl.retrialProc == null) {
            return ON;
        }
        for (int ist = 0; ist < ON.getNumRows(); ist++) {
            if (ist >= snl.stations.size()) {
                continue;
            }
            Map<JobClass, MatrixCell> procs = snl.retrialProc.get(snl.stations.get(ist));
            if (procs == null) {
                continue;
            }
            for (int r = 0; r < ON.getNumCols(); r++) {
                MatrixCell proc = procs.get(snl.jobclasses.get(r));
                if (proc == null) {
                    continue; // not a retrial station for this class: no orbit
                }
                double rate = snl.rates.get(ist, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    continue;
                }
                double inService = TN.get(ist, r) / rate; // Little's law on the servers
                ON.set(ist, r, Math.max(0.0, QN.get(ist, r) - inService));
            }
        }
        return ON;
    }

    /**
     * Table of the mean orbit length of every retrial station-class pair, with the
     * station population and the in-service population it decomposes into.
     *
     * <p>Reported as a separate table rather than as an extra column of getAvgTable
     * so that the average table keeps its shape for models without retrials.</p>
     *
     * @return the orbit table
     */
    public NetworkAvgOrbitTable getAvgOrbitTable() {
        Matrix ON = getAvgOrbit();
        Matrix QN = getAvgQLen();
        Matrix TN = getAvgTput();
        NetworkStruct snl = this.model.getStruct(false);

        List<String> stationNames = new ArrayList<String>();
        List<String> classNames = new ArrayList<String>();
        List<Double> qlen = new ArrayList<Double>();
        List<Double> inservice = new ArrayList<Double>();
        List<Double> orbit = new ArrayList<Double>();

        if (snl.retrialProc != null) {
            for (int ist = 0; ist < ON.getNumRows(); ist++) {
                if (ist >= snl.stations.size()) {
                    continue;
                }
                Map<JobClass, MatrixCell> procs = snl.retrialProc.get(snl.stations.get(ist));
                if (procs == null) {
                    continue;
                }
                for (int r = 0; r < ON.getNumCols(); r++) {
                    if (procs.get(snl.jobclasses.get(r)) == null) {
                        continue;
                    }
                    stationNames.add(snl.stations.get(ist).getName());
                    classNames.add(snl.jobclasses.get(r).getName());
                    qlen.add(QN.get(ist, r));
                    double rate = snl.rates.get(ist, r);
                    if (Double.isFinite(rate) && rate > 0) {
                        inservice.add(TN.get(ist, r) / rate);
                    } else {
                        inservice.add(0.0);
                    }
                    orbit.add(ON.get(ist, r));
                }
            }
        }

        NetworkAvgOrbitTable table = new NetworkAvgOrbitTable(qlen, inservice, orbit);
        table.setStationNames(stationNames);
        table.setClassNames(classNames);
        table.setOptions(this.options);
        return table;
    }

    /**
     * Returns a table of loss (drop) metrics for every station-class pair that
     * receives offered traffic: offered arrival rate, carried throughput, loss
     * rate (ArvR - Tput) and loss ratio (LossRate / ArvR). Only pairs with
     * ArvR &gt; 0 are listed, which excludes the Source. Port of MATLAB
     * getAvgLossTable.m.
     *
     * @return the loss table
     */
    public NetworkLossTable getAvgLossTable() {
        Matrix AN = getAvgArvR();
        Matrix TN = getAvgTput();
        NetworkStruct snl = this.model.getStruct(false);

        List<String> stationNames = new ArrayList<String>();
        List<String> classNames = new ArrayList<String>();
        List<Double> arvr = new ArrayList<Double>();
        List<Double> tput = new ArrayList<Double>();
        List<Double> lossRate = new ArrayList<Double>();
        List<Double> lossRatio = new ArrayList<Double>();

        // Fork-Join quorum sibling-drop rate (LDES only): at a synchronizing Join
        // the station identity LossRate = ArvR - Tput does not hold (Tput is in
        // parent units, discarded siblings in sibling units), so on Join rows the
        // explicit drop rate replaces ArvR - Tput. ArvR is already the offered
        // sibling rate (flow balance over all forked siblings), so LossRatio =
        // drop / ArvR.
        Matrix dropJoin = this.result == null ? null : this.result.DropRateJoin;

        for (int ist = 0; ist < AN.getNumRows(); ist++) {
            if (ist >= snl.stations.size()) {
                continue;
            }
            for (int r = 0; r < AN.getNumCols(); r++) {
                double a = AN.get(ist, r);
                if (!Double.isFinite(a) || a <= 0) {
                    continue;
                }
                double t = TN.get(ist, r);
                double d = (dropJoin != null && ist < dropJoin.getNumRows()
                        && r < dropJoin.getNumCols()) ? dropJoin.get(ist, r) : 0.0;
                double lr;
                double lc;
                if (Double.isFinite(d) && d > 0) {
                    lr = d;
                    lc = d / a;
                } else {
                    lr = a - t;
                    lc = (a - t) / a;
                }
                stationNames.add(snl.stations.get(ist).getName());
                classNames.add(snl.jobclasses.get(r).getName());
                arvr.add(a);
                tput.add(t);
                lossRate.add(lr);
                lossRatio.add(lc);
            }
        }

        NetworkLossTable table = new NetworkLossTable(arvr, tput, lossRate, lossRatio);
        table.setStationNames(stationNames);
        table.setClassNames(classNames);
        table.setOptions(this.options);
        return table;
    }

    /**
     * Returns a table of loss (drop) metrics per finite-capacity region and
     * class, for regions that drop jobs (DROP rule). Each row reports the
     * offered arrival rate (carried Tput plus drop rate), the carried
     * throughput, the loss rate (region drop rate) and the loss ratio
     * (LossRate / ArvR). Only regions with offered traffic are listed; empty
     * for solvers that do not track region drops (only the LDES simulation
     * populates the FCR drop rate).
     *
     * @return the region loss table
     */
    public NetworkLossTable getAvgRegionLossTable() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        List<String> regionNames = new ArrayList<String>();
        List<String> classNames = new ArrayList<String>();
        List<Double> arvr = new ArrayList<Double>();
        List<Double> tput = new ArrayList<Double>();
        List<Double> lossRate = new ArrayList<Double>();
        List<Double> lossRatio = new ArrayList<Double>();

        Matrix TNfcr = this.result == null ? null : this.result.TNfcr;
        Matrix DR = this.result == null ? null : this.result.DropRateNfcr;
        if (TNfcr != null && DR != null) {
            for (int f = 0; f < DR.getNumRows(); f++) {
                for (int r = 0; r < DR.getNumCols(); r++) {
                    double t = TNfcr.get(f, r);
                    double d = DR.get(f, r);
                    double a = t + d;
                    if (!Double.isFinite(a) || a <= 0) {
                        continue;
                    }
                    regionNames.add("FCRegion" + (f + 1));
                    classNames.add(snClassName(r));
                    arvr.add(a);
                    tput.add(t);
                    lossRate.add(d);
                    lossRatio.add(d / a);
                }
            }
        }

        NetworkLossTable table = new NetworkLossTable(arvr, tput, lossRate, lossRatio);
        table.setStationNames(regionNames);
        table.setClassNames(classNames);
        table.setFirstColumnName("Region");
        table.setOptions(this.options);
        return table;
    }

    private String snClassName(int r) {
        NetworkStruct snl = this.model.getStruct(false);
        if (snl.jobclasses != null && r < snl.jobclasses.size()) {
            return snl.jobclasses.get(r).getName();
        }
        return "Class" + (r + 1);
    }

    /**
     * Returns average queue lengths aggregated by job chains.
     *
     * @return matrix of queue lengths [stations x chains]
     */
    public Matrix getAvgQLenChain() {
        int C = sn.nchains;
        Matrix QNclass = getAvgQLen();
        Matrix QN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    QN.set(i, c, QN.get(i, c) + QNclass.get(i, k));
                }
            }
        }
        return QN;
    }

    /**
     * Returns the average queue length metric handles.
     *
     * @return handles for queue length metrics
     */
    public AvgHandle getAvgQLenHandles() {
        return this.avgHandles.Q;
    }

    // ========== Kotlin-style Alias Methods ==========
    // Aliases for get* methods following Kotlin naming conventions

    /**
     * Computes and returns average residence times in queue (including service).
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average residence times [stations x classes]
     */
    public Matrix getAvgResidT() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.WN;
    }

    /**
     * Returns average residence times aggregated by job chains.
     *
     * @return matrix of residence times [stations x chains]
     */
    public Matrix getAvgResidTChain() {
        int C = sn.nchains;
        Matrix WNclass = getAvgResidT();
        Matrix WN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                double sum = 0.0;
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    sum += WNclass.get(i, k);
                }
                WN.set(i, c, sum);
            }
        }
        return WN;
    }

    /**
     * Returns the average residence time metric handles.
     *
     * @return handles for residence time metrics
     */
    public AvgHandle getAvgResidTHandles() {
        return this.avgHandles.W;
    }

    /**
     * Computes and returns average response times at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average response times [stations x classes]
     */
    public Matrix getAvgRespT() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.RN;
    }

    /**
     * Returns average response times aggregated by job chains.
     *
     * @return matrix of response times [stations x chains]
     */
    public Matrix getAvgRespTChain() {
        int C = sn.nchains;
        Matrix RNclass = getAvgRespT();
        Matrix RN = new Matrix(sn.nstations, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    RN.set(i, c, RN.get(i, c) + RNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return RN;
    }

    /**
     * Returns the average response time metric handles.
     *
     * @return handles for response time metrics
     */
    public AvgHandle getAvgRespTHandles() {
        return this.avgHandles.R;
    }

    /**
     * Computes average system-level metrics at steady state.
     * This includes system response times and throughputs aggregated across all chains.
     */
    public void getAvgSys() {

        this.sn = model.getStruct(true);

        this.getAvg();

        // Fork-join support for system response time computation
        if (this.model.hasFork() || this.model.hasJoin()) {
            if (this.model.hasOpenClasses()) {
                line_error(mfilename(new Object() {
                        }),
                        "System response time computation not yet supported with open classes in the presence of fork/join nodes.");
                this.result.RN.fill(Double.NaN);
            } else {
                // For closed networks with fork-join, ensure response times are properly adjusted
                // Fork-join networks require special handling of synchronization delays
                // The join node response time should include synchronization delay
                if (sn.fj != null && !sn.fj.isEmpty()) {
                    // Apply fork-join transformations to response times
                    // This is a simplified implementation - full fork-join support would require
                    // detailed path analysis and synchronization delay computation
                    for (int i = 0; i < sn.nnodes; i++) {
                        if (sn.nodetype.get(i) == NodeType.Join) {
                            // Join nodes get synchronization delay - placeholder implementation
                            // In full implementation, this would be computed using order statistics
                            // based on parallel path response times
                            for (int r = 0; r < sn.nclasses; r++) {
                                if (this.result.RN != null && !this.result.RN.isEmpty()) {
                                    int stationIndex = (int) sn.nodeToStation.get(i);
                                    double currentRT = this.result.RN.get(stationIndex, r);
                                    if (!Double.isNaN(currentRT) && currentRT > 0) {
                                        // Apply basic synchronization delay factor
                                        this.result.RN.set(stationIndex, r, currentRT * 1.2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        boolean[] completes = new boolean[sn.nclasses];
        for (int idx = 0; idx < sn.nclasses; idx++) {
            completes[idx] = model.getClasses().get(idx).getCompletes();
        }

        // Optimization: Check if the model has any open classes (infinite jobs)
        // This could be optimised by computing the statistics only for open chains
        boolean hasOpenClasses = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (isInf(sn.njobs.get(k))) {
                hasOpenClasses = true;
                break;
            }
        }
        // Note: Optimization for open classes not yet implemented - would skip closed chain calculations

        // Compute chain visits
        Matrix alpha = new Matrix(sn.nstations, sn.nclasses);
        Matrix CNclass = new Matrix(1, sn.nclasses);
        if (!this.model.hasJoin() && !this.model.hasFork()) {
            for (int c = 0; c < sn.nchains; c++) {
                Matrix inchain = sn.inchain.get(c);
                for (int i = 0; i < inchain.length(); i++) {
                    int r = (int) inchain.get(i);
                    for (int j = 0; j < sn.nstations; j++) {
                        // Not empty and not a source
                        if (this.result.RN != null && !this.result.RN.isEmpty() &&
                                (!(isInf(this.sn.njobs.get(r)) && j == sn.refstat.get(r)))) {
                            CNclass.set(
                                    0,
                                    r,
                                    CNclass.get(0, r)
                                            + sn.visits.get(c).get((int) sn.stationToStateful.get(j), r)
                                            * this.result.RN.get(j, r)
                                            / sn.visits
                                            .get(c)
                                            .get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), r));
                        }
                    }
                }
            }
        }


        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null);
            for (int i = 0; i < completingClasses.length(); i++) {
                if (!completes[i]) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            for (int i = 0; i < sn.nstations; i++) {
                if (sn.refclass.get(c) >= 0) {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                    for (double value : intersection) {
                        int k = (int) value;
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits +=
                                        sn.visits
                                                .get(c)
                                                .get(
                                                        (int) sn.stationToStateful.get((int) sn.refstat.get(k)),
                                                        idx);
                            }
                        }
                        alpha.set(
                                i,
                                k,
                                alpha.get(i, k)
                                        + sn.visits.get(c).get((int) sn.stationToStateful.get(i), k) / sumVisits);
                    }
                } else {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    for (int j = 0; j < inchain.length(); j++) {
                        int k = (int) inchain.get(j);
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits += sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(k)), idx);
                            }
                        }
                        alpha.set(i, k, alpha.get(i, k) + sn.visits.get(c).get((int) sn.stationToStateful.get(i), k) / sumVisits);
                    }
                }
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (isInf(alpha.get(i, k)) || Double.isNaN(alpha.get(i, k))) {
                    alpha.set(i, k, 0.0);
                }
            }
        }

        // Compute average chain metrics
        this.result.CN = new Matrix(1, sn.nchains);
        this.result.XN = new Matrix(1, sn.nchains);

        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null).find();
            completingClasses = completingClasses.transpose();
            for (int i = 0; i < inchain.length(); i++) {
                int classIndex = (int) inchain.get(i);
                if (!this.model.getClasses().get(classIndex).getCompletes()) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            if (!result.TN.isEmpty()) {
                // All classes in same chain must share the same refstation, so we use the first one
                int ref = (int) sn.refstat.get((int) inchain.get(0));
                // We now compute the incoming system throughput to the reference station from completing
                // classes
                for (int i = 0; i < sn.nstations; i++) {
                    for (int j = 0; j < completingClasses.length(); j++) {
                        int r = (int) completingClasses.get(j);
                        if (completingClasses.get(j) >= 0) {
                            List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                            for (double value : intersection) {
                                int s = (int) value;
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.result.XN.set(
                                            0,
                                            c,
                                            this.result.XN.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                            for (int k = 0; k < inchain.length(); k++) {
                                int s = (int) inchain.get(k);
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.result.XN.set(
                                            0,
                                            c,
                                            this.result.XN.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                        }
                    }
                }
            }

            // If this is a closed chain we simply apply Little's law
            int nJobsChain = 0;
            for (int i = 0; i < sn.chains.getNumCols(); i++) {
                if (sn.chains.get(c, i) > 0) {
                    nJobsChain += sn.njobs.get(i);
                }
            }

            if (this.model.hasFork() && this.model.hasJoin()) {
                // In this case, CN is unreliable as it sums the contribution across all stations,
                // which would include also forked tasks, we use Little's law instead
                this.result.CN.set(0, c, nJobsChain / this.result.XN.get(0, c));
            } else {
                // Standard chain response time computation for non-fork-join networks
                if (isInf(nJobsChain)) {
                    if (inchain.length() != completingClasses.length()) {
                        throw new RuntimeException(
                                "Edge-based chain definition not yet supported for open queueing networks.");
                    }
                }
                double sumFinite = 0;
                for (int i = 0; i < inchain.length(); i++) {
                    double value = alpha.get((int) sn.refstat.get((int) inchain.get(0)), (int) inchain.get(i))
                            * CNclass.get((int) inchain.get(i));
                    if (!isInf(value) && !Double.isNaN(value))
                        sumFinite += value;
                }

                this.result.CN.set(0, c, sumFinite);
            }
        }
    }

    /**
     * Computes average system-level metrics at steady state using specified handles.
     * This includes system response times and throughputs aggregated across all chains.
     *
     * @param avgHandles custom handles for performance metrics
     */
    public void getAvgSys(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            getAvgSys();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Computes average system-level metrics at steady state using individual handles.
     * This includes system response times and throughputs aggregated across all chains.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     */
    public void getAvgSys(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        getAvgSys(customHandles);
    }

    /**
     * Returns average system response times at steady state.
     *
     * @return matrix of system response times by chain
     */
    public Matrix getAvgSysRespT() {
        if (!this.hasResults()) {
            this.getAvgSys();
        }
        return this.result.CN;
    }

    /**
     * Returns a table of average system-level metrics.
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable() {

        this.getAvgSys();

        NetworkAvgSysTable avgSysTable = new NetworkAvgSysTable(this.result.CN.toList1D(), this.result.XN.toList1D(), this.options);


        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int c = 0; c < this.sn.nchains; c++) {
            chainNames.add("Chain" + (c + 1));
            Matrix inchain = sn.inchain.get(c);
            String chainMembers = "(";
            for (int i = 0; i < inchain.length(); i++) {
                int r = (int) inchain.get(i);
                if (i == 0) {
                    chainMembers = chainMembers.concat(sn.classnames.get(r));
                } else {
                    chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                }
            }
            inChainNames.add(chainMembers + ")");
        }

        avgSysTable.setChainNames(chainNames);
        avgSysTable.setInChainNames(inChainNames);

        return avgSysTable;
    }

    /**
     * Returns a table of average system-level metrics using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgSysTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average system-level metrics using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgSysTable(customHandles);
    }

    /**
     * Alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT() {
        return getAvgSysTable();
    }

    /**
     * Alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT() {
        return getAvgSysTable();
    }

    /**
     * Alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Short alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST() {
        return getAvgSysTable();
    }

    /**
     * Short alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Short alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Returns average system throughputs at steady state.
     *
     * @return matrix of system throughputs by chain
     */
    public Matrix getAvgSysTput() {
        if (!this.hasResults()) {
            this.getAvgSys();
        }
        return this.result.XN;
    }

    /**
     * Returns a table of average station metrics organized by job classes.
     *
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (Double.isFinite(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgTable: " + e.getMessage());
            return null;
        }

        // Check if result is null or incomplete before accessing fields
        if (this.result == null) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgTable.");
            return null;
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix WN = this.result.WN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;

        // For transient analysis, RN/WN/AN may not be computed (e.g., FLD passes empty WN).
        // Match MATLAB getAvgTable.m line 114: [RN, WN, AN] = deal(zeros(size(QN)))
        if (Double.isFinite(options.timespan[1]) && QN != null && !QN.isEmpty()) {
            int rows = QN.getNumRows();
            int cols = QN.getNumCols();
            if (RN == null || RN.isEmpty() || RN.getNumRows() != rows || RN.getNumCols() != cols) {
                RN = new Matrix(rows, cols);
            }
            if (WN == null || WN.isEmpty() || WN.getNumRows() != rows || WN.getNumCols() != cols) {
                WN = new Matrix(rows, cols);
            }
            if (AN == null || AN.isEmpty() || AN.getNumRows() != rows || AN.getNumCols() != cols) {
                AN = new Matrix(rows, cols);
            }
        }

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        if (!keepDisabled) {
            Matrix V = new Matrix(sn.nstateful, K);
            for (int i = 0; i < sn.visits.size(); i++) {
                V = V.add(1, sn.visits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.visits.keySet()) {
                    Matrix chainVisits = sn.visits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    if (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0) {
                        int c = -1;
                        for (int row = 0; row < sn.chains.getNumRows(); row++) {
                            if (sn.chains.get(row, k) > 0) {
                                c = row;
                                break;
                            }
                        }
                        Qval.add(QN.get(i, k));
                        Uval.add(UN.get(i, k));
                        Rval.add(RN.get(i, k));
                        ArvR.add(AN.get(i, k));
                        Tval.add(TN.get(i, k));
                        className.add(model.getClasses().get(k).getName());
                        stationName.add(this.model.getStations().get(i).getName());
                        Residval.add(WN.get(i, k));
                    }
                }
            }
            NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setStationNames(stationName);
            return avgTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            Matrix V = new Matrix(sn.nstateful, K);
            for (int i = 0; i < sn.visits.size(); i++) {
                V = V.add(1, sn.visits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.visits.keySet()) {
                    Matrix chainVisits = sn.visits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setStationNames(stationName);
            return avgTable;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgTable(customHandles);
    }

    /**
     * Returns a table of average station metrics organized by job classes with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgTable: " + e.getMessage());
            return null;
        }

        if (this.result == null) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgTable.");
            return null;
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;


        // For transient analysis, RN/WN/AN may not be computed (e.g., FLD passes empty WN).
        // Match MATLAB getAvgTable.m line 114: [RN, WN, AN] = deal(zeros(size(QN)))
        if (Double.isFinite(options.timespan[1]) && QN != null && !QN.isEmpty()) {
            int rows = QN.getNumRows();
            int cols = QN.getNumCols();
            if (RN == null || RN.isEmpty() || RN.getNumRows() != rows || RN.getNumCols() != cols) {
                RN = new Matrix(rows, cols);
            }
            if (WN == null || WN.isEmpty() || WN.getNumRows() != rows || WN.getNumCols() != cols) {
                WN = new Matrix(rows, cols);
            }
            if (AN == null || AN.isEmpty() || AN.getNumRows() != rows || AN.getNumCols() != cols) {
                AN = new Matrix(rows, cols);
            }
        }
        if (QN == null || QN.isEmpty()) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to print AvgTable.");
            return null;
        }

        Matrix V = new Matrix(sn.nstateful, K);
        for (int i = 0; i < sn.visits.size(); i++) {
            V = V.add(1, sn.visits.get(i));
        }
        if (V.isEmpty()) { // SSA
            // Implementation for SSA: sum all chain visits across all chains
            for (Integer chainIndex : sn.visits.keySet()) {
                Matrix chainVisits = sn.visits.get(chainIndex);
                if (chainVisits != null) {
                    V = V.add(1.0, chainVisits);
                }
            }
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0)) {
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                    Residval.add(WN.get(i, k));
                }
            }
        }
        NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgTable.setOptions(this.options);
        avgTable.setClassNames(className);
        avgTable.setStationNames(stationName);
        return avgTable;
    }

    /**
     * Returns a table of average station metrics organized by job classes using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgTable(customHandles, keepDisabled);
    }

    /**
     * Performance sensitivities with respect to service rates.
     *
     * <p>Returns a table with one row per (Station, JobClass) giving the derivative of
     * that row's mean performance measures with respect to that station-class service
     * RATE: dTput_dRate, dRespT_dRate, dQLen_dRate and dUtil_dRate. Every derivative
     * is with respect to the row's OWN rate, i.e. the diagonal of the full parameter
     * Jacobian.</p>
     *
     * <p>Two branches produce the derivatives, selected automatically:</p>
     *
     * <p>"exact": analytic differentiation of a product-form recursion, exact to
     * machine precision and cheaper than a single extra solve. Closed networks use
     * {@code pfqn_sens} (differentiated MVA), open networks the closed-form BCMP
     * sensitivities, their stations being decoupled. Rate derivatives follow the chain
     * rule d(.)/d(rate) = -(L/rate) d(.)/dL, since L(i,r) = visits(i,r)/rate(i,r).
     * Available only on the solvers that evaluate that recursion, {@code SolverMVA}
     * and {@code SolverNC}, and only for single-server queues plus an optional delay,
     * with mixed (open+closed) models excluded.</p>
     *
     * <p>"fd": forward or central finite differences on the CALLING solver's own
     * predictions: the service process at (station,class) is rate-scaled by (1+h), the
     * same solver with the same options is re-run, and the difference quotient is
     * formed. This costs 1+M*R solves (forward) or 2*M*R (central), and it is the only
     * branch that applies to non-product-form models, so it is what every solver other
     * than {@code SolverMVA} and {@code SolverNC} uses.</p>
     *
     * <p>The branch actually taken is reported by
     * {@link NetworkSensitivityTable#getMethod()}. The raw {@code pfqn_sens} result,
     * which carries the full Jacobian that this table shows only the diagonal of, is
     * on {@link NetworkSensitivityTable#getSens()}; it is the second output that the
     * MATLAB reference returns alongside the table, and is null on the open exact
     * branch and on the finite-difference branch.</p>
     *
     * @return table of per-class rate sensitivities
     * @see #getSensitivityTable(String, double, String)
     * @see NetworkSensitivityTable
     * @see #getMomentTable(int[])
     */
    public NetworkSensitivityTable getSensitivityTable() {
        return getSensitivityTable("auto", Double.NaN, "forward");
    }

    /**
     * Performance sensitivities with respect to service rates, with an explicit
     * branch, step and difference scheme. Counterpart of MATLAB's
     * {@code getSensitivityTable(self, 'method', M, 'step', H, 'scheme', S)}.
     *
     * @param method one of "auto" (exact where available and in scope, finite
     *               differences otherwise), "exact", "fd"
     * @param step   relative step of the rate perturbation; NaN selects the default,
     *               1e-4 for the deterministic solvers and 1e-2 for the simulators,
     *               whose Monte Carlo error would otherwise dominate the difference
     *               quotient
     * @param scheme "forward" or "central"
     * @return table of per-class rate sensitivities
     * @throws RuntimeException if the arguments are invalid, or if "exact" is asked of
     *                          a solver or a model that is out of the analytic scope
     * @see #getSensitivityTable()
     */
    public NetworkSensitivityTable getSensitivityTable(String method, double step, String scheme) {
        String meth = (method == null) ? "auto" : method.toLowerCase();
        String sch = (scheme == null) ? "forward" : scheme.toLowerCase();
        if (!"auto".equals(meth) && !"exact".equals(meth) && !"fd".equals(meth)) {
            line_error(mfilename(new Object() {
            }), "The method must be one of 'auto', 'exact', 'fd'.");
        }
        if (!"forward".equals(sch) && !"central".equals(sch)) {
            line_error(mfilename(new Object() {
            }), "The scheme must be 'forward' or 'central'.");
        }

        NetworkStruct sn = model.getStruct();
        int R = sn.nclasses;

        List<Integer> queueIndices = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Queue) {
                queueIndices.add(Integer.valueOf(i));
            }
        }
        int Mq = queueIndices.size();

        boolean exactAvailable = supportsExactSensitivity();
        String scopeMsg = exactSensitivityScope(sn);
        boolean useExact;
        if ("exact".equals(meth)) {
            if (!exactAvailable) {
                line_error(mfilename(new Object() {
                }), "Exact analytic sensitivities differentiate a product-form recursion and are "
                        + "available on SolverMVA and SolverNC only; " + this.getClass().getSimpleName()
                        + " must use 'fd'.");
            }
            if (scopeMsg != null) {
                line_error(mfilename(new Object() {
                }), scopeMsg);
            }
            useExact = true;
        } else if ("fd".equals(meth)) {
            useExact = false;
        } else {
            useExact = exactAvailable && scopeMsg == null;
        }

        SensitivityBlocks blk;
        String methodUsed;
        if (useExact) {
            blk = sensitivityExact(sn, queueIndices, R);
            methodUsed = "exact";
        } else {
            blk = sensitivityFD(sn, queueIndices, R, step, "central".equals(sch));
            methodUsed = "fd";
        }

        List<String> stationName = new ArrayList<String>();
        List<String> className = new ArrayList<String>();
        List<Double> dTput = new ArrayList<Double>();
        List<Double> dRespT = new ArrayList<Double>();
        List<Double> dQLen = new ArrayList<Double>();
        List<Double> dUtil = new ArrayList<Double>();
        for (int ist = 0; ist < Mq; ist++) {
            int node = queueIndices.get(ist).intValue();
            for (int r = 0; r < R; r++) {
                if (!blk.mask[ist][r]) {
                    continue;
                }
                stationName.add(sn.nodenames.get(node));
                className.add(sn.classnames.get(r));
                dTput.add(Double.valueOf(blk.dTput.get(ist, r)));
                dRespT.add(Double.valueOf(blk.dRespT.get(ist, r)));
                dQLen.add(Double.valueOf(blk.dQLen.get(ist, r)));
                dUtil.add(Double.valueOf(blk.dUtil.get(ist, r)));
            }
        }

        NetworkSensitivityTable table = new NetworkSensitivityTable(dTput, dRespT, dQLen, dUtil);
        table.setOptions(this.options);
        table.setStationNames(stationName);
        table.setClassNames(className);
        table.setSens(blk.sens);
        table.setMethod(methodUsed);
        return table;
    }

    /**
     * True when the solver evaluates a product-form recursion that
     * {@link #getSensitivityTable()} can differentiate analytically. False here, so
     * that a solver reaching this base implementation obtains its sensitivities by
     * finite differences on its own predictions. Overridden by {@code SolverMVA} and
     * {@code SolverNC}.
     *
     * @return false
     */
    public boolean supportsExactSensitivity() {
        return false;
    }

    /**
     * The four derivative blocks of a sensitivity computation, plus the row mask and
     * the raw differentiated-MVA result. Java has no multiple returns; this is the
     * MATLAB helpers' output list.
     */
    private static class SensitivityBlocks {
        Matrix dTput;
        Matrix dRespT;
        Matrix dQLen;
        Matrix dUtil;
        boolean[][] mask;
        Ret.pfqnSens sens;

        SensitivityBlocks(int Mq, int R) {
            this.dTput = new Matrix(Mq, R);
            this.dRespT = new Matrix(Mq, R);
            this.dQLen = new Matrix(Mq, R);
            this.dUtil = new Matrix(Mq, R);
            this.mask = new boolean[Mq][R];
            this.sens = null;
        }
    }

    /**
     * Scope of the analytic branch: single-server stations (a delay, with infinite
     * servers, is allowed) and not a mixed open-and-closed model. Class switching is
     * supported: the branch aggregates classes into chains before differentiating.
     *
     * @param sn the network structure
     * @return null when the model is in scope, otherwise MATLAB's rejection message
     */
    private String exactSensitivityScope(NetworkStruct sn) {
        Ret.snGetProductFormParams pf;
        try {
            pf = snGetProductFormParams(sn);
        } catch (Exception e) {
            return "getSensitivityTable could not extract the product-form parameters of this model.";
        }
        Matrix S = pf.S;
        // Multiserver stations are out of scope: pfqn_sens and the closed-form BCMP
        // branch are both single-server statements. An infinite S is a delay, which
        // is allowed. Mirrors MATLAB's any(S(isfinite(S)) > 1).
        if (S != null) {
            for (int i = 0; i < S.length(); i++) {
                double s = S.get(i);
                if (!Double.isNaN(s) && !isInf(s) && s > 1) {
                    return "getSensitivityTable supports single-server stations only.";
                }
            }
        }
        Matrix N = sn.njobs;
        boolean isOpen = false;
        boolean anyFinite = false;
        for (int r = 0; r < sn.nclasses; r++) {
            if (isInf(N.get(r))) {
                isOpen = true;
            } else {
                anyFinite = true;
            }
        }
        if (isOpen && anyFinite) {
            return "getSensitivityTable does not yet support mixed (open+closed) networks.";
        }
        return null;
    }

    /**
     * Analytic branch: differentiated MVA (closed) or closed-form BCMP (open), both
     * evaluated at CHAIN level and then disaggregated back to the classes.
     *
     * <p>A product-form model is solved per chain, not per class: a chain carries the
     * population and the arrival rate, and a class is a share of it. With class
     * switching the two differ, the whole chain population sitting on the reference
     * class, so the recursion must see the chain demands
     * Dc(i,c) = sum_{r in c} D(i,r), Nc(c) = sum_{r in c} N(r), Zc likewise. The
     * parameter of the table is still the per-class rate mu(i,r), which enters exactly
     * one chain demand, giving the chain rule
     * d(.)/dmu(i,r) = -(D(i,r)/mu(i,r)) d(.)/dDc(i,c). The per-class measures are then
     * composed from the chain ones, which is also what fixes the visit ratios: a class
     * throughput at a station is X_c*v(i,r) and a per-visit response time is
     * Q(i,r)/T(i,r), whereas the chain quantities are per chain and per visit-chain
     * respectively.</p>
     *
     * @param sn           the network structure
     * @param queueIndices the node indices of the queues, in table order
     * @param R            the number of classes
     * @return the derivative blocks
     */
    private SensitivityBlocks sensitivityExact(NetworkStruct sn, List<Integer> queueIndices, int R) {
        Ret.snGetProductFormParams pf = snGetProductFormParams(sn);
        Matrix lambda = pf.lambda;
        Matrix D = pf.D;
        Matrix Np = pf.N;
        Matrix Z = pf.Z;
        Matrix N = sn.njobs;

        boolean isOpen = false;
        for (int r = 0; r < R; r++) {
            if (isInf(N.get(r))) {
                isOpen = true;
            }
        }
        int Mq = queueIndices.size();
        int C = sn.nchains;
        SensitivityBlocks blk = new SensitivityBlocks(Mq, R);
        double[][] rates = new double[Mq][R];
        for (int ist = 0; ist < Mq; ist++) {
            int sIdx = (int) sn.nodeToStation.get(queueIndices.get(ist).intValue());
            for (int r = 0; r < R; r++) {
                rates[ist][r] = sn.rates.get(sIdx, r);
                blk.mask[ist][r] = isFiniteRate(rates[ist][r]) && rates[ist][r] > 0
                        && D.get(ist, r) > 0;
            }
        }

        // chainOf[r] is the chain class r belongs to; Dc, Zc, Nc and lambdac
        // aggregate the per-class quantities over the classes of each chain.
        Matrix Zrow = sumColumns(Z, R);
        int[] chainOf = new int[R];
        Matrix Dc = new Matrix(Mq, C);
        double[] Zc = new double[C];
        double[] Nc = new double[C];
        double[] lambdac = new double[C];
        for (int c = 0; c < C; c++) {
            for (int r = 0; r < R; r++) {
                if (sn.chains.get(c, r) == 0) {
                    continue;
                }
                chainOf[r] = c;
                for (int ist = 0; ist < Mq; ist++) {
                    Dc.set(ist, c, Dc.get(ist, c) + D.get(ist, r));
                }
                Zc[c] += Zrow.get(0, r);
                if (!isInf(Np.get(r)) && !Double.isNaN(Np.get(r))) {
                    Nc[c] += Np.get(r);
                }
                lambdac[c] += lambda.get(r);
            }
        }

        if (!isOpen) {
            // ---- closed branch: differentiated MVA at chain level ------------
            Matrix Ncm = new Matrix(1, C);
            Matrix Zcm = new Matrix(1, C);
            for (int c = 0; c < C; c++) {
                Ncm.set(0, c, Nc[c]);
                Zcm.set(0, c, Zc[c]);
            }
            blk.sens = pfqn_sens(Dc, Ncm, Zcm);
            for (int ist = 0; ist < Mq; ist++) {
                for (int r = 0; r < R; r++) {
                    if (!blk.mask[ist][r]) {
                        continue;
                    }
                    int c = chainOf[r];
                    if (Dc.get(ist, c) <= 0) {
                        continue;
                    }
                    double rate = rates[ist][r];
                    double Dir = D.get(ist, r);
                    double visits = Dir * rate;       // chain-normalized visit ratio
                    int p = ist * C + c;
                    double chain = -Dir / rate;       // dDc(i,c)/dmu(i,r)
                    double Xc = blk.sens.X.get(c);
                    double Qc = blk.sens.Q.get(ist, c);
                    double dXc = blk.sens.dX.get(c, p) * chain;
                    double dQc = blk.sens.dQ[p].get(ist, c) * chain;
                    // Class share of the chain queue at this station, and its own
                    // dependence on the rate.
                    double alpha = Dir / Dc.get(ist, c);
                    double dalpha = chain * (Dc.get(ist, c) - Dir)
                            / (Dc.get(ist, c) * Dc.get(ist, c));
                    double Qir = alpha * Qc;
                    double dQir = dalpha * Qc + alpha * dQc;
                    double Tir = Xc * visits;
                    double dTir = dXc * visits;
                    blk.dTput.set(ist, r, dTir);
                    blk.dQLen.set(ist, r, dQir);
                    blk.dUtil.set(ist, r, dXc * Dir + Xc * chain);
                    if (Tir > 0) {
                        // Per-visit response time by Little's law, R = Q/T.
                        blk.dRespT.set(ist, r, (dQir * Tir - Qir * dTir) / (Tir * Tir));
                    }
                }
            }
        } else {
            // ---- open branch: exact closed-form BCMP -------------------------
            // rho(i,r) = lambdac(c)*D(i,r) with c the chain of r; U(i) = sum_r
            // rho(i,r). The stations decouple, so only the own service rate mu(i,r)
            // moves the measures at (i,r). The throughput T(i,r) = lambdac(c)*v(i,r)
            // is fixed by the arrival rate, hence its rate derivative is exactly zero.
            Matrix rho = new Matrix(Mq, R);
            for (int ist = 0; ist < Mq; ist++) {
                for (int r = 0; r < R; r++) {
                    if (D.get(ist, r) > 0) {
                        rho.set(ist, r, lambdac[chainOf[r]] * D.get(ist, r));
                    }
                }
            }
            for (int ist = 0; ist < Mq; ist++) {
                double Ui = 0.0;
                for (int r = 0; r < R; r++) {
                    Ui += rho.get(ist, r);
                }
                double denom = 1 - Ui;
                for (int r = 0; r < R; r++) {
                    if (!blk.mask[ist][r]) {
                        continue;
                    }
                    double rate = rates[ist][r];
                    double svct = 1.0 / rate;         // per-visit service time
                    double drho = -rho.get(ist, r) / rate;
                    double dU = drho;                 // own class only
                    double dsvct = -svct / rate;
                    blk.dRespT.set(ist, r, (dsvct * denom + svct * dU) / (denom * denom));
                    blk.dQLen.set(ist, r, (drho * denom + rho.get(ist, r) * dU) / (denom * denom));
                    blk.dUtil.set(ist, r, drho);
                    blk.dTput.set(ist, r, 0.0);   // open throughput = lambda*visits (fixed)
                }
            }
        }
        return blk;
    }

    /**
     * Finite-difference branch: re-run the calling solver on rate-perturbed copies of
     * the model. The perturbation is a pure time scaling of the service process, so
     * the shape of the distribution, and in particular its SCV, is preserved and only
     * the rate moves.
     *
     * @param sn           the network structure
     * @param queueIndices the node indices of the queues, in table order
     * @param R            the number of classes
     * @param step         the relative rate step, or NaN for the default
     * @param central      true for central differences, false for forward
     * @return the derivative blocks
     */
    private SensitivityBlocks sensitivityFD(NetworkStruct sn, List<Integer> queueIndices, int R,
                                            double step, boolean central) {
        int Mq = queueIndices.size();
        boolean isSimulation = (this instanceof SolverSSA) || (this instanceof SolverJMT)
                || (this instanceof SolverLDES);

        double h = step;
        if (Double.isNaN(h)) {
            h = isSimulation ? 1e-2 : 1e-4;
        }
        if (Double.isInfinite(h) || h <= 0 || h >= 1) {
            line_error(mfilename(new Object() {
            }), "The finite-difference step must be a scalar in (0,1).");
        }
        if (isSimulation && this.options != null && this.options.seed <= 0) {
            // Common random numbers: the base and perturbed runs must share a seed,
            // otherwise the difference quotient measures Monte Carlo noise. Every run
            // reuses this solver's options object, so pinning it here pins it for all.
            this.options.seed = 23000;
        }

        // The rates are read before any perturbation: refreshProcesses rewrites
        // sn.rates in place, and the difference quotient needs the unperturbed value.
        Matrix rates0 = sn.rates.copy();
        boolean[][] visited = visitMask(sn);
        SensitivityBlocks blk = new SensitivityBlocks(Mq, R);
        for (int ist = 0; ist < Mq; ist++) {
            int sIdx = (int) sn.nodeToStation.get(queueIndices.get(ist).intValue());
            for (int r = 0; r < R; r++) {
                double rate = rates0.get(sIdx, r);
                blk.mask[ist][r] = isFiniteRate(rate) && rate > 0 && visited[sIdx][r];
            }
        }

        SolverResult base = solveOnceForSensitivity();

        List<Station> stations = model.getStations();
        List<JobClass> classes = model.getClasses();
        for (int ist = 0; ist < Mq; ist++) {
            int sIdx = (int) sn.nodeToStation.get(queueIndices.get(ist).intValue());
            Station station = stations.get(sIdx);
            for (int r = 0; r < R; r++) {
                if (!blk.mask[ist][r]) {
                    continue;
                }
                JobClass jobclass = classes.get(r);
                double rate = rates0.get(sIdx, r);
                Distribution baseDist = station.getServiceProcess(jobclass);
                SolverResult plus;
                SolverResult minus;
                double denom;
                try {
                    setServiceForSensitivity(station, jobclass,
                            DistributionScaling.scaleRate(baseDist, 1 + h));
                    plus = solveOnceForSensitivity();
                    if (central) {
                        setServiceForSensitivity(station, jobclass,
                                DistributionScaling.scaleRate(baseDist, 1 - h));
                        minus = solveOnceForSensitivity();
                        denom = 2 * rate * h;
                    } else {
                        minus = base;
                        denom = rate * h;
                    }
                } finally {
                    setServiceForSensitivity(station, jobclass, baseDist);
                }
                blk.dTput.set(ist, r, (plus.TN.get(sIdx, r) - minus.TN.get(sIdx, r)) / denom);
                blk.dRespT.set(ist, r, (plus.RN.get(sIdx, r) - minus.RN.get(sIdx, r)) / denom);
                blk.dQLen.set(ist, r, (plus.QN.get(sIdx, r) - minus.QN.get(sIdx, r)) / denom);
                blk.dUtil.set(ist, r, (plus.UN.get(sIdx, r) - minus.UN.get(sIdx, r)) / denom);
            }
        }
        return blk;
    }

    /**
     * Solves with a fresh instance of the calling solver class, carrying its options
     * over so that seeds, tolerances and method selection are those of the caller. A
     * fresh instance is used because a solver caches its results, and a reset alone
     * would not discard a warm start.
     *
     * @return the mean performance measures under this solver's method
     */
    private SolverResult solveOnceForSensitivity() {
        try {
            // the equivalent of MATLAB's feval(class(self), model, self.options):
            // every Network solver exposes a (Network, SolverOptions) constructor
            Constructor<? extends NetworkSolver> ctor =
                    this.getClass().getConstructor(Network.class, SolverOptions.class);
            NetworkSolver solver = ctor.newInstance(this.model, getOptions());
            return solver.getAvg();
        } catch (Exception e) {
            throw new RuntimeException("the finite-difference sensitivity path could not re-run "
                    + this.getClass().getSimpleName() + ": it must expose a "
                    + "(Network, SolverOptions) constructor", e);
        }
    }

    /**
     * Installs a service process at a station-class and refreshes the model, so that
     * the cached structure carries the perturbed rate.
     *
     * @param station  the station
     * @param jobclass the job class
     * @param distrib  the service process to install
     */
    private void setServiceForSensitivity(Station station, JobClass jobclass, Distribution distrib) {
        station.setService(jobclass, distrib);
        model.refreshProcesses();
    }

    /**
     * A (nstations x nclasses) mask of the station-class pairs that carry visits, used
     * in place of the demand matrix, which a non-product-form model may not admit.
     *
     * @param sn the network structure
     * @return the visit mask
     */
    private static boolean[][] visitMask(NetworkStruct sn) {
        boolean[][] visited = new boolean[sn.nstations][sn.nclasses];
        for (int c = 0; c < sn.nchains; c++) {
            Matrix V = sn.visits.get(Integer.valueOf(c));
            if (V == null || V.isEmpty()) {
                continue;
            }
            for (int i = 0; i < sn.nstations; i++) {
                for (int r = 0; r < sn.nclasses; r++) {
                    if (V.get(i, r) > GlobalConstants.Zero) {
                        visited[i][r] = true;
                    }
                }
            }
        }
        return visited;
    }

    /**
     * MATLAB's isfinite for a rate: false for both NaN and an infinite rate.
     *
     * @param rate the rate
     * @return true if the rate is a finite number
     */
    private static boolean isFiniteRate(double rate) {
        return !Double.isNaN(rate) && !isInf(rate);
    }

    /**
     * Exact higher moments of the per-class performance measures, up to the second
     * moment.
     *
     * @return table of per-class queue-length and response-time moments
     * @see #getMomentTable(int[])
     */
    public NetworkMomentTable getMomentTable() {
        return getMomentTable(2);
    }

    /**
     * Exact higher moments of the per-class performance measures, up to the given
     * moment order.
     *
     * <p>A scalar order k is read as the set 1..k, "everything up to order k", so
     * that {@code getMomentTable(2)} means "up to the second moment" and not "the
     * second moment alone". Use {@link #getMomentTable(int[])} to select an explicit
     * set of orders.</p>
     *
     * @param order the highest moment order to report, an integer in 1..3
     * @return table of per-class queue-length and response-time moments
     * @see #getMomentTable(int[])
     */
    public NetworkMomentTable getMomentTable(int order) {
        return getMomentTable(momentOrderUpTo(order, 3));
    }

    /**
     * Exact higher moments of the per-class performance measures.
     *
     * <p>Returns a table with one row per (Station, JobClass) giving, in addition to
     * the means that {@link #getAvgTable()} reports, the higher moments of that
     * row's queue length and response time. {@code order} is a SET of moment orders,
     * taken literally: order 1 contributes QLen and RespT, order 2 contributes
     * QLenVar, QLenSCV, RespTVar and RespTSCV, order 3 contributes QLenSkew and
     * RespTSkew. So {@code new int[]{1,2}} is the default and is the same as the
     * scalar 2, while {@code new int[]{2,3}} selects the second moments and the
     * skewness without the means.</p>
     *
     * <p>Order 3 also adds QLenSkew, the skewness of the per-class queue length.
     * That quantity is reachable because the generating parameter need not scale a
     * whole demand column: scaling L(i,r) alone is Theorem 1 of Akyildiz and Strelen
     * with the class subset T = {r}, and it generates the moments of n(i,r) itself.
     * QLenSkew is available only for closed single-server models, which is the scope
     * of {@code Pfqn_sens_mom}; it is NaN otherwise.</p>
     *
     * <p>All of it is exact, not simulated and not approximated. The queue-length
     * moments come from the product-form identity Cov[n(i,r),n(j,s)] = L(j,s)
     * dQ(i,r)/dL(j,s), evaluated by the {@code Pfqn_sens_*} family.</p>
     *
     * <p>RESPONSE-TIME MOMENTS ARE FCFS-ONLY. RespTVar, RespTSCV and RespTSkew are
     * NaN at any station that is not FCFS, and in open or mixed models. This is a
     * limitation of the theory, not of the implementation: the sojourn-time
     * distribution at a processor-sharing or LCFS center is not known in general
     * (Strelen 1990, Section 4), so there is no correct value to report and a wrong
     * one is worse than a blank. The mean RespT is always reported, since it needs no
     * distributional result.</p>
     *
     * <p>Scope by model type: closed single-server uses {@code pfqn_sens_mva};
     * closed multiserver and mixed open/closed use {@code pfqn_sens_mvaldmx}; purely
     * open single-server uses the exact BCMP closed form; purely open multiserver is
     * not supported.</p>
     *
     * <p>The raw results, including the full covariance matrices that this table
     * shows only the diagonal of, are on {@link NetworkMomentTable#getMoments()}.</p>
     *
     * @param order the set of moment orders to report; each must be an integer in
     *              1..3. Null or empty selects the default, up to the second moment.
     * @return table of per-class queue-length and response-time moments
     * @see #getMomentStationTable(int[])
     * @see #getAvgTable()
     */
    public NetworkMomentTable getMomentTable(int[] order) {
        if (order == null || order.length == 0) {
            order = momentOrderUpTo(2, 3);
        }
        order = validateMomentOrder(order, 3);
        int maxOrder = order[order.length - 1];

        NetworkStruct snl = model.getStruct();
        int R = snl.nclasses;
        Matrix N = snl.njobs;

        Ret.snGetProductFormParams pf = snGetProductFormParams(snl);
        Matrix lambda = pf.lambda;
        Matrix D = pf.D;
        Matrix Np = pf.N;
        Matrix Z = pf.Z;
        Matrix mu = pf.mu;
        Matrix Ssrv = pf.S;

        List<Integer> queueIndices = new ArrayList<Integer>();
        for (int i = 0; i < snl.nodetype.size(); i++) {
            if (snl.nodetype.get(i) == NodeType.Queue) {
                queueIndices.add(Integer.valueOf(i));
            }
        }
        int Mq = queueIndices.size();
        Matrix Ztot = sumColumns(Z, R);

        boolean isOpen = false;
        boolean isClosed = false;
        for (int r = 0; r < R; r++) {
            double nr = N.get(r);
            if (isInf(nr)) {
                isOpen = true;
            } else if (nr > 0) {
                isClosed = true;
            }
        }
        boolean isMixed = isOpen && isClosed;

        NetworkMomentResult mom = new NetworkMomentResult();
        Matrix QLen = new Matrix(Mq, R);
        Matrix QLenVar = new Matrix(Mq, R);

        // ---- queue-length moments ------------------------------------------
        if (!isOpen) {
            boolean allSingleServer = true;
            for (int i = 0; i < Ssrv.length(); i++) {
                if (Ssrv.get(i) != 1.0) {
                    allSingleServer = false;
                    break;
                }
            }
            if (allSingleServer) {
                mom.qlenMva = pfqn_sens_mva(D, Np, Ztot);
                mom.X = mom.qlenMva.X;
                mom.QCov = mom.qlenMva.QCov;
                QLen = mom.qlenMva.Q;
                QLenVar = mom.qlenMva.QVar;
            } else {
                mom.qlenLdmx = pfqn_sens_mvaldmx(new Matrix(1, R), D, Np, Ztot, mu, Ssrv);
                mom.X = mom.qlenLdmx.X;
                mom.QCov = mom.qlenLdmx.QCov;
                QLen = mom.qlenLdmx.Q;
                QLenVar = mom.qlenLdmx.QVar;
            }
        } else if (isMixed) {
            mom.qlenLdmx = pfqn_sens_mvaldmx(lambda, D, Np, Ztot, mu, Ssrv);
            mom.X = mom.qlenLdmx.X;
            mom.QCov = mom.qlenLdmx.QCov;
            QLen = mom.qlenLdmx.Q;
            QLenVar = mom.qlenLdmx.QVar;
        } else {
            // Purely open. The moment analysis of the Pfqn_sens_* family needs at
            // least one closed class to have a population lattice to recurse on, but
            // a purely open BCMP single-server station has a closed-form joint law
            // that needs no recursion at all: with rho(i,r) = lambda(r)*D(i,r) and
            // rho_i = sum_r rho(i,r) < 1,
            //   P(n_i) = (1-rho_i)*(sum_r n(i,r))!/prod_r n(i,r)! * prod_r rho(i,r)^n(i,r)
            // so the total n_i is geometric with parameter rho_i and, conditionally
            // on n_i, the classes are multinomial with p_r = rho(i,r)/rho_i.
            // Compounding,
            //   Cov[n(i,r),n(i,s)] = E[n_i]*(delta_rs*p_r - p_r*p_s) + p_r*p_s*Var[n_i]
            for (int i = 0; i < Ssrv.length(); i++) {
                if (Ssrv.get(i) > 1) {
                    line_error(mfilename(new Object() {
                    }), "getMomentTable does not support multiserver stations in a purely open "
                            + "model: the queue-length law is not geometric there. Add a closed "
                            + "class, or use a single-server model.");
                }
            }
            Matrix rho = new Matrix(Mq, R);
            for (int ist = 0; ist < Mq; ist++) {
                for (int r = 0; r < R; r++) {
                    if (isInf(N.get(r))) {
                        rho.set(ist, r, lambda.get(r) * D.get(ist, r));
                    }
                }
            }
            for (int ist = 0; ist < Mq; ist++) {
                double ri = 0.0;
                for (int r = 0; r < R; r++) {
                    ri += rho.get(ist, r);
                }
                if (ri >= 1) {
                    line_error(mfilename(new Object() {
                    }), String.format("Station %s is unstable (utilization %.4f >= 1); its "
                                    + "queue-length moments do not exist.",
                            snl.nodenames.get(queueIndices.get(ist).intValue()), ri));
                }
                if (ri <= 0) {
                    continue;
                }
                double En = ri / (1 - ri);
                double Vn = ri / ((1 - ri) * (1 - ri));
                for (int r = 0; r < R; r++) {
                    double pr = rho.get(ist, r) / ri;
                    QLen.set(ist, r, pr * En);
                    QLenVar.set(ist, r, En * (pr - pr * pr) + pr * pr * Vn);
                }
            }
        }
        mom.Q = QLen;
        mom.QVar = QLenVar;

        // ---- per-class queue-length skewness (closed, single-server only) ----
        // Pfqn_sens_mom with groups = 1..R scales one class at a time, which is the
        // class-subset parameter T = {r} of Akyildiz and Strelen's Theorem 1, so it
        // yields the moments of n(i,r) rather than of the station total.
        Matrix QLenSkew = new Matrix(Mq, R);
        for (int ist = 0; ist < Mq; ist++) {
            for (int r = 0; r < R; r++) {
                QLenSkew.set(ist, r, Double.NaN);
            }
        }
        boolean allSingleServerSkew = true;
        for (int i = 0; i < Ssrv.length(); i++) {
            if (Ssrv.get(i) != 1.0) {
                allSingleServerSkew = false;
                break;
            }
        }
        if (hasOrder(order, 3) && !isOpen && allSingleServerSkew) {
            mom.qlenmom = pfqn_sens_mom(D, Np, Ztot, Matrix.ones(1, Mq),
                    Pfqn_sens_mom.perClassGroups(R));
            QLenSkew = mom.qlenmom.Skew;
        }

        // ---- response-time moments (FCFS only) ------------------------------
        Matrix RespT = new Matrix(Mq, R);
        Matrix RespTVar = new Matrix(Mq, R);
        Matrix RespTSkew = new Matrix(Mq, R);
        for (int ist = 0; ist < Mq; ist++) {
            for (int r = 0; r < R; r++) {
                RespTVar.set(ist, r, Double.NaN);
                RespTSkew.set(ist, r, Double.NaN);
            }
        }
        // Pfqn_sens_respt carries one service time S(i) per station, but S enters its
        // queue-length recursion ONLY through the product rho(i,r) = S(i)*V(i,r) =
        // the demand. The rate mu(i) = 1/S(i) is read directly by equation (4.5)
        // alone, and (4.5) is only evaluated at FCFS stations. So a non-FCFS station
        // may keep class-dependent demands: give it any nominal S and let V absorb
        // the rest. Only the FCFS stations must have class-independent rates, which
        // BCMP requires of them anyway.
        boolean respAvail = !isOpen && fcfsRatesAreClassIndependent(snl, queueIndices, R);
        if (respAvail) {
            Matrix Ssvc = Matrix.ones(Mq, 1);
            Matrix Vq = new Matrix(Mq, R);
            for (int ist = 0; ist < Mq; ist++) {
                int sIdx = (int) snl.nodeToStation.get(queueIndices.get(ist).intValue());
                if (snl.sched.get(snl.stations.get(sIdx)) == SchedStrategy.FCFS) {
                    double st = serviceTimeOf(snl, sIdx, R);
                    if (st > 0) {
                        Ssvc.set(ist, 0, st);   // the true per-visit rate; (4.5) needs it
                    }
                }
                for (int r = 0; r < R; r++) {
                    Vq.set(ist, r, D.get(ist, r) / Ssvc.get(ist, 0));
                }
            }
            mom.respt = pfqn_sens_respt(Ssvc, Vq, Np, Ztot, Ssrv, maxOrder);
        }
        for (int ist = 0; ist < Mq; ist++) {
            int sIdx = (int) snl.nodeToStation.get(queueIndices.get(ist).intValue());
            for (int r = 0; r < R; r++) {
                if (D.get(ist, r) <= 0) {
                    continue;
                }
                // Mean response time per visit, by Little's law at the station: the
                // arrival rate of class r to station i is X(r) times the visit ratio
                // V(i,r) = D(i,r)/S(i,r) = D(i,r)*rate(i,r). This needs no
                // distributional result, so it is always reported.
                double xr = throughputOf(mom, lambda, N, r);
                double Vir = D.get(ist, r) * snl.rates.get(sIdx, r);
                if (xr > 0 && Vir > 0) {
                    RespT.set(ist, r, QLen.get(ist, r) / (xr * Vir));
                }
            }
            // The variance needs the sojourn-time distribution, which is known only
            // at FCFS centers; elsewhere RespTVar stays NaN.
            if (mom.respt != null && snl.sched.get(snl.stations.get(sIdx)) == SchedStrategy.FCFS) {
                for (int r = 0; r < R; r++) {
                    if (D.get(ist, r) > 0) {
                        RespT.set(ist, r, mom.respt.W.get(ist, r));
                        if (maxOrder >= 2) {
                            RespTVar.set(ist, r, mom.respt.WVar.get(ist, r));
                        }
                        if (maxOrder >= 3) {
                            RespTSkew.set(ist, r, mom.respt.WSkew.get(ist, r));
                        }
                    }
                }
            }
        }

        // ---- assemble -------------------------------------------------------
        List<String> stationName = new ArrayList<String>();
        List<String> className = new ArrayList<String>();
        List<Double> QLenv = new ArrayList<Double>();
        List<Double> QLenVarv = new ArrayList<Double>();
        List<Double> QLenSCVv = new ArrayList<Double>();
        List<Double> RespTv = new ArrayList<Double>();
        List<Double> RespTVarv = new ArrayList<Double>();
        List<Double> RespTSCVv = new ArrayList<Double>();
        List<Double> RespTSkewv = new ArrayList<Double>();
        List<Double> QLenSkewv = new ArrayList<Double>();
        for (int ist = 0; ist < Mq; ist++) {
            int node = queueIndices.get(ist).intValue();
            for (int r = 0; r < R; r++) {
                if (D.get(ist, r) <= 0) {
                    continue;   // class r does not visit this station
                }
                stationName.add(snl.nodenames.get(node));
                className.add(snl.classnames.get(r));
                QLenv.add(Double.valueOf(QLen.get(ist, r)));
                QLenVarv.add(Double.valueOf(QLenVar.get(ist, r)));
                QLenSCVv.add(Double.valueOf(momentScv(QLenVar.get(ist, r), QLen.get(ist, r))));
                RespTv.add(Double.valueOf(RespT.get(ist, r)));
                RespTVarv.add(Double.valueOf(RespTVar.get(ist, r)));
                RespTSCVv.add(Double.valueOf(momentScv(RespTVar.get(ist, r), RespT.get(ist, r))));
                RespTSkewv.add(Double.valueOf(RespTSkew.get(ist, r)));
                QLenSkewv.add(Double.valueOf(QLenSkew.get(ist, r)));
            }
        }

        List<String> names = new ArrayList<String>();
        List<List<Double>> cols = new ArrayList<List<Double>>();
        if (hasOrder(order, 1)) {
            names.add("QLen");
            cols.add(QLenv);
        }
        if (hasOrder(order, 2)) {
            names.add("QLenVar");
            cols.add(QLenVarv);
            names.add("QLenSCV");
            cols.add(QLenSCVv);
        }
        if (hasOrder(order, 3)) {
            names.add("QLenSkew");
            cols.add(QLenSkewv);
        }
        if (hasOrder(order, 1)) {
            names.add("RespT");
            cols.add(RespTv);
        }
        if (hasOrder(order, 2)) {
            names.add("RespTVar");
            cols.add(RespTVarv);
            names.add("RespTSCV");
            cols.add(RespTSCVv);
        }
        if (hasOrder(order, 3)) {
            names.add("RespTSkew");
            cols.add(RespTSkewv);
        }

        NetworkMomentTable table = new NetworkMomentTable(names, cols);
        table.setOptions(this.options);
        table.setStationNames(stationName);
        table.setClassNames(className);
        table.setMoments(mom);
        return table;
    }

    /**
     * Higher moments of the total queue length, up to the second moment.
     *
     * @return table of per-station total queue-length moments
     * @see #getMomentStationTable(int[])
     */
    public NetworkMomentStationTable getMomentStationTable() {
        return getMomentStationTable(momentOrderUpTo(2, 3));
    }

    /**
     * Higher moments of the total queue length, up to the given moment order.
     *
     * @param order the highest moment order to report, an integer in 1..3
     * @return table of per-station total queue-length moments
     * @see #getMomentStationTable(int[])
     */
    public NetworkMomentStationTable getMomentStationTable(int order) {
        return getMomentStationTable(momentOrderUpTo(order, 3));
    }

    /**
     * Higher moments of the total queue length.
     *
     * <p>Returns a table with one row per Station giving the moments of the TOTAL
     * queue length at that station, Q_i = sum_r n(i,r). {@code order} is a SET of
     * moment orders, taken literally: order 1 contributes QLen, order 2 contributes
     * QLenVar and QLenSCV, order 3 contributes QLenM3 and QLenSkew. So
     * {@code new int[]{1,3}} gives the mean and the third moment without the
     * variance.</p>
     *
     * <p>Why this is a separate table from {@link #getMomentTable(int[])}. The
     * generating parameter here is x_i, the reciprocal of the capacity of station i,
     * which scales the service times of ALL classes at that station at once, so the
     * moments it produces are those of the station total. The per-class and
     * per-chain groupings of the same recursion are in
     * {@link #getMomentTable(int[])} and {@link #getMomentChainTable(int[])}
     * respectively; all three are Theorem 1 of Akyildiz and Strelen under different
     * class subsets T, of which Strelen's whole-column x_i is the case T = all
     * classes. The three are consistent: Var[Q_i] here equals the sum of the
     * per-class covariances of {@link #getMomentTable(int[])} at station i over all
     * class pairs.</p>
     *
     * <p>The ALGORITHM is chosen by the solver's method, set at construction, not by
     * an argument here: a method is a property of the solver object, so passing one
     * per call would let a single solver answer with two different algorithms.</p>
     *
     * <ul>
     *   <li>{@code new SolverMVA(model)} -&gt; {@code pfqn_sens_mom}, exact, but it
     *       walks the whole population lattice at a cost of prod(N+1), so it is
     *       unusable once the populations are large.</li>
     *   <li>{@code new SolverMVA(model,"method","lin")} -&gt;
     *       {@code pfqn_sens_linearizer}, approximate and polynomial-time. The
     *       reference reports relative errors below 2.1% on E[Q], 4.1% on E[Q^2] and
     *       6.2% on E[Q^3].</li>
     * </ul>
     *
     * <p>Any Linearizer-family method ({@code "lin"}, {@code "amva.lin"},
     * {@code "egflin"}, {@code "gflin"}) takes the approximate path; every other
     * method takes the exact one.</p>
     *
     * <p>Restricted to closed single-server models. Mixed and open second moments
     * are available per class from {@link #getMomentTable(int[])}; the higher moments
     * of the reference are stated for closed load-independent networks only.</p>
     *
     * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
     * its Linearizer", Performance Evaluation 11:127-142, 1990, Theorem 3.1 and
     * equation (3.2).</p>
     *
     * @param order the set of moment orders to report; each must be an integer in
     *              1..3. Null or empty selects the default, up to the second moment.
     * @return table of per-station total queue-length moments
     * @see #getMomentTable(int[])
     */
    public NetworkMomentStationTable getMomentStationTable(int[] order) {
        if (order == null || order.length == 0) {
            order = momentOrderUpTo(2, 3);
        }
        order = validateMomentOrder(order, 3);
        // The algorithm is a property of the solver, not of this call: it comes from
        // the method set at construction, e.g. new SolverMVA(model,"method","lin").
        // Taking a per-call method argument here would have given the same solver
        // object two disagreeing methods.
        String method = getOptions().method;

        NetworkStruct snl = model.getStruct();
        int R = snl.nclasses;
        Matrix N = snl.njobs;

        Ret.snGetProductFormParams pf = snGetProductFormParams(snl);
        Matrix D = pf.D;
        Matrix Np = pf.N;
        Matrix Z = pf.Z;
        Matrix Ssrv = pf.S;

        List<Integer> queueIndices = new ArrayList<Integer>();
        for (int i = 0; i < snl.nodetype.size(); i++) {
            if (snl.nodetype.get(i) == NodeType.Queue) {
                queueIndices.add(Integer.valueOf(i));
            }
        }
        int Mq = queueIndices.size();
        Matrix Ztot = sumColumns(Z, R);

        // The moments rest on the product-form identity Cov = L dQ/dL, which is a
        // theorem about the product form: it follows from a d/da log G = E[n], which
        // needs the state distribution to be exponential-family in the demands.
        // Outside product form, L dQ/dL remains computable and is simply NOT a
        // covariance, so differentiating a non-product-form solver would return a
        // confident wrong number. Refuse rather than do that; exact moments of such a
        // model come from SolverCTMC's stationary distribution, or from LDES via
        // setReward.
        if (!snHasProductForm(snl)) {
            line_error(mfilename(new Object() {
            }), "getMomentStationTable requires a product-form model: the moment identity "
                    + "Cov[n,n] = L dQ/dL holds only under product form, so no correct value "
                    + "exists here. Use SolverCTMC (exact distribution) or SolverLDES with "
                    + "setReward for the moments of a non-product-form model.");
        }
        for (int r = 0; r < R; r++) {
            if (isInf(N.get(r))) {
                line_error(mfilename(new Object() {
                }), "getMomentStationTable supports closed models only. Per-class second moments "
                        + "of an open or mixed model are available from getMomentTable.");
            }
        }
        for (int i = 0; i < Ssrv.length(); i++) {
            if (Ssrv.get(i) > 1) {
                line_error(mfilename(new Object() {
                }), "getMomentStationTable supports single-server stations only: the "
                        + "higher-moment recursion of the reference is stated for "
                        + "load-independent stations. Per-class second moments of a multiserver "
                        + "model are available from getMomentTable.");
            }
        }

        // Three paths exist, and the solver's method selects between them: the two
        // hand-differentiated algorithms where they apply, and the general
        // numerical-derivative oracle everywhere else. The solver has already rejected
        // methods it does not support via listValidMethods, so no further validation
        // is needed here.
        NetworkMomentStationResult mom = new NetworkMomentStationResult();
        if (isLinearizerMethod(method)) {
            mom.linearizer = pfqn_sens_linearizer(D, Np, Ztot);
            mom.m = mom.linearizer.m;
            mom.Var = mom.linearizer.Var;
            mom.Cov = mom.linearizer.Cov;
            mom.M2 = mom.linearizer.M2;
            mom.M3 = mom.linearizer.M3;
            mom.Skew = mom.linearizer.Skew;
            mom.CovAsym = mom.linearizer.CovAsym;
        } else if (isExactMvaMethod(method)) {
            mom.exact = pfqn_sens_mom(D, Np, Ztot);
            mom.m = mom.exact.m;
            mom.Var = mom.exact.Var;
            mom.Cov = mom.exact.Cov;
            mom.CovG = mom.exact.CovG;
            mom.M2 = mom.exact.M2;
            mom.M3 = mom.exact.M3;
            mom.Skew = mom.exact.Skew;
            mom.CovAsym = mom.exact.CovAsym;
        } else {
            // A method with no hand-differentiated counterpart. The identity does not
            // care HOW the mean queue lengths were obtained, so THIS solver is used as
            // a mean-value oracle and differentiated numerically. This is exactly the
            // move the Linearizer makes analytically (Strelen Sec. 5: differentiate the
            // approximate fixed point, then apply the exact identity), generalized to
            // any product-form solver and method. The moments inherit the accuracy of
            // that method's means, so an exact oracle gives exact moments.
            mom = momentsByFiniteDifference(snl, queueIndices);
        }

        List<String> stationName = new ArrayList<String>();
        List<Double> QLen = new ArrayList<Double>();
        List<Double> QLenVar = new ArrayList<Double>();
        List<Double> QLenSCV = new ArrayList<Double>();
        List<Double> QLenM3 = new ArrayList<Double>();
        List<Double> QLenSkew = new ArrayList<Double>();
        for (int ist = 0; ist < Mq; ist++) {
            boolean visited = false;
            for (int r = 0; r < R; r++) {
                if (D.get(ist, r) > 0) {
                    visited = true;
                    break;
                }
            }
            if (!visited) {
                continue;   // no class visits this station
            }
            stationName.add(snl.nodenames.get(queueIndices.get(ist).intValue()));
            QLen.add(Double.valueOf(mom.m.get(ist)));
            QLenVar.add(Double.valueOf(mom.Var.get(ist)));
            if (mom.m.get(ist) > 0) {
                QLenSCV.add(Double.valueOf(mom.Var.get(ist) / (mom.m.get(ist) * mom.m.get(ist))));
            } else {
                QLenSCV.add(Double.valueOf(Double.NaN));
            }
            QLenM3.add(Double.valueOf(mom.M3.get(ist)));
            QLenSkew.add(Double.valueOf(mom.Skew.get(ist)));
        }

        List<String> names = new ArrayList<String>();
        List<List<Double>> cols = new ArrayList<List<Double>>();
        if (hasOrder(order, 1)) {
            names.add("QLen");
            cols.add(QLen);
        }
        if (hasOrder(order, 2)) {
            names.add("QLenVar");
            cols.add(QLenVar);
            names.add("QLenSCV");
            cols.add(QLenSCV);
        }
        if (hasOrder(order, 3)) {
            names.add("QLenM3");
            cols.add(QLenM3);
            names.add("QLenSkew");
            cols.add(QLenSkew);
        }

        NetworkMomentStationTable table = new NetworkMomentStationTable(names, cols);
        table.setOptions(this.options);
        table.setStationNames(stationName);
        table.setMoments(mom);
        return table;
    }

    /**
     * Higher moments of the per-chain queue length, up to the second moment.
     *
     * @return table of per-(Station, Chain) queue-length moments
     * @see #getMomentChainTable(int[])
     */
    public NetworkMomentChainTable getMomentChainTable() {
        return getMomentChainTable(momentOrderUpTo(2, 3));
    }

    /**
     * Higher moments of the per-chain queue length, up to the given moment order.
     *
     * @param order the highest moment order to report, an integer in 1..3
     * @return table of per-(Station, Chain) queue-length moments
     * @see #getMomentChainTable(int[])
     */
    public NetworkMomentChainTable getMomentChainTable(int order) {
        return getMomentChainTable(momentOrderUpTo(order, 3));
    }

    /**
     * Higher moments of the per-chain queue length.
     *
     * <p>Returns a table with one row per (Station, Chain) giving the moments of the
     * queue length of that chain at that station, Q_(i,c) = sum_(r in chain c)
     * n(i,r). This is the chain-level analogue of {@link #getAvgChainTable()}, and it
     * sits between the two other moment tables: {@link #getMomentTable(int[])} is per
     * class, {@link #getMomentStationTable(int[])} is per station total, and
     * this one is per chain, i.e. per group of classes that circulate together.</p>
     *
     * <p>Unlike the per-class table, order 3 is fully available here. All three
     * tables are the same recursion under different groupings of the classes: the
     * generating parameter scales the service times of a class subset T at a station,
     * and the moments it produces are those of sum_(r in T) n(i,r). T = {r} gives
     * {@link #getMomentTable(int[])}, T = chain gives this table, T = all classes
     * gives {@link #getMomentStationTable(int[])}. That is Theorem 1 of
     * Akyildiz and Strelen; Strelen's own x_i is the last case.</p>
     *
     * <p>The ALGORITHM is chosen by the solver's method, set at construction, not by
     * an argument here; see {@link #getMomentStationTable(int[])}. A
     * Linearizer-family method approximates the per-station totals only, so it cannot
     * express a per-chain grouping and is rejected here unless every class already
     * sits in one chain, in which case the chain IS the station total.</p>
     *
     * <p>Restricted to closed, single-server models, which is the scope of
     * {@code Pfqn_sens_mom}.</p>
     *
     * <p>Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
     * Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
     * Communications 39(6):828-832, 1991, Theorem 1; J. C. Strelen, "Moment Analysis
     * for Closed Queuing Networks and its Linearizer", Performance Evaluation
     * 11:127-142, 1990, equation (3.2).</p>
     *
     * @param order the set of moment orders to report; each must be an integer in
     *              1..3. Null or empty selects the default, up to the second moment.
     * @return table of per-(Station, Chain) queue-length moments
     * @see #getMomentTable(int[])
     * @see #getMomentStationTable(int[])
     */
    public NetworkMomentChainTable getMomentChainTable(int[] order) {
        if (order == null || order.length == 0) {
            order = momentOrderUpTo(2, 3);
        }
        order = validateMomentOrder(order, 3);
        // The algorithm is a property of the solver, not of this call; see
        // getMomentStationTable.
        String method = getOptions().method;

        NetworkStruct snl = model.getStruct();
        int R = snl.nclasses;
        Matrix N = snl.njobs;

        Ret.snGetProductFormParams pf = snGetProductFormParams(snl);
        Matrix D = pf.D;
        Matrix Np = pf.N;
        Matrix Z = pf.Z;
        Matrix Ssrv = pf.S;

        List<Integer> queueIndices = new ArrayList<Integer>();
        for (int i = 0; i < snl.nodetype.size(); i++) {
            if (snl.nodetype.get(i) == NodeType.Queue) {
                queueIndices.add(Integer.valueOf(i));
            }
        }
        int Mq = queueIndices.size();
        Matrix Ztot = sumColumns(Z, R);

        // See getMomentStationTable: the moment identity Cov = L dQ/dL is a theorem
        // about the product form, so outside it L dQ/dL is not a covariance and no
        // correct value exists to return.
        if (!snHasProductForm(snl)) {
            line_error(mfilename(new Object() {
            }), "getMomentChainTable requires a product-form model: the moment identity "
                    + "Cov[n,n] = L dQ/dL holds only under product form, so no correct value "
                    + "exists here. Use SolverCTMC (exact distribution) or SolverLDES with "
                    + "setReward for the moments of a non-product-form model.");
        }
        for (int r = 0; r < R; r++) {
            if (isInf(N.get(r))) {
                line_error(mfilename(new Object() {
                }), "getMomentChainTable supports closed models only. Per-class second moments "
                        + "of an open or mixed model are available from getMomentTable.");
            }
        }
        for (int i = 0; i < Ssrv.length(); i++) {
            if (Ssrv.get(i) > 1) {
                line_error(mfilename(new Object() {
                }), "getMomentChainTable supports single-server stations only: the higher-moment "
                        + "recursion of the reference is stated for load-independent stations.");
            }
        }

        // the chain of each class; sn.chains is (nchains x nclasses)
        int[] rawChain = new int[R];
        for (int r = 0; r < R; r++) {
            int c = -1;
            for (int j = 0; j < snl.chains.getNumRows(); j++) {
                if (snl.chains.get(j, r) != 0.0) {
                    c = j;
                    break;
                }
            }
            if (c < 0) {
                line_error(mfilename(new Object() {
                }), String.format("class %s belongs to no chain", snl.classnames.get(r)));
            }
            rawChain[r] = c;
        }
        // Pfqn_sens_mom requires the group labels to be consecutive from 1, so drop
        // any chain that holds no class rather than leaving a hole in the numbering
        TreeSet<Integer> usedSet = new TreeSet<Integer>();
        for (int r = 0; r < R; r++) {
            usedSet.add(Integer.valueOf(rawChain[r]));
        }
        int Cg = usedSet.size();
        int[] used = new int[Cg];
        int jj = 0;
        for (Integer v : usedSet) {
            used[jj++] = v.intValue();
        }
        Matrix groups = new Matrix(1, R);
        int[] grp = new int[R];
        for (int r = 0; r < R; r++) {
            int g = 0;
            while (used[g] != rawChain[r]) {
                g++;
            }
            grp[r] = g;
            groups.set(0, r, g + 1);
        }

        // See getMomentStationTable: the solver's method selects the algorithm. The
        // Linearizer approximates the per-station totals only, so it cannot express a
        // per-chain grouping unless every class already sits in one chain, in which
        // case the chain IS the station total.
        NetworkMomentStationResult mom = new NetworkMomentStationResult();
        if (isLinearizerMethod(method)) {
            if (Cg > 1) {
                line_error(mfilename(new Object() {
                }), String.format("the solver method '%s' approximates the per-station totals, "
                        + "so it cannot produce a per-chain grouping of %d chains. Use an exact "
                        + "method, or getMomentStationTable.", method, Cg));
            }
            mom.linearizer = pfqn_sens_linearizer(D, Np, Ztot);
            mom.m = mom.linearizer.m;
            mom.Var = mom.linearizer.Var;
            mom.Cov = mom.linearizer.Cov;
            mom.M2 = mom.linearizer.M2;
            mom.M3 = mom.linearizer.M3;
            mom.Skew = mom.linearizer.Skew;
            mom.CovAsym = mom.linearizer.CovAsym;
        } else if (isExactMvaMethod(method)) {
            mom.exact = pfqn_sens_mom(D, Np, Ztot, Matrix.ones(1, Mq), groups);
            mom.m = mom.exact.m;
            mom.Var = mom.exact.Var;
            mom.Cov = mom.exact.Cov;
            mom.CovG = mom.exact.CovG;
            mom.M2 = mom.exact.M2;
            mom.M3 = mom.exact.M3;
            mom.Skew = mom.exact.Skew;
            mom.CovAsym = mom.exact.CovAsym;
        } else {
            // See getMomentStationTable: this solver supplies the means and the
            // derivatives are taken numerically. Here the parameter scales the demands
            // of one CHAIN's classes at one station, which is the class subset T of
            // Akyildiz-Strelen Theorem 1, so the moments it generates are those of that
            // chain's queue length.
            mom = chainMomentsByFiniteDifference(snl, queueIndices, grp, Cg);
        }

        List<String> stationName = new ArrayList<String>();
        List<String> chainName = new ArrayList<String>();
        List<Double> QLen = new ArrayList<Double>();
        List<Double> QLenVar = new ArrayList<Double>();
        List<Double> QLenSCV = new ArrayList<Double>();
        List<Double> QLenM3 = new ArrayList<Double>();
        List<Double> QLenSkew = new ArrayList<Double>();
        for (int ist = 0; ist < Mq; ist++) {
            for (int g = 0; g < Cg; g++) {
                boolean visited = false;
                for (int r = 0; r < R; r++) {
                    if (grp[r] == g && D.get(ist, r) > 0) {
                        visited = true;
                        break;
                    }
                }
                if (!visited) {
                    continue;   // no class of this chain visits this station
                }
                stationName.add(snl.nodenames.get(queueIndices.get(ist).intValue()));
                chainName.add("Chain" + (used[g] + 1));
                QLen.add(Double.valueOf(mom.m.get(ist, g)));
                QLenVar.add(Double.valueOf(mom.Var.get(ist, g)));
                if (mom.m.get(ist, g) > 0) {
                    QLenSCV.add(Double.valueOf(mom.Var.get(ist, g)
                            / (mom.m.get(ist, g) * mom.m.get(ist, g))));
                } else {
                    QLenSCV.add(Double.valueOf(Double.NaN));
                }
                QLenM3.add(Double.valueOf(mom.M3.get(ist, g)));
                QLenSkew.add(Double.valueOf(mom.Skew.get(ist, g)));
            }
        }

        List<String> names = new ArrayList<String>();
        List<List<Double>> cols = new ArrayList<List<Double>>();
        if (hasOrder(order, 1)) {
            names.add("QLen");
            cols.add(QLen);
        }
        if (hasOrder(order, 2)) {
            names.add("QLenVar");
            cols.add(QLenVar);
            names.add("QLenSCV");
            cols.add(QLenSCV);
        }
        if (hasOrder(order, 3)) {
            names.add("QLenM3");
            cols.add(QLenM3);
            names.add("QLenSkew");
            cols.add(QLenSkew);
        }

        NetworkMomentChainTable table = new NetworkMomentChainTable(names, cols);
        table.setOptions(this.options);
        table.setStationNames(stationName);
        table.setChainNames(chainName);
        table.setMoments(mom);
        return table;
    }

    /**
     * True for the Linearizer family of solver methods. Everything else maps to the
     * exact recursion: the moment analysis has only these two algorithms, and the
     * exact one is the right default for a method with no approximate counterpart.
     *
     * @param method the solver's method token
     * @return true if the method selects the Linearizer approximation
     */
    private static boolean isLinearizerMethod(String method) {
        if (method == null) {
            return false;
        }
        return method.equalsIgnoreCase("lin")
                || method.equalsIgnoreCase("amva.lin")
                || method.equalsIgnoreCase("egflin")
                || method.equalsIgnoreCase("gflin");
    }

    /**
     * Methods whose means are the exact MVA recursion, so {@code Pfqn_sens_mom}'s
     * analytic derivatives apply directly.
     *
     * @param method the solver's method token
     * @return true if the method selects the exact recursion
     */
    private static boolean isExactMvaMethod(String method) {
        if (method == null) {
            return false;
        }
        return method.equalsIgnoreCase("default")
                || method.equalsIgnoreCase("mva")
                || method.equalsIgnoreCase("exact");
    }

    /**
     * Mean queue lengths of a perturbed structure, under this solver's own method.
     *
     * <p>This is the mean-value oracle behind the numerical-derivative path of
     * {@link #getMomentChainTable(int[])} and {@link #getMomentStationTable(int[])}.
     * The moment identity Cov[n,n] = L dQ/dL does not care HOW the mean queue lengths
     * were obtained, only that they are the means of a product-form model as a
     * function of its demands. So rather than hand-differentiating each algorithm, the
     * solver is re-run on a perturbed structure and differentiated numerically.</p>
     *
     * <p>Running THIS solver, rather than one chosen algorithm, is what makes the path
     * general: it covers every method of every solver, including the
     * normalizing-constant methods of {@code SolverNC} (comom, ca, le, mom, ...) and
     * the summation methods of {@code SolverMVA} (sum, esum), none of which any
     * hand-differentiated implementation reaches. Restricting the oracle to one
     * analyzer would restrict the moments to that analyzer's methods for no
     * mathematical reason.</p>
     *
     * <p>The perturbed sn is injected by copying the model and overwriting its cached
     * structure, then constructing a fresh solver of this class with this solver's
     * options. The copy matters: {@link Network} is a reference type, so writing sn on
     * the caller's model would corrupt it; {@code setHasStruct(true)} matters so the
     * copy does not regenerate the struct and discard the perturbation; and the fresh
     * solver matters because a solver caches its results and would otherwise return
     * the unperturbed answer.</p>
     *
     * @param sn the perturbed network structure
     * @return the mean queue lengths (stations x classes) under this solver's method
     */
    private Matrix solveMeansForStruct(NetworkStruct sn) {
        Network m2 = this.model.copy();
        m2.setStruct(sn);
        m2.setHasStruct(true);   // the copy must not regenerate and discard the perturbation
        try {
            // the equivalent of MATLAB's feval(class(self), model, self.getOptions):
            // every Network solver exposes a (Network, SolverOptions) constructor
            Constructor<? extends NetworkSolver> ctor =
                    this.getClass().getConstructor(Network.class, SolverOptions.class);
            NetworkSolver solver = ctor.newInstance(m2, getOptions());
            return solver.getAvgQLen();
        } catch (Exception e) {
            throw new RuntimeException("the numerical-derivative moment path could not "
                    + "re-run " + this.getClass().getSimpleName() + " as a mean-value "
                    + "oracle: it must expose a (Network, SolverOptions) constructor", e);
        }
    }

    /**
     * Moments of the per-station totals from any solver and method, by central
     * differences of that method's OWN mean queue lengths.
     *
     * <p>The identity does not care how the means were obtained, so this solver is
     * used as a mean-value oracle via {@link #solveMeansForStruct(NetworkStruct)}. The
     * parameter is y_i, a scaling of
     * station i's whole demand column, which is Strelen's x_i; at y = 1 the
     * y-derivatives are the scaled x-derivatives that (3.2) asks for. Because
     * D(i,r) = visits(i,r)/rate(i,r), scaling station i's rates by 1/f scales its
     * whole demand column by f, which is exactly the column perturbation wanted.</p>
     *
     * <p>Unlike MATLAB, whose NetworkStruct is a plain struct and therefore copies by
     * value, the JAR's {@link NetworkStruct} is a reference type, so the perturbation
     * is applied to a deep {@link jline.lang.Copyable#copy()} and the caller's sn is
     * left untouched.</p>
     *
     * <p>Cost: 2*M extra solves. Accuracy: the moments inherit the accuracy of the
     * method's means, and the second derivative inherits the usual h^2 truncation.</p>
     *
     * @param sn           the (unperturbed) network struct
     * @param queueIndices node indices of the queueing stations
     * @return the moment result built from numerically obtained derivatives
     */
    private NetworkMomentStationResult momentsByFiniteDifference(NetworkStruct sn,
            List<Integer> queueIndices) {
        int M = queueIndices.size();
        double h = 1e-4;
        int[] qst = new int[M];
        for (int i = 0; i < M; i++) {
            qst[i] = (int) sn.nodeToStation.get(queueIndices.get(i).intValue());
        }
        Matrix m0 = solveTotals(sn, qst);
        Matrix dm = new Matrix(M, M);
        Matrix d2m = new Matrix(M, 1);
        for (int hcol = 0; hcol < M; hcol++) {
            Matrix mp = solveTotals(scaleDemandColumn(sn, qst[hcol], 1 + h), qst);
            Matrix mm = solveTotals(scaleDemandColumn(sn, qst[hcol], 1 - h), qst);
            for (int i = 0; i < M; i++) {
                dm.set(i, hcol, (mp.get(i, 0) - mm.get(i, 0)) / (2 * h));
            }
            d2m.set(hcol, 0,
                    (mp.get(hcol, 0) - 2 * m0.get(hcol, 0) + mm.get(hcol, 0)) / (h * h));
        }
        return packFiniteDifference(m0, dm, d2m);
    }

    /**
     * Scales one station's whole demand column by {@code factor}, via its service
     * rates, on a deep copy of the struct.
     *
     * @param sn      the struct to perturb; not modified
     * @param station station index whose demand column is scaled
     * @param factor  the demand scaling factor
     * @return a copy of {@code sn} with the station's rates divided by {@code factor}
     */
    private static NetworkStruct scaleDemandColumn(NetworkStruct sn, int station,
            double factor) {
        NetworkStruct sn2 = sn.copy();
        for (int r = 0; r < sn.rates.getNumCols(); r++) {
            sn2.rates.set(station, r, sn.rates.get(station, r) / factor);
        }
        return sn2;
    }

    /**
     * Station totals under this solver's own method.
     *
     * @param sn  the struct to solve
     * @param qst station indices of the queueing stations
     * @return M x 1 vector of per-station total mean queue lengths
     */
    private Matrix solveTotals(NetworkStruct sn, int[] qst) {
        Matrix QN = solveMeansForStruct(sn);
        int M = qst.length;
        Matrix m = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double s = 0;
            for (int r = 0; r < QN.getNumCols(); r++) {
                s += QN.get(qst[i], r);
            }
            m.set(i, 0, s);
        }
        return m;
    }

    /**
     * Equation (3.2), applied to numerically obtained derivatives.
     *
     * @param m   M x 1 means
     * @param dm  M x M first derivatives
     * @param d2m M x 1 second derivatives
     * @return the assembled moment result
     */
    private static NetworkMomentStationResult packFiniteDifference(Matrix m, Matrix dm,
            Matrix d2m) {
        int M = m.getNumRows();
        NetworkMomentStationResult mom = new NetworkMomentStationResult();
        mom.m = m;
        Matrix Cov = new Matrix(M, M);
        double covAsym = 0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                Cov.set(i, j, (dm.get(i, j) + dm.get(j, i)) / 2);
                covAsym = Math.max(covAsym, Math.abs(dm.get(i, j) - dm.get(j, i)));
            }
        }
        mom.Cov = Cov;
        mom.CovAsym = covAsym;
        Matrix Var = new Matrix(M, 1);
        Matrix M2 = new Matrix(M, 1);
        Matrix M3 = new Matrix(M, 1);
        Matrix Skew = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double mi = m.get(i, 0);
            double d1 = dm.get(i, i);
            Var.set(i, 0, d1);
            M2.set(i, 0, d1 + mi * mi);
            double m3 = d2m.get(i, 0) + (1 + 3 * mi) * d1 + mi * mi * mi;
            M3.set(i, 0, m3);
            double mu3 = m3 - 3 * mi * M2.get(i, 0) + 2 * mi * mi * mi;
            if (d1 > 0) {
                Skew.set(i, 0, mu3 / Math.pow(d1, 1.5));
            } else {
                Skew.set(i, 0, Double.NaN);
            }
        }
        mom.Var = Var;
        mom.M2 = M2;
        mom.M3 = M3;
        mom.Skew = Skew;
        return mom;
    }

    /**
     * Per-chain moments from any solver and method, by central differences.
     *
     * <p>The parameter y_(i,g) scales the demands of group g's classes at station i,
     * which is the class-subset parameter T of Akyildiz-Strelen Theorem 1; the moments
     * it generates are those of Q_(i,g) = sum_(r in g) n(i,r). Since
     * D(i,r) = visits(i,r)/rate(i,r), the perturbation is applied to the rates of that
     * group's classes at that station alone, leaving the other classes' demands at the
     * same station untouched. That per-class granularity is what separates this from
     * {@link #momentsByFiniteDifference(NetworkStruct, List)}, which scales a whole
     * column.</p>
     *
     * <p>Cost: 2*M*Cg extra solves.</p>
     *
     * @param sn           the (unperturbed) network struct
     * @param queueIndices node indices of the queueing stations
     * @param grp          zero-based group label of each class
     * @param Cg           number of groups (chains) that hold at least one class
     * @return the moment result, with per-(station, group) entries
     */
    private NetworkMomentStationResult chainMomentsByFiniteDifference(NetworkStruct sn,
            List<Integer> queueIndices, int[] grp, int Cg) {
        int M = queueIndices.size();
        double h = 1e-4;
        int[] qst = new int[M];
        for (int i = 0; i < M; i++) {
            qst[i] = (int) sn.nodeToStation.get(queueIndices.get(i).intValue());
        }
        Matrix m0 = solveGroupTotals(sn, qst, grp, Cg);
        // dm is the 4-index derivative d m(i,g) / d y(hi,hg), held flat under the
        // bijection p(i,g) = i*Cg + g. Any consistent bijection gives the same
        // symmetrization and the same asymmetry residual, so the choice is free.
        Matrix dm = new Matrix(M * Cg, M * Cg);
        Matrix d2m = new Matrix(M, Cg);
        for (int hi = 0; hi < M; hi++) {
            for (int hg = 0; hg < Cg; hg++) {
                Matrix mp = solveGroupTotals(scaleGroupDemands(sn, qst[hi], grp, hg, 1 + h),
                        qst, grp, Cg);
                Matrix mm = solveGroupTotals(scaleGroupDemands(sn, qst[hi], grp, hg, 1 - h),
                        qst, grp, Cg);
                for (int i = 0; i < M; i++) {
                    for (int g = 0; g < Cg; g++) {
                        dm.set(i * Cg + g, hi * Cg + hg,
                                (mp.get(i, g) - mm.get(i, g)) / (2 * h));
                    }
                }
                d2m.set(hi, hg,
                        (mp.get(hi, hg) - 2 * m0.get(hi, hg) + mm.get(hi, hg)) / (h * h));
            }
        }
        return packChainFiniteDifference(m0, dm, d2m, Cg);
    }

    /**
     * Scales the demands of group g's classes at one station by {@code factor}, via
     * rates, on a deep copy of the struct.
     *
     * @param sn      the struct to perturb; not modified
     * @param station station index whose rates are scaled
     * @param grp     zero-based group label of each class
     * @param g       the group whose classes are scaled
     * @param factor  the demand scaling factor
     * @return a copy of {@code sn} with that group's rates at that station scaled
     */
    private static NetworkStruct scaleGroupDemands(NetworkStruct sn, int station, int[] grp,
            int g, double factor) {
        NetworkStruct sn2 = sn.copy();
        for (int r = 0; r < grp.length; r++) {
            if (grp[r] == g) {
                sn2.rates.set(station, r, sn.rates.get(station, r) / factor);
            }
        }
        return sn2;
    }

    /**
     * Per-(station, group) mean queue lengths under the solver's own method.
     *
     * @param sn  the struct to solve
     * @param qst station indices of the queueing stations
     * @param grp zero-based group label of each class
     * @param Cg  number of groups
     * @return M x Cg matrix of per-(station, group) mean queue lengths
     */
    private Matrix solveGroupTotals(NetworkStruct sn, int[] qst, int[] grp, int Cg) {
        Matrix QN = solveMeansForStruct(sn);
        int M = qst.length;
        Matrix mg = new Matrix(M, Cg);
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < Cg; g++) {
                double s = 0;
                for (int r = 0; r < grp.length; r++) {
                    if (grp[r] == g) {
                        s += QN.get(qst[i], r);
                    }
                }
                mg.set(i, g, s);
            }
        }
        return mg;
    }

    /**
     * Equation (3.2), applied to numerically obtained derivatives, per (station, group).
     *
     * @param m   M x Cg means
     * @param dm  (M*Cg) x (M*Cg) flat first derivatives under p(i,g) = i*Cg + g
     * @param d2m M x Cg second derivatives
     * @param Cg  number of groups
     * @return the assembled moment result
     */
    private static NetworkMomentStationResult packChainFiniteDifference(Matrix m, Matrix dm,
            Matrix d2m, int Cg) {
        int M = m.getNumRows();
        int P = M * Cg;
        NetworkMomentStationResult mom = new NetworkMomentStationResult();
        mom.m = m;
        Matrix sym = new Matrix(P, P);
        double covAsym = 0;
        for (int a = 0; a < P; a++) {
            for (int b = 0; b < P; b++) {
                sym.set(a, b, (dm.get(a, b) + dm.get(b, a)) / 2);
                covAsym = Math.max(covAsym, Math.abs(dm.get(a, b) - dm.get(b, a)));
            }
        }
        mom.CovAsym = covAsym;
        if (Cg == 1) {
            // one group per station: the chain IS the station total, so the covariance
            // collapses to the M x M cross-station form the station table carries
            mom.Cov = sym;
        } else {
            // more than one group: cross-chain and cross-station covariances both
            // exist, so only the general (M x Cg x M x Cg) form can hold them
            Matrix[][] CovG = new Matrix[M][Cg];
            for (int i = 0; i < M; i++) {
                for (int g = 0; g < Cg; g++) {
                    Matrix blk = new Matrix(M, Cg);
                    for (int j = 0; j < M; j++) {
                        for (int g2 = 0; g2 < Cg; g2++) {
                            blk.set(j, g2, sym.get(i * Cg + g, j * Cg + g2));
                        }
                    }
                    CovG[i][g] = blk;
                }
            }
            mom.CovG = CovG;
        }
        Matrix Var = new Matrix(M, Cg);
        Matrix M2 = new Matrix(M, Cg);
        Matrix M3 = new Matrix(M, Cg);
        Matrix Skew = new Matrix(M, Cg);
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < Cg; g++) {
                double mig = m.get(i, g);
                double d1 = dm.get(i * Cg + g, i * Cg + g);
                Var.set(i, g, d1);
                M2.set(i, g, d1 + mig * mig);
                double m3 = d2m.get(i, g) + (1 + 3 * mig) * d1 + mig * mig * mig;
                M3.set(i, g, m3);
                double mu3 = m3 - 3 * mig * M2.get(i, g) + 2 * mig * mig * mig;
                if (d1 > 0) {
                    Skew.set(i, g, mu3 / Math.pow(d1, 1.5));
                } else {
                    Skew.set(i, g, Double.NaN);
                }
            }
        }
        mom.Var = Var;
        mom.M2 = M2;
        mom.M3 = M3;
        mom.Skew = Skew;
        return mom;
    }

    /**
     * The set 1..order, the meaning a scalar moment order carries: "everything up to
     * order k", not "order k alone".
     */
    private static int[] momentOrderUpTo(int order, int maxorder) {
        if (order < 1 || order > maxorder) {
            line_error(mfilename(new Object() {
            }), String.format("order must be an integer in 1..%d, or a vector of such integers.",
                    maxorder));
        }
        int[] out = new int[order];
        for (int k = 0; k < order; k++) {
            out[k] = k + 1;
        }
        return out;
    }

    /**
     * Validates a set of moment orders and returns it sorted and deduplicated.
     *
     * <p>A set of more than one order is taken literally. A ONE-ELEMENT set takes the
     * scalar path instead, i.e. {@code {k}} expands to 1..k: MATLAB's
     * {@code isscalar([2])} is true, so {@code getMomentTable([2])} means "up to the
     * second moment" there, and the Python port reproduces it. Diverging here would
     * make the same argument mean different things in different codebases.</p>
     */
    private static int[] validateMomentOrder(int[] order, int maxorder) {
        if (order.length == 1) {
            return momentOrderUpTo(order[0], maxorder);
        }
        TreeSet<Integer> set = new TreeSet<Integer>();
        for (int k = 0; k < order.length; k++) {
            if (order[k] < 1 || order[k] > maxorder) {
                line_error(mfilename(new Object() {
                }), String.format("order must be an integer in 1..%d, or a vector of such "
                        + "integers.", maxorder));
            }
            set.add(Integer.valueOf(order[k]));
        }
        if (set.isEmpty()) {
            line_error(mfilename(new Object() {
            }), String.format("order must be an integer in 1..%d, or a vector of such integers.",
                    maxorder));
        }
        int[] out = new int[set.size()];
        int j = 0;
        for (Integer v : set) {
            out[j++] = v.intValue();
        }
        return out;
    }

    /**
     * Whether a moment order belongs to a validated (sorted) order set.
     */
    private static boolean hasOrder(int[] order, int k) {
        for (int j = 0; j < order.length; j++) {
            if (order[j] == k) {
                return true;
            }
        }
        return false;
    }

    /**
     * Column sums of a matrix, as a 1 x R row.
     */
    private static Matrix sumColumns(Matrix A, int R) {
        Matrix out = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double acc = 0.0;
            for (int i = 0; i < A.getNumRows(); i++) {
                acc += A.get(i, r);
            }
            out.set(0, r, acc);
        }
        return out;
    }

    /**
     * Squared coefficient of variation. NaN when the mean is zero, since the SCV is
     * then undefined rather than infinite in any useful sense.
     */
    private static double momentScv(double variance, double meanv) {
        if (Double.isNaN(variance) || meanv <= 0) {
            return Double.NaN;
        }
        return variance / (meanv * meanv);
    }

    /**
     * An FCFS station in a BCMP network must serve every class at the same
     * exponential rate; that is the precondition for reading a per-visit rate off
     * it. Non-FCFS stations are unconstrained here, see the caller.
     */
    private static boolean fcfsRatesAreClassIndependent(NetworkStruct snl,
                                                        List<Integer> queueIndices, int R) {
        for (int ist = 0; ist < queueIndices.size(); ist++) {
            int sIdx = (int) snl.nodeToStation.get(queueIndices.get(ist).intValue());
            if (snl.sched.get(snl.stations.get(sIdx)) != SchedStrategy.FCFS) {
                continue;
            }
            double ref = -1;
            for (int r = 0; r < R; r++) {
                double rate = snl.rates.get(sIdx, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    continue;
                }
                if (ref < 0) {
                    ref = rate;
                } else if (Math.abs(rate - ref) > GlobalConstants.FineTol * Math.max(1, ref)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * The common service time of a station, i.e. the reciprocal of the
     * class-independent rate. Zero if no class is served here.
     */
    private static double serviceTimeOf(NetworkStruct snl, int sIdx, int R) {
        for (int r = 0; r < R; r++) {
            double rate = snl.rates.get(sIdx, r);
            if (Double.isFinite(rate) && rate > 0) {
                return 1 / rate;
            }
        }
        return 0;
    }

    /**
     * The throughput of class r: the arrival rate for an open class, the computed
     * throughput of the queue-length recursion for a closed one.
     */
    private static double throughputOf(NetworkMomentResult mom, Matrix lambda, Matrix N, int r) {
        if (isInf(N.get(r))) {
            return lambda.get(r);
        } else if (mom.hasQlen()) {
            return mom.X.get(r);
        }
        return 0;
    }

    /**
     * Alias for getAvgTable().
     *
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT() {
        return getAvgTable();
    }

    /**
     * Alias for getAvgTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    /**
     * Alias for getAvgTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    /**
     * Alias for getAvgTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // ===== Table -> T shorthand aliases for the auxiliary result tables =====
    // (moment/sensitivity/cache/item/orbit), mirroring aT/avgT/getAvgT.

    // getMomentTable
    public NetworkMomentTable momentT() { return getMomentTable(); }
    public NetworkMomentTable momentT(int order) { return getMomentTable(order); }
    public NetworkMomentTable momentT(int[] order) { return getMomentTable(order); }
    public NetworkMomentTable mT() { return getMomentTable(); }
    public NetworkMomentTable mT(int order) { return getMomentTable(order); }
    public NetworkMomentTable mT(int[] order) { return getMomentTable(order); }
    public NetworkMomentTable getMomentT() { return getMomentTable(); }
    public NetworkMomentTable getMomentT(int order) { return getMomentTable(order); }
    public NetworkMomentTable getMomentT(int[] order) { return getMomentTable(order); }

    // getMomentChainTable
    public NetworkMomentChainTable momentChainT() { return getMomentChainTable(); }
    public NetworkMomentChainTable momentChainT(int order) { return getMomentChainTable(order); }
    public NetworkMomentChainTable momentChainT(int[] order) { return getMomentChainTable(order); }
    public NetworkMomentChainTable mCT() { return getMomentChainTable(); }
    public NetworkMomentChainTable mCT(int order) { return getMomentChainTable(order); }
    public NetworkMomentChainTable mCT(int[] order) { return getMomentChainTable(order); }
    public NetworkMomentChainTable getMomentChainT() { return getMomentChainTable(); }
    public NetworkMomentChainTable getMomentChainT(int order) { return getMomentChainTable(order); }
    public NetworkMomentChainTable getMomentChainT(int[] order) { return getMomentChainTable(order); }

    // getMomentStationTable
    public NetworkMomentStationTable momentStationT() { return getMomentStationTable(); }
    public NetworkMomentStationTable momentStationT(int order) { return getMomentStationTable(order); }
    public NetworkMomentStationTable momentStationT(int[] order) { return getMomentStationTable(order); }
    public NetworkMomentStationTable mST() { return getMomentStationTable(); }
    public NetworkMomentStationTable mST(int order) { return getMomentStationTable(order); }
    public NetworkMomentStationTable mST(int[] order) { return getMomentStationTable(order); }
    public NetworkMomentStationTable getMomentStationT() { return getMomentStationTable(); }
    public NetworkMomentStationTable getMomentStationT(int order) { return getMomentStationTable(order); }
    public NetworkMomentStationTable getMomentStationT(int[] order) { return getMomentStationTable(order); }

    // getSensitivityTable
    public NetworkSensitivityTable sensitivityT() { return getSensitivityTable(); }
    public NetworkSensitivityTable sensitivityT(String method, double step, String scheme) { return getSensitivityTable(method, step, scheme); }
    public NetworkSensitivityTable sT() { return getSensitivityTable(); }
    public NetworkSensitivityTable sT(String method, double step, String scheme) { return getSensitivityTable(method, step, scheme); }
    public NetworkSensitivityTable getSensitivityT() { return getSensitivityTable(); }
    public NetworkSensitivityTable getSensitivityT(String method, double step, String scheme) { return getSensitivityTable(method, step, scheme); }

    // getAvgCacheTable
    public NetworkAvgCacheTable cacheAvgT() { return getAvgCacheTable(); }
    public NetworkAvgCacheTable aCaT() { return getAvgCacheTable(); }
    public NetworkAvgCacheTable getAvgCacheT() { return getAvgCacheTable(); }

    // getAvgItemTable
    public NetworkAvgItemTable itemAvgT() { return getAvgItemTable(); }
    public NetworkAvgItemTable aIT() { return getAvgItemTable(); }
    public NetworkAvgItemTable getAvgItemT() { return getAvgItemTable(); }

    // getAvgOrbitTable
    public NetworkAvgOrbitTable orbitAvgT() { return getAvgOrbitTable(); }
    public NetworkAvgOrbitTable aOT() { return getAvgOrbitTable(); }
    public NetworkAvgOrbitTable getAvgOrbitT() { return getAvgOrbitTable(); }

    // getAvgLossTable
    public NetworkLossTable lossAvgT() { return getAvgLossTable(); }
    public NetworkLossTable aLT() { return getAvgLossTable(); }
    public NetworkLossTable getAvgLossT() { return getAvgLossTable(); }

    // getAvgRegionLossTable
    public NetworkLossTable regionLossAvgT() { return getAvgRegionLossTable(); }
    public NetworkLossTable aRLT() { return getAvgRegionLossTable(); }
    public NetworkLossTable getAvgRegionLossT() { return getAvgRegionLossTable(); }

    /**
     * Returns a table of deadline-related metrics (response time and tardiness) organized by station and job class.
     *
     * @return table containing response time and tardiness metrics, or null if tardiness data is not available
     */
    public NetworkAvgTable getDeadlineTable() {
        this.sn = model.getStruct(true);

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (Double.isFinite(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            line_error(mfilename(new Object() {
            }), "Unable to compute results and therefore unable to get deadline table: " + e.getMessage());
            return null;
        }

        if (this.result == null || this.result.TardN == null || this.result.SysTardN == null) {
            return null;
        }

        Matrix RN = this.result.RN;
        Matrix TardN = this.result.TardN;
        Matrix SysTardN = this.result.SysTardN;

        if (TardN.isEmpty() || SysTardN.isEmpty()) {
            return null;
        }

        List<Double> Rval = new ArrayList<>();
        List<Double> Tardval = new ArrayList<>();
        List<Double> SysTardval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (RN.get(i, k) > 0 || TardN.get(i, k) > 0 || SysTardN.get(0, k) > 0) {
                    Rval.add(RN.get(i, k));
                    Tardval.add(TardN.get(i, k));
                    SysTardval.add(SysTardN.get(0, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
        }

        NetworkAvgTable deadlineTable = new NetworkAvgTable(
                new ArrayList<>(), new ArrayList<>(), Rval,
                new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                Tardval, SysTardval);
        deadlineTable.setOptions(this.options);
        deadlineTable.setClassNames(className);
        deadlineTable.setStationNames(stationName);
        return deadlineTable;
    }

    /**
     * Computes and returns average throughputs at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average throughputs [stations x classes]
     */
    /**
     * Third-party libraries this solver will use on the given model, named once
     * per session by {@link #showLibraryAttribution}.
     *
     * <p>LINE bundles ports of BUTools, Q-MAM, MAMSolver, rmf_tool and others
     * under {@code jline/lib}, and credits them so a user knows whose
     * algorithms produced the numbers. This is distinct from
     * {@code InputOutput.line_ack}, which credits an external TOOL a wrapper
     * solver shells out to (JMT, LQNS, qnsolver). Solvers that use no bundled
     * library inherit the empty list.
     *
     * @param sn the model structure, may be null
     * @param options the solver options, may be null
     * @return the library names, possibly empty
     */
    public List<String> getLibrariesUsed(NetworkStruct sn, SolverOptions options) {
        return new ArrayList<String>();
    }

    /**
     * Bibliographic references for the algorithms this solver used.
     *
     * <p>The references follow what the run actually did: the method the
     * analyzer resolved to (not merely the one requested), the fork-join
     * transformation if the model has forks, and the percentile method of the
     * last getPerctRespT call. Mirrors MATLAB @NetworkSolver/citations.m and
     * python NetworkSolver.citations.
     *
     * <p>The bibliography key on each entry is internal shorthand for
     * doc/latex/biblio.bib and is not meant to be displayed.
     *
     * @return one entry per distinct paper
     */
    public List<LineCitations.Citation> citations() {
        String family = "";
        String cls = this.getClass().getSimpleName();
        if (cls.equals("SolverCTMC")) { family = "ctmc"; }
        else if (cls.equals("SolverSSA")) { family = "ssa"; }
        else if (cls.equals("SolverFluid") || cls.equals("SolverFLD")) { family = "fld"; }
        else if (cls.equals("SolverLN")) { family = "ln"; }
        else if (cls.equals("SolverJMT")) { family = "jmt"; }
        else if (cls.equals("SolverLQNS")) { family = "lqns"; }
        else if (cls.equals("SolverENV")) { family = "env"; }
        else if (cls.equals("SolverMVA")) { family = "mva"; }
        else if (cls.equals("SolverNC")) { family = "nc"; }
        else if (cls.equals("SolverMAM")) { family = "mam"; }
        else if (cls.equals("SolverBA")) { family = "ba"; }

        List<String> tokens = new ArrayList<String>();
        if (family.equals("ctmc") || family.equals("ssa") || family.equals("fld")
                || family.equals("ln") || family.equals("jmt") || family.equals("lqns")
                || family.equals("env")) {
            tokens.add(family);
        }
        if (this.options != null && this.options.method != null) {
            addMethodToken(tokens, family, this.options.method);
        }
        if (this.result != null && this.result.method != null) {
            // the reported name is 'default/<actual>' when dispatch chose
            for (String part : this.result.method.split("/")) {
                addMethodToken(tokens, family, part);
            }
        }
        try {
            if (this.model != null && this.model.hasFork()) {
                String fjm = (this.options != null && this.options.config != null
                        && this.options.config.fork_join != null)
                        ? this.options.config.fork_join.toLowerCase() : "mmt";
                if (fjm.equals("default") || fjm.equals("fjt") || fjm.isEmpty()) {
                    fjm = "mmt";
                } else if (fjm.equals("heidelberger-trivedi")) {
                    fjm = "ht";
                }
                tokens.add(fjm);
            }
        } catch (Exception e) {
            // a model that cannot report its topology contributes no token
        }
        if (this.lastPerctMethod != null && !this.lastPerctMethod.isEmpty()) {
            tokens.add(this.lastPerctMethod);
        }
        return LineCitations.citationsFor(tokens);
    }

    private static void addMethodToken(List<String> tokens, String family, String method) {
        if (method == null) {
            return;
        }
        String m = method.trim().toLowerCase();
        if (m.isEmpty() || m.equals("default")) {
            return;
        }
        tokens.add(family.isEmpty() ? m : family + "." + m);
    }

    /**
     * Third-party libraries this solver will use on its own model, without
     * printing anything.
     *
     * <p>Attribution in LINE is pull-based, in the spirit of Sage's
     * {@code sage.misc.citation.get_systems}: nothing is written to the console
     * during a solve, and the user asks for the list when citing. The wrapper
     * acknowledgement ({@code InputOutput.line_ack}) follows the same rule: it
     * prints only at {@link VerboseLevel#DEBUG}.
     *
     * @return the library names, possibly empty
     */
    public List<String> libraries() {
        NetworkStruct sn = null;
        try {
            sn = this.model.getStruct(false);
        } catch (Exception e) {
            sn = null;
        }
        return getLibrariesUsed(sn, this.options);
    }

    /**
     * Prints the library attribution once per session, honouring the verbosity.
     * Not called automatically: see {@link #libraries()}.
     *
     * @param sn the model structure, may be null
     * @param options the solver options, may be null
     */
    public void showLibraryAttribution(NetworkStruct sn, SolverOptions options) {
        if (options != null && options.verbose == VerboseLevel.SILENT) {
            return;
        }
        if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
            return;
        }
        if (jline.lang.GlobalConstants.isLibraryAttributionShown()) {
            return;
        }
        List<String> libs = getLibrariesUsed(sn, options);
        if (libs == null || libs.isEmpty()) {
            return;
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < libs.size(); i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(libs.get(i));
        }
        System.out.println("The solver will leverage " + sb + ".");
        jline.lang.GlobalConstants.setLibraryAttributionShown(true);
    }

    /**
     * Predicts fork-join request response time percentiles with the ForkTail
     * approximation.
     *
     * <p>ForkTail needs only the mean and variance of the per-branch task
     * response times, so it applies to heterogeneous branches and mixed service
     * laws. The branch arrival rate is read from the solved throughputs and the
     * service moments from the node distributions; see
     * {@link jline.api.fj.FJ_tail_forktail}.
     *
     * <p>The topology gate is strict, because the approximation is defined on a
     * per-branch M/G/1 task response time: exactly one fork with a matching
     * join, and every branch a single station feeding that join.
     *
     * <p>Being a heavy-traffic result, it under-predicts at low load, the more
     * so the more variable the service; a warning is issued when the busiest
     * branch is below 0.5 utilization.
     *
     * @param percentiles percentiles, fractions in (0,1) or percentages in (0,100)
     * @param method extraction method, only "forktail" is supported here
     * @return matrix [classes x percentiles] of predicted request response times,
     *         with NaN rows for classes that do not traverse the fork
     */
    /** Percentile extraction method of the last getPerctRespT call, for citations(). */
    protected String lastPerctMethod = "";

    public Matrix getPerctRespT(double[] percentiles, String method) {
        this.lastPerctMethod = method == null ? "" : method.toLowerCase();
        if (method == null || !"forktail".equalsIgnoreCase(method)) {
            throw new IllegalArgumentException(
                "Unsupported percentile method '" + method + "'; only 'forktail' is available here.");
        }
        NetworkStruct sn = this.model.getStruct(false);

        java.util.List<Integer> forks = new java.util.ArrayList<Integer>();
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == jline.lang.constant.NodeType.Fork) {
                forks.add(i);
            }
        }
        if (forks.isEmpty()) {
            throw new RuntimeException("The forktail method requires a model with a Fork node.");
        }
        if (forks.size() > 1) {
            throw new RuntimeException(
                "The forktail method supports a single fork-join pair; this model has " + forks.size() + " forks.");
        }
        int f = forks.get(0);
        int joinIdx = -1;
        for (int j = 0; j < sn.nnodes; j++) {
            if (sn.fj != null && sn.fj.get(f, j) > 0) {
                joinIdx = j;
                break;
            }
        }
        if (joinIdx < 0) {
            throw new RuntimeException("The fork node has no matching join; the request response time is undefined.");
        }

        java.util.List<Integer> branches = new java.util.ArrayList<Integer>();
        for (int j = 0; j < sn.nnodes; j++) {
            if (sn.connmatrix.get(f, j) > 0) {
                if (sn.isstation.get(j) == 0) {
                    throw new RuntimeException("Branch node " + sn.nodenames.get(j)
                        + " is not a station; the forktail method needs one queueing station per branch.");
                }
                if (sn.connmatrix.get(j, joinIdx) == 0) {
                    throw new RuntimeException("Branch station " + sn.nodenames.get(j)
                        + " does not feed the join directly; the forktail method needs one station per branch.");
                }
                branches.add(j);
            }
        }

        Matrix TN = this.getAvgTput();
        Matrix UN = this.getAvgUtil();
        int K = sn.nclasses;
        Matrix out = new Matrix(K, percentiles.length);
        for (int r = 0; r < K; r++) {
            double[] ET = new double[branches.size()];
            double[] VT = new double[branches.size()];
            double maxrho = 0.0;
            boolean traverses = true;
            for (int bi = 0; bi < branches.size(); bi++) {
                int ist = (int) sn.nodeToStation.get(branches.get(bi));
                double lambda = TN.get(ist, r);
                if (lambda <= GlobalConstants.FineTol) {
                    traverses = false;
                    break;
                }
                Distribution svc = ((jline.lang.nodes.ServiceStation) this.model.getNodes().get(branches.get(bi)))
                    .getServiceProcess(this.model.getClassByIndex(r));
                double ES = svc.getMean();
                double VS = svc.getVar();
                double ES2 = VS + ES * ES;
                double ES3 = svc.getSkewness() * Math.pow(VS, 1.5) + 3 * ES * ES2 - 2 * ES * ES * ES;
                FJ_tail_forktail.ResptMoments m = FJ_tail_forktail.fj_mg1_respt_moments(lambda, ES, ES2, ES3);
                ET[bi] = m.mean;
                VT[bi] = m.variance;
                maxrho = Math.max(maxrho, UN.get(ist, r));
            }
            for (int pi = 0; pi < percentiles.length; pi++) {
                out.set(r, pi, traverses
                    ? FJ_tail_forktail.fj_tail_forktail(ET, VT, percentiles[pi])
                    : Double.NaN);
            }
            if (traverses && maxrho < 0.5) {
                line_warning(mfilename(new Object(){}), String.format(
                    "ForkTail is a heavy-traffic approximation; the busiest branch is at utilization %.2f, "
                    + "so the tail is likely under-predicted.\n", maxrho));
            }
        }
        return out;
    }

    public Matrix getAvgTput() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.TN;
    }

    /**
     * Returns average throughputs aggregated by job chains.
     *
     * @return matrix of throughputs [stations x chains]
     */
    public Matrix getAvgTputChain() {
        int C = sn.nchains;
        Matrix TNclass = getAvgTput();
        Matrix TN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    TN.set(i, c, TN.get(i, c) + TNclass.get(i, k));
                }
            }
        }
        return TN;
    }

    /**
     * Returns the average throughput metric handles.
     *
     * @return handles for throughput metrics
     */
    public AvgHandle getAvgTputHandles() {
        return this.avgHandles.T;
    }

    /**
     * Computes and returns average server utilizations at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average utilizations [stations x classes]
     */
    public Matrix getAvgUtil() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.UN;
    }

    /**
     * Returns average server utilizations aggregated by job chains.
     *
     * @return matrix of utilizations [stations x chains]
     */
    public Matrix getAvgUtilChain() {
        int C = sn.nchains;
        Matrix UNclass = getAvgUtil();
        Matrix UN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    UN.set(i, c, UN.get(i, c) + UNclass.get(i, k));
                }
            }
        }
        return UN;
    }

    /**
     * Returns the average utilization metric handles.
     *
     * @return handles for utilization metrics
     */
    public AvgHandle getAvgUtilHandles() {
        return this.avgHandles.U;
    }

    /**
     * Computes and returns average waiting times in queue excluding service time.
     * Waiting time = Response time - Service time (1/rate)
     *
     * @return matrix of average waiting times [stations x classes]
     */
    public Matrix getAvgWaitT() {
        if (!this.hasResults()) {
            this.getAvg();
        }

        Matrix RN = this.result.RN;
        if (RN == null || RN.isEmpty()) {
            return new Matrix(0, 0);
        }

        // Get service rates for waiting time calculation
        this.sn = model.getStruct(true);
        Matrix WT = RN.copy();

        // Calculate waiting time = response time - service time (1/rate)
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                double serviceTime = 1.0 / sn.rates.get(i, k);
                WT.set(i, k, RN.get(i, k) - serviceTime);

                // Set waiting time to 0 for source nodes (they don't have queues)
                if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
                    WT.set(i, k, 0.0);
                }

                // Ensure non-negative waiting times (numerical safety)
                if (WT.get(i, k) < 0) {
                    WT.set(i, k, 0.0);
                }
            }
        }

        return WT;
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfPassT(AvgHandle R) {
        throw new RuntimeException("getCdfPassT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfPassT() {
        return getCdfPassT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * Uses an exponential approximation based on average response times.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for response times [stations x classes]
     */
    public DistributionResult getCdfRespT(AvgHandle R) {
        NetworkStruct sn = this.sn;
        if (sn == null) {
            sn = model.getStruct(false);
        }
        int M = sn.nstations;
        int K = sn.nclasses;

        DistributionResult dr = new DistributionResult(M, K, "response_time");

        if (!this.hasAvgResults()) {
            getAvg();
        }

        if (this.result == null || this.result.RN == null) {
            return dr;
        }

        int n = 100;
        for (int i = 0; i < M; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) != NodeType.Source) {
                for (int c = 0; c < K; c++) {
                    double respT = this.result.RN.get(i, c);
                    if (Double.isFinite(respT) && respT > 0) {
                        double lambda = 1.0 / respT;
                        Matrix cdf = new Matrix(n, 2);
                        for (int q = 0; q < n; q++) {
                            double quantile = 0.001 + (0.999 - 0.001) * q / (n - 1);
                            cdf.set(q, 0, quantile);
                            cdf.set(q, 1, -Math.log(1.0 - quantile) / lambda);
                        }
                        dr.setCdf(i, c, cdf);
                    } else {
                        Matrix cdf = new Matrix(1, 2);
                        cdf.set(0, 0, 1.0);
                        cdf.set(0, 1, 0.0);
                        dr.setCdf(i, c, cdf);
                    }
                }
            }
        }
        return dr;
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for response times [stations x classes]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfRespT() {
        return getCdfRespT(getAvgRespTHandles());
    }

    /**
     * Returns the queueing network model being solved.
     *
     * @return the network model
     */
    public Network getModel() {
        return model;
    }

    /**
     * Warm-start the solver from the steady-state solution of an auxiliary
     * solver. The auxiliary solver's steady-state distribution decides an
     * integer job placement (see {@link WarmStart#warmStartPlacement}): with
     * SolverCTMC the mode of the exact aggregate stationary distribution, with
     * any other solver the rounded mean queue lengths conserving each
     * closed-class population. The placement is applied as the model initial
     * state via initFromMarginal, which the state-driven solvers honor:
     * SolverFluid starts the ODE integration from it, SolverSSA starts the
     * simulated trajectory from it, and SolverJMT preloads the stations with
     * it. Note that this modifies the initial state of the model object shared
     * with any other solver instance.
     *
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @return this solver, for chaining
     */
    public NetworkSolver initFromSolver(NetworkSolver initSolver) {
        NetworkStruct snModel = this.model.getStruct(true);
        Matrix placement = WarmStart.warmStartPlacement(initSolver, snModel);
        this.model.initFromMarginal(placement);
        // Re-sync the solver's cached struct: getStruct(true) re-populates
        // sn.state from the node states just assigned by initFromMarginal,
        // which analyzers reading this.sn (e.g. SolverSSA) would otherwise
        // miss having snapshotted the default state at construction time.
        this.sn = this.model.getStruct(true);
        return this;
    }

    /**
     * Sets the queueing network model to be solved.
     *
     * @param model the network model to set
     */
    public void setModel(Network model) {
        this.model = model;
    }

    /**
     * Returns marginal state probabilities for a specific node and state.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node  the node index for which to compute probabilities
     * @param state the state vector to query (optional, null for all states)
     * @return result containing marginal state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProb(int node, Matrix state) {
        throw new RuntimeException("getProb is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns marginal state probabilities for a specific node (all states).
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node the node index for which to compute probabilities
     * @return result containing marginal state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProb(int node) {
        return getProb(node, null);
    }

    /**
     * Probability of a SPECIFIC per-class job distribution at a station.
     * Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
     *
     * <p>Compare with {@link #getProbMarg}: returns queue-length distribution for a
     * single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).</p>
     *
     * @param node    the node index for which to compute probabilities
     * @param state_a per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
     * @return scalar probability in [0,1]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbAggr(int node, Matrix state_a) {
        throw new RuntimeException("getProbAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Probability of a SPECIFIC per-class job distribution at a station (current state).
     * Returns P(n1 jobs of class 1, n2 jobs of class 2, ...).
     *
     * <p>Compare with {@link #getProbMarg}: returns queue-length distribution for a
     * single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).</p>
     *
     * @param node the node index for which to compute probabilities
     * @return scalar probability in [0,1]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbAggr(int node) {
        return getProbAggr(node, null);
    }

    /**
     * Probability distribution for queue length of a SINGLE class at a station.
     * Returns P(n jobs of class r) for n=0,1,...,N(r).
     *
     * <p>Compare with {@link #getProbAggr}: returns probability of a specific per-class
     * distribution, e.g., P(2 class-1, 1 class-2) as a scalar.</p>
     *
     * @param node     the node index for which to compute probabilities
     * @param jobclass the job class index for marginalization
     * @param state_m  specific states to query, or null for all
     * @return vector where element n+1 = P(n jobs of this class)
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbMarg(int node, int jobclass, Matrix state_m) {
        throw new RuntimeException("getProbMarg is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Probability distribution for queue length of a SINGLE class at a station (all states).
     * Returns P(n jobs of class r) for n=0,1,...,N(r).
     *
     * <p>Compare with {@link #getProbAggr}: returns probability of a specific per-class
     * distribution, e.g., P(2 class-1, 1 class-2) as a scalar.</p>
     *
     * @param node     the node index for which to compute probabilities
     * @param jobclass the job class index for marginalization
     * @return vector where element n+1 = P(n jobs of this class)
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbMarg(int node, int jobclass) {
        return getProbMarg(node, jobclass, null);
    }

    /**
     * Returns the logarithm of the normalizing constant of state probabilities.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing the log normalizing constant
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbNormConstAggr() {
        throw new RuntimeException("getProbNormConstAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns joint state probabilities for the entire system.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing joint state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbSys() {
        throw new RuntimeException("getProbSys is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns aggregated joint state probabilities for the entire system.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing aggregated joint state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbSysAggr() {
        throw new RuntimeException("getProbSysAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns a table of average stage metrics organized by job classes.
     * For non-environment models, this returns the same as getAvgTable()
     * since there is only one implicit stage.
     *
     * @return table containing stage-level metrics for each class
     */
    public NetworkAvgTable getStageTable() {
        return getAvgTable();
    }

    /**
     * Returns a table of average stage metrics organized by job classes with keepDisabled option.
     * For non-environment models, this returns the same as getAvgTable(keepDisabled)
     * since there is only one implicit stage.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing stage-level metrics for each class
     */
    public NetworkAvgTable getStageTable(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    /**
     * Computes transient average station metrics over the specified time interval.
     * The timespan is defined in the solver options.
     */
    public void getTranAvg() {

        this.tranHandles = model.getTranHandles();

        // NOTE: This was in LINE, but I believe is legacy as 'matrix' method can provide tran results
    /*    if (!Objects.equals(options.method, "default")) {
      System.err.println(
          "getTranAvg is not offered by the specified method. Setting the solution method to \"closing\".");
      resetResults();
    }
    options.method = "closing";*/

        sn = model.getStruct(true);
        double minRate = sn.rates.elementMin();
        if (!hasTranResults()) {
            if (isInf(options.timespan[0]) && isInf(options.timespan[1])) {
                options.timespan[0] = 0;
                options.timespan[1] = 30 / minRate;
                if (options.verbose != VerboseLevel.SILENT) {
                    System.out.format(
                            "Timespan of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                            options.timespan[1]);
                }
            } else if (isInf(options.timespan[0])) {
                options.timespan[0] = 0;
                if (options.verbose != VerboseLevel.SILENT) {
                    System.out.format(
                            "Start time of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                            options.timespan[1]);
                }
            } else if (isInf(options.timespan[1])) {
                options.timespan[1] = 30 / minRate;
                if (options.verbose != VerboseLevel.SILENT) {
                    System.out.format(
                            "End time of transient analysis unspecified, setting the timespan option to [%f, %f].\n",
                            options.timespan[0], options.timespan[1]);
                }
            }
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                line_error(mfilename(new Object(){}), "IllegalAccessException upon running runAnalyzer()");
            } catch (ParserConfigurationException e) {
                line_error(mfilename(new Object(){}), "ParserConfigurationException upon running runAnalyzer()");
            } catch (IOException e) {
                line_error(mfilename(new Object(){}), "IOException upon running runAnalyzer()");
            }
        }

        // Note: Metrics storage functionality not yet implemented
        // This would involve storing computed transient metrics for later retrieval
        // Implementation depends on the specific metrics storage system to be defined
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing transient CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfPassT(AvgHandle R) {
        throw new RuntimeException("getTranCdfPassT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * Uses default response time handles.
     *
     * @return result containing transient CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfPassT() {
        return getTranCdfPassT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing transient CDFs for response times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfRespT(AvgHandle R) {
        throw new RuntimeException("getTranCdfRespT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * Uses default response time handles.
     *
     * @return result containing transient CDFs for response times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfRespT() {
        return getTranCdfRespT(getAvgRespTHandles());
    }

    /**
     * Returns the transient performance metric handles.
     *
     * @return the transient handles object
     */
    public SolverTranHandles getTranHandles() {
        return this.tranHandles;
    }

    /**
     * Sets the transient performance metric handles.
     *
     * @param handles the transient handles to set
     */
    public void setTranHandles(SolverTranHandles handles) {
        this.tranHandles = handles;
    }

    /**
     * Checks if the solver has computed steady-state average metrics.
     *
     * @return true if steady-state results are available, false otherwise
     */
    protected boolean hasAvgResults() {
        return !(result == null) && (!((result.QN == null || result.QN.isEmpty()) &&
                (result.UN == null || result.UN.isEmpty()) &&
                (result.RN == null || result.RN.isEmpty()) &&
                (result.TN == null || result.TN.isEmpty()) &&
                (result.CN == null || result.CN.isEmpty()) &&
                (result.XN == null || result.XN.isEmpty())));
    }

    /**
     * Checks if the solver has computed steady-state distribution metrics.
     *
     * @return true if distribution results are available
     */
    public boolean hasDistribResults() {
        // Check if any distribution methods have been implemented and can return results
        try {
            DistributionResult cdfResult = getCdfRespT();
            return cdfResult != null && cdfResult.cdfData != null && !cdfResult.cdfData.isEmpty();
        } catch (RuntimeException e) {
            // Distribution methods not implemented for this solver
            return false;
        }
    }

    /**
     * Helper method to check if a matrix array contains NaN values.
     *
     * @param matrices 2D array of matrices to check
     * @return true if any matrix contains NaN values, false otherwise
     */
    private boolean hasNaNValues(Matrix[][] matrices) {
        if (matrices == null) {
            return false;
        }
        for (Matrix[] matrixRow : matrices) {
            if (matrixRow == null) {
                continue;
            }
            for (Matrix matrix : matrixRow) {
                if (matrix != null && matrix.hasNaN()) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Checks if the solver has computed transient average metrics.
     *
     * @return true if transient results are available, false otherwise
     */
    public boolean hasTranResults() {
        if (this.hasResults() && result.QNt != null) {
            return result.QNt.length > 0 && result.QNt[0].length > 0
                    && result.QNt[0][0] != null && !result.QNt[0][0].isEmpty();
        }
        return false;
    }

    /**
     * Initializes performance metric handles from the model.
     * This method retrieves and sets both average and transient handles.
     */
    protected void initHandles() {
        this.avgHandles = model.getAvgHandles();
        this.tranHandles = model.getTranHandles();
        // Force model to refresh struct if needed
        this.sn = model.getStruct(true);
    }

    // Model
    public Network model() {
        return getModel();
    }

    /**
     * Detailed print function that displays all contents of the NetworkSolver.
     * This function prints actual values, not reference addresses, and includes ALL fields
     * from both NetworkSolver and its parent Solver class for comprehensive comparison.
     */
    public void print() {
        // NetworkSolver fields
        System.out.println("name: " + (this.name != null ? "\"" + this.name + "\"" : "null"));
        System.out.println("model: " + (this.model != null ? "\"" + this.model.getName() + "\"" : "null"));
        System.out.println("enableChecks: " + this.enableChecks);
        System.out.println("options: " + (this.options != null ? this.options.toString() : "null"));
        System.out.println("result: " + (this.result != null ? "SolverResult[" + this.result.getClass().getSimpleName() + "]" : "null"));
        System.out.println("avgHandles: " + (this.avgHandles != null ? "SolverAvgHandles[" + this.avgHandles.getClass().getSimpleName() + "]" : "null"));
        System.out.println("tranHandles: " + (this.tranHandles != null ? "SolverTranHandles[" + this.tranHandles.getClass().getSimpleName() + "]" : "null"));

        // NetworkStruct (sn) information
        if (this.sn == null) {
            System.out.println("sn: null");
        } else {
            this.sn.print();
        }
    }

    // Probabilities
    public ProbabilityResult prob(int node) {
        return getProb(node);
    }

    public ProbabilityResult prob(int node, Matrix state_a) {
        return getProb(node, state_a);
    }

    public ProbabilityResult probAggr(int node) {
        return getProbAggr(node);
    }

    public ProbabilityResult probAggr(int node, Matrix state_a) {
        return getProbAggr(node, state_a);
    }

    public ProbabilityResult probMarg(int node, int jobclass) {
        return getProbMarg(node, jobclass);
    }

    public ProbabilityResult probMarg(int node, int jobclass, Matrix state_m) {
        return getProbMarg(node, jobclass, state_m);
    }

    public ProbabilityResult probNormConstAggr() {
        return getProbNormConstAggr();
    }

    public ProbabilityResult probSys() {
        return getProbSys();
    }

    public ProbabilityResult probSysAggr() {
        return getProbSysAggr();
    }

    /**
     * Validates model compatibility and method support before analysis.
     *
     * @param options solver options containing method specification
     * @throws RuntimeException if model contains unsupported features or method is invalid
     */
    public void runAnalyzerChecks(SolverOptions options) {
        // Propagate solver verbose level to global so that model-level
        // messages (e.g., priority info in refreshStruct) respect it
        if (options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        // Basic model validation - check for empty model
        if (model == null) {
            throw new RuntimeException("Model cannot be null");
        }

        if (model.getNumberOfNodes() == 0) {
            throw new RuntimeException("Model must contain at least one node");
        }

        // Basic options validation
        if (options == null) {
            throw new RuntimeException("SolverOptions cannot be null");
        }

        // Additional model structure validation
        NetworkStruct sn = model.getStruct(false);
        if (sn == null) {
            throw new RuntimeException("Unable to obtain model structure");
        }

        // Check for finite timespan values if applicable
        if (options.timespan != null && options.timespan.length >= 2) {
            if (options.timespan[0] < 0 || (Double.isFinite(options.timespan[1]) && options.timespan[1] <= options.timespan[0])) {
                throw new RuntimeException("Invalid timespan configuration: start time must be non-negative and end time must be greater than start time");
            }
        }
    }

    /**
     * Samples state trajectories for a specific node.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node      the node index to sample from
     * @param numEvents the number of events to sample
     * @return result containing sampled state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sample(int node, int numEvents) {
        throw new RuntimeException("sample is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples aggregated state trajectories for a specific node.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node      the node index to sample from
     * @param numEvents the number of events to sample
     * @return result containing sampled aggregated state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleAggr(int node, int numEvents) {
        throw new RuntimeException("sampleAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples joint system state trajectories.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param numEvents the number of events to sample
     * @return result containing sampled joint system state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleSys(int numEvents) {
        throw new RuntimeException("sampleSys is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples aggregated joint system state trajectories.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param numEvents the number of events to sample
     * @return result containing sampled aggregated joint system state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleSysAggr(int numEvents) {
        throw new RuntimeException("sampleSysAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Stores computed average metrics at steady-state in the solver result.
     *
     * @param Q       queue length matrix [stations x classes]
     * @param U       utilization matrix [stations x classes]
     * @param R       response time matrix [stations x classes]
     * @param T       throughput matrix [stations x classes]
     * @param A       arrival rate matrix [stations x classes]
     * @param W       residence time matrix [stations x classes]
     * @param C       system response time vector [chains]
     * @param X       system throughput vector [chains]
     * @param runtime computation time in seconds
     * @param method  solution algorithm used
     * @param iter    number of iterations performed
     */
    public void setAvgResults(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix A, Matrix W, Matrix C, Matrix X,
                              double runtime, String method, int iter) {
        this.result.solver = this.getName();
        this.result.QN = Q.copy();
        this.result.UN = U.copy();
        this.result.RN = R.copy();
        this.result.TN = T.copy();
        this.result.AN = A.copy();
        this.result.WN = W.copy();
        this.result.CN = C.copy();
        this.result.XN = X.copy();
        this.result.runtime = runtime;
        this.result.method = method;
        // The iteration count was previously used only for the printout below,
        // leaving SolverResult.iter at 0 for API callers. MATLAB stores it as
        // result.Avg.iter; store it here too so the count is readable, not just
        // printable -- it is what distinguishes a converged AMVA solve from one
        // that exhausted options.iter_max.
        this.result.iter = iter;

        if (this.options.verbose != VerboseLevel.SILENT) {
            if (iter <= 1) {
                System.out.printf(
                        "%s analysis [method: %s, lang: %s, env: %s] completed in %fs.\n",
                        this.name.replaceFirst("^Solver", ""),   // solver name with prefix stripped
                        this.result.method,                      // algorithm/method
                        "java",                                  // language label
                        System.getProperty("java.version"),      // actual JVM version in use
                        this.result.runtime                      // elapsed time in seconds
                );
            } else {
                System.out.printf(
                        "%s analysis [method: %s, lang: %s, env: %s] completed in %fs. Iterations: %d.\n",
                        this.name.replaceFirst("^Solver", ""),   // solver name with prefix stripped
                        this.result.method,                      // algorithm/method
                        "java",                                  // language label
                        System.getProperty("java.version"),      // actual JVM version in use
                        this.result.runtime,                     // elapsed time in seconds
                        iter                                     // iteration count
                );
            }
            System.out.flush();
        }
    }

    /**
     * Stores distribution metrics at steady-state.
     *
     * @param RD      distribution data [stations x classes] containing CDFs
     * @param runtime computation time for distribution analysis
     */
    protected final void setDistribResults(Matrix RD, double runtime) {
        this.result.solver = this.getName();
        //this.result.RD = RD.copy();
        this.result.runtime = runtime;
    }

    /**
     * Sets the language for solver execution.
     * This method configures the solver to use either MATLAB or Java backend.
     */
    protected void setLang() {
        // Java implementation always uses Java backend
        // This method is provided for consistency with MATLAB interface
    }

    /**
     * Stores computed transient average metrics.
     *
     * @param Qt       transient queue length matrices [time][stations x classes]
     * @param Ut       transient utilization matrices [time][stations x classes]
     * @param Rt       transient response time matrices [time][stations x classes]
     * @param Tt       transient throughput matrices [time][stations x classes]
     * @param Ct       transient system response time matrices [time][chains]
     * @param Xt       transient system throughput matrices [time][chains]
     * @param runtimet computation time for transient analysis
     */
    public void setTranAvgResults(Matrix[][] Qt, Matrix[][] Ut, Matrix[][] Rt, Matrix[][] Tt, Matrix[][]
            Ct, Matrix[][] Xt, double runtimet) {
        this.result.solver = getName();
        this.result.method = getOptions().method;
        // NaN values in transient metrics are normal (e.g., throughput at time 0
        // when no jobs exist at a station), so store the arrays as-is
        this.result.QNt = Qt != null ? Qt.clone() : new Matrix[0][0];
        this.result.UNt = Ut != null ? Ut.clone() : new Matrix[0][0];
        this.result.RNt = Rt != null ? Rt.clone() : new Matrix[0][0];
        this.result.TNt = Tt != null ? Tt.clone() : new Matrix[0][0];
        this.result.CNt = Ct != null ? Ct.clone() : new Matrix[0][0];
        this.result.XNt = Xt != null ? Xt.clone() : new Matrix[0][0];
        this.result.runtime = runtimet;
    }

    /**
     * Stores transient probability distributions.
     *
     * @param t        time points vector
     * @param pi_t     transient probability distributions
     * @param SS       steady-state distribution
     * @param runtimet computation time for transient analysis
     */
    protected final void setTranProb(Matrix t, Matrix pi_t, Matrix SS, double runtimet) {
        this.result.solver = getName();
        this.result.method = getOptions().method;
        this.result.t = t.copy();
        this.result.pi_t = pi_t.copy();
        this.result.SS = SS.copy();
        this.result.runtime = runtimet;
    }

    // Stage table
    public Object stageTable() {
        return getStageTable();
    }

    // stageTable -> stageT aliases
    public Object stageT() {
        return stageTable();
    }

    public NetworkAvgTable getStageT() {
        return getStageTable();
    }

    public NetworkAvgTable getStageT(boolean keepDisabled) {
        return getStageTable(keepDisabled);
    }

    // Transient analysis
    public void tranAvg() {
        getTranAvg();
    }

    public DistributionResult tranCdfPassT() {
        return getTranCdfPassT();
    }

    public DistributionResult tranCdfPassT(AvgHandle R) {
        return getTranCdfPassT(R);
    }

    public DistributionResult tranCdfRespT() {
        return getTranCdfRespT();
    }

    public DistributionResult tranCdfRespT(AvgHandle R) {
        return getTranCdfRespT(R);
    }

    public SolverTranHandles tranHandles() {
        return getTranHandles();
    }

    /**
     * Shared structural gate for finite station capacity (setCapacity) and
     * finite per-class buffers (setClassCap), used by the product-form solvers
     * (MVA, NC). A product-form solver has no representation of a finite
     * buffer, so without this gate it silently returns the UNCONSTRAINED
     * answer (e.g. QLen=4 instead of the M/M/1/2 value 0.8525). There is no
     * registry feature name for plain capacity, hence the structural test;
     * this mirrors the MATLAB NetworkSolver.checkBindingCapacity and the
     * native-Python check in solvers/solver_mva/solver_mva.py.
     * <p>
     * The test reads the node-level capacity set by the user, NOT sn.cap /
     * sn.classcap: refreshCapacity derives a FINITE sn.classcap (the chain
     * population) for every closed model, so an sn-level test would reject
     * every closed model.
     * <p>
     * Only a capacity that can actually BIND is rejected. A closed model whose
     * station capacity is at least the total population can never block a job,
     * so the declaration is a no-op and the product-form answer stays exact (a
     * common idiom: setCap(N) on an order-independent station of an N-job
     * closed model). njobs is Inf for an open class, so any finite capacity
     * reachable by an open class binds.
     * <p>
     * Cache models are exempt: Cache builds retrieval queues that legitimately
     * carry a per-class capacity of 1, and MVA/NC solve those through their
     * dedicated cache/retrieval analyzers rather than as a buffer constraint.
     *
     * @param model      - the network model
     * @param sn         - the network structure (for njobs)
     * @param solverName - solver name used in the returned message
     * @return null if no capacity can bind, otherwise the reason naming the
     *         offending station
     */
    public static String bindingCapacityReason(Network model, NetworkStruct sn, String solverName) {
        List<Node> nodes = model.getNodes();
        for (Node node : nodes) {
            if (node instanceof Cache) {
                return null;
            }
        }
        List<JobClass> jobClasses = model.getClasses();
        double totalJobs = 0; // Inf as soon as one class is open
        for (int r = 0; r < sn.nclasses; r++) {
            totalJobs += sn.njobs.get(0, r);
        }
        for (Node node : nodes) {
            if (!(node instanceof Station) || node instanceof Source || node instanceof Sink) {
                continue;
            }
            Station station = (Station) node;
            // hasFiniteCap() decodes the three "unbounded" encodings (MAX_VALUE, Inf, and
            // JMT2LINE's negative sentinel) in one place; see Station.hasFiniteCap.
            if (station.hasFiniteCap() && station.getCap() < totalJobs) {
                return "Finite station capacity (setCapacity=" + (long) station.getCap() + ") at station '"
                        + station.getName() + "' is not supported by " + solverName
                        + ". Use SolverCTMC, SolverJMT or SolverLDES.";
            }
            for (int r = 0; r < Math.min(jobClasses.size(), sn.nclasses); r++) {
                JobClass jobClass = jobClasses.get(r);
                double classCap = station.getClassCap(jobClass);
                if (classCap > 0 && classCap < Integer.MAX_VALUE && classCap < sn.njobs.get(0, r)) {
                    return "Finite per-class capacity (setClassCap=" + (long) classCap + " for class '"
                            + jobClass.getName() + "') at station '" + station.getName()
                            + "' is not supported by " + solverName
                            + ". Use SolverCTMC, SolverJMT or SolverLDES.";
                }
            }
        }
        return null;
    }

}
