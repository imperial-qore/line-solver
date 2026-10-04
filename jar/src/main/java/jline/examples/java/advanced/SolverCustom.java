/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

import java.util.Arrays;
import java.util.List;

/**
 * Template for a user-written LINE solver.
 *
 * <p>The MATLAB twin is a {@code @SolverCustom} class folder holding
 * SolverCustom.m and runAnalyzer.m; the two free functions the analyzer calls
 * stay in {@link SolverCustomAlgorithm} and {@link SolverCustomAnalyzer},
 * exactly as they do there.</p>
 *
 * <p>To write a real solver: fill in
 * {@link SolverCustomAlgorithm#solver_custom}, narrow {@link #getFeatureSet()}
 * to the features the algorithm actually supports, and list the method names it
 * accepts in {@link #listValidMethods()}.</p>
 */
public class SolverCustom extends NetworkSolver {

    public SolverCustom(Network model) {
        this(model, defaultOptions());
    }

    public SolverCustom(Network model, SolverOptions options) {
        super(model, "SolverCustom", options);
        // every concrete solver owns its result container, as SolverMVA does:
        // the base class declares the field but never allocates it
        this.result = new SolverResult();
    }

    /** The data structure summarizing the model (no initial state needed). */
    public NetworkStruct getStruct() {
        return this.model.getStruct(false);
    }

    /** The method names this solver accepts. */
    public List<String> listValidMethods() {
        return Arrays.asList("default");
    }

    /**
     * The model features this solver claims to support.
     *
     * <p>Kept identical to the MATLAB template's list, which is deliberately
     * broad: narrow it to what the algorithm can really answer, or
     * {@link #supports(Network)} will accept a model that
     * {@link SolverCustomAlgorithm#solver_custom} then gets wrong.</p>
     *
     * @return the declared feature set
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source", "Router", "ClassSwitch", "DelayStation", "Queue",
                "Fork", "Join", "Forker", "Joiner", "Logger",
                "Coxian", "Cox2", "APH", "Erlang", "Exp", "HyperExp", "Det",
                "Gamma", "Lognormal", "MAP", "MMPP2", "Normal", "PH", "Pareto",
                "Weibull", "Replayer", "Uniform",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer",
                "Buffer", "Dispatcher", "Server", "JobSink", "RandomSource",
                "ServiceTunnel", "LogTunnel", "Linkage",
                "Enabling", "Timing", "Firing", "Storage", "Place", "Transition",
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_DPS",
                "SchedStrategy_FCFS", "SchedStrategy_GPS", "SchedStrategy_SIRO",
                "SchedStrategy_HOL", "SchedStrategy_LCFS", "SchedStrategy_LCFSPR",
                "SchedStrategy_SEPT", "SchedStrategy_LEPT", "SchedStrategy_SJF",
                "SchedStrategy_LJF",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "RoutingStrategy_RROBIN", "RoutingStrategy_WRROBIN",
                "RoutingStrategy_SQ", "SchedStrategy_EXT",
                "ClosedClass", "OpenClass"});
        return featSupported;
    }

    /**
     * Whether this solver declares every feature the model uses.
     *
     * @param model the model to check
     * @return true if the model falls inside {@link #getFeatureSet()}
     */
    @Override
    public boolean supports(Network model) {
        return FeatureSet.supports(getFeatureSet(), model.getUsedLangFeatures());
    }

    /** The options this solver starts from. */
    public static SolverOptions defaultOptions() {
        SolverOptions options = Solver.defaultOptions();
        options.timespan = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        return options;
    }

    /** Runs the solver: read the structure, call the analyzer, publish results. */
    @Override
    public void runAnalyzer() {
        NetworkStruct sn = getStruct();
        SolverResult res = SolverCustomAnalyzer.solver_custom_analyzer(sn, this.options);
        // AN and WN have no separate template slot: the reference returns six
        // matrices, so the arrival rate and waiting time follow TN and RN.
        setAvgResults(res.QN, res.UN, res.RN, res.TN, res.TN, res.RN, res.CN, res.XN,
                res.runtime, this.options.method == null ? "default" : this.options.method, 1);
    }

    /** Builds a two-station closed model and runs the template over it. */
    public static void main(String[] args) {
        Network model = new Network("customExample");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", 2, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));

        SolverCustom solver = new SolverCustom(model);
        solver.runAnalyzer();
        Matrix QN = solver.result.QN;
        System.out.printf("QN is %dx%d and all zeros: %s%n",
                QN.getNumRows(), QN.getNumCols(), QN.elementMax() == 0);
    }
}
