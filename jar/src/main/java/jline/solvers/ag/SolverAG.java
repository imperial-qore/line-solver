/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.FeatureSet;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.ag.analyzers.Solver_ag_analyzer;
import jline.util.matrix.Matrix;

import jline.solvers.AvgHandle;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.api.sn.SnGetResidTFromRespT.snGetResidTFromRespT;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Agent-based (RCAT) solver.
 *
 * <p>Solves a network by the Reversed Compound Agent Theorem: every
 * (station, class) pair becomes an isolated agent, and the agents are coupled
 * ONLY through the reversed rates of the synchronizing actions. Agent k carries
 *
 * <pre>Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c</pre>
 *
 * and publishes, for every action it is active on, a scalar read off its own
 * stationary vector. The fixed point over that scalar vector is the whole
 * analysis.</p>
 *
 * <p>Because the coupling is that thin, the sweep parallelises exactly rather
 * than approximately: see {@link AgExec} for the {@code threads} and
 * {@code cluster} execution backends and why they walk the same iterates as the
 * serial loop.</p>
 *
 * <p>Methods: {@code inap}, {@code inapplus}, {@code inapinf}, and the vestigial
 * {@code exact} alias which warns and falls back to {@code inap}. Per Marin,
 * Rota Bulo and Balsamo, "A Numerical Algorithm for the Decomposition of
 * Cooperating Structured Markov Processes", MASCOTS 2012.</p>
 *
 * <p>These methods used to live in SolverMAM. They are the only algorithms in
 * LINE that read {@code sn.issignal}, so the G-network feature names belong here
 * and to no other solver.</p>
 */
public class SolverAG extends NetworkSolver {

    private static final String[] RCAT_METHODS =
            new String[]{"default", "inap", "inapplus", "inapinf", "exact"};

    public SolverAG(Network model) {
        super(model, "SolverAG", new AGOptions());
        this.result = new AGResult();
    }

    public SolverAG(Network model, String method) {
        super(model, "SolverAG", defaultOptions().method(method));
        this.result = new AGResult();
    }

    public SolverAG(Network model, Object... varargin) {
        this(model, defaultOptions());
        Solver.parseOptions(this.options, varargin);
        this.result = new AGResult();
    }

    public SolverAG(Network model, SolverOptions options) {
        super(model, "SolverAG", options);
        this.result = new AGResult();
    }

    public static AGOptions defaultOptions() {
        return new AGOptions();
    }

    public NetworkStruct getStruct() {
        return model.getStruct(true);
    }

    public List<String> listValidMethods() {
        return new ArrayList<String>(Arrays.asList(RCAT_METHODS));
    }

    public List<String> listValidMethods(Network model) {
        return listValidMethods();
    }

    /**
     * The RCAT feature envelope. The G-network names live here because
     * {@link jline.solvers.ag.handlers.Solver_ag_build} is the only code in LINE
     * that reads {@code sn.issignal}.
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "Fork", "Join", "Forker", "Joiner",
                "Delay", "DelayStation", "Queue",
                "APH", "Coxian", "Erlang", "Exp", "HyperExp", "MAP", "MMPP2",
                "Det", "Gamma", "Lognormal", "Pareto", "Uniform", "Weibull",
                "StatelessClassSwitcher", "InfiniteServer",
                "SharedServer", "Buffer", "Dispatcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_FCFS",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ClosedClass",
                "OpenClass",
                "OpenSignal", "ClosedSignal",
                "SignalType_NEGATIVE", "SignalType_CATASTROPHE",
                "SignalBatchRemoval"
        });
        return featSupported;
    }

    /**
     * Every AG method is the same decomposition differing only in how the
     * reversed rate is read off an agent, so they share one envelope; the genuine
     * restrictions are structural and applied in
     * {@link #supportsModelMethod(String)}.
     */
    public static FeatureSet methodFeatureSet(String method) {
        return SolverAG.getFeatureSet();
    }

    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        return SolverAG.methodFeatureSet(method);
    }

    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        return FeatureSet.supports(SolverAG.getFeatureSet(), featUsed);
    }

    /**
     * True when the RCAT analyzers can give a process of this type a phase
     * dimension: after snNonmarkovToPh (which AG runs with phfit='ph' and
     * preserveDet=false) it must hold a genuine (D0,D1) pair with non-negative
     * off-diagonal rates and a single arrival per epoch.
     *
     * <p>The list is an ALLOW-list on purpose: a process type nobody has checked
     * against this construction must be refused, not answered. Compare BY NAME
     * across codebases, never by the raw ordinal.</p>
     */
    public static boolean rcatSupportsProcess(ProcessType pt) {
        if (pt == null) {
            return true;   // nothing to serve
        }
        return pt == ProcessType.EXP || pt == ProcessType.ERLANG || pt == ProcessType.HYPEREXP
                || pt == ProcessType.PH || pt == ProcessType.APH || pt == ProcessType.COXIAN
                || pt == ProcessType.COX2 || pt == ProcessType.MAP || pt == ProcessType.MMPP2
                || pt == ProcessType.DET || pt == ProcessType.UNIFORM || pt == ProcessType.GAMMA
                || pt == ProcessType.PARETO || pt == ProcessType.WEIBULL
                || pt == ProcessType.LOGNORMAL || pt == ProcessType.REPLAYER
                || pt == ProcessType.IMMEDIATE || pt == ProcessType.DISABLED;
    }

    @Override
    public String supportsModelMethod(String method) {
        NetworkStruct sn = getStruct();
        if (sn.procid == null) {
            return "";
        }
        for (int i = 0; i < sn.nstations; i++) {
            Station st = sn.stations.get(i);
            Map<JobClass, ProcessType> byClass = sn.procid.get(st);
            if (byClass == null) {
                continue;
            }
            boolean isSource = st instanceof Source;
            for (int r = 0; r < sn.nclasses; r++) {
                JobClass jc = sn.jobclasses.get(r);
                ProcessType pt = byClass.get(jc);
                boolean isSignal = sn.issignal != null && sn.issignal.get(r, 0) > 0;
                // A signal is a trigger with no service, so its service entry is
                // never read; only its Source arrival rate is.
                if (isSignal && !isSource) {
                    continue;
                }
                double rate = sn.rates.get(i, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    continue;
                }
                if (isSignal) {
                    if (pt == null || pt == ProcessType.EXP) {
                        continue;
                    }
                    return "The " + method + " method needs an exponential signal arrival "
                            + "process (a removal signal is folded into the agent as a scalar "
                            + "rate), but station " + (i + 1) + " class " + (r + 1) + " is "
                            + ProcessType.toText(pt) + ". Use SolverMAM (method 'dec.source') "
                            + "for such models.";
                }
                if (rcatSupportsProcess(pt)) {
                    continue;
                }
                return "The " + method + " method supports processes with a Markovian "
                        + "(D0,D1) representation only (RCAT builds a phase dimension per "
                        + "agent out of it), but station " + (i + 1) + " class " + (r + 1)
                        + " is " + ProcessType.toText(pt) + ". Use SolverMAM (method "
                        + "'dec.source') for such models.";
            }
        }
        // RCAT models every station single-server; see _kb/06-solver-catalog.md
        for (int i = 0; i < sn.nstations; i++) {
            double c = sn.nservers.get(i);
            if (Double.isFinite(c) && c > 1) {
                return "The " + method + " method supports single-server stations only "
                        + "(RCAT does not model sn.nservers, so a multiserver station is "
                        + "driven at rho = lambda/mu instead of lambda/(c*mu)), but station "
                        + (i + 1) + " has " + (int) c + " servers. Use SolverMAM (method "
                        + "'dec.source') for multiserver models.";
            }
        }
        // A FINITE BUFFER IS NOT SOMETHING RCAT CAN CARRY, and it was being
        // answered rather than refused: nothing under solvers/ag reads sn.cap or
        // sn.classcap, so a capped station was decomposed as an unbounded one and
        // the table reported the UNCONSTRAINED figures (a closed 2-job tandem with
        // cap 1 on the second queue returned the same numbers with and without the
        // cap, 1.09 jobs in a buffer of 1). Lowering the component's level bound
        // (nlev = njobs(r)+1) to the buffer would not fix it: the top-level
        // boundary is a self-loop, which LOSES the arrival, whereas a closed job
        // refused at a full buffer must BLOCK the upstream departure, and that
        // coupling is exactly the independence RCAT assumes.
        // see _kb/06-solver-catalog.md
        String capReason = NetworkSolver.bindingCapacityReason(this.model, sn, "SolverAG");
        if (capReason != null) {
            return capReason;
        }
        return "";
    }

    @Override
    public void runAnalyzer() {
        if (this.options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        if (this.enableChecks) {
            String reason = this.supportsModelMethod(this.options.method);
            if (!reason.isEmpty()) {
                line_error(mfilename(new Object() {
                }), "This model contains features not supported by the solver. " + reason);
                return;
            }
        }

        // Finite Capacity Region: the RCAT decomposition has no aggregate
        // per-region job limit and would silently return the unconstrained answer.
        if (this.model.getStruct(false).nregions > 0) {
            throw new RuntimeException("This model uses a Finite Capacity Region (addRegion), "
                    + "which is not supported by SolverAG. Use SolverCTMC, SolverJMT, "
                    + "SolverSSA or SolverLDES, or setCapacity for a single-station limit.");
        }

        double start = System.nanoTime();
        NetworkStruct sn = getStruct();

        AGResult res = Solver_ag_analyzer.solver_ag_analyzer(sn, this.options);

        AvgHandle T = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(sn, res.TN, T);
        AvgHandle W = getAvgResidTHandles();
        Matrix WN = snGetResidTFromRespT(sn, res.RN, W);
        double finish = System.nanoTime();
        setAvgResults(res.QN, res.UN, res.RN, res.TN, AN, WN, res.CN, res.XN,
                (finish - start) / 1.0e9, res.method, res.iter);

        AGResult stored = new AGResult();
        stored.QN = res.QN;
        stored.UN = res.UN;
        stored.RN = res.RN;
        stored.TN = res.TN;
        stored.CN = res.CN;
        stored.XN = res.XN;
        stored.iter = res.iter;
        stored.method = res.method;
        stored.actionRates = res.actionRates;
        stored.equilibrium = res.equilibrium;
        stored.generators = res.generators;
        stored.runtime = res.runtime;
        this.result = stored;
    }

    /**
     * The converged reversed rates and the agents they induce, in addition to the
     * mean measures returned by getAvg.
     *
     * <p>Mean values alone hide the objects the decomposition is built on. The
     * reversed rates ARE the coupling between agents -- one scalar per
     * synchronizing action -- so a run cannot be checked against a published
     * derivation, or against a product-form condition, without them.</p>
     */
    public AGResult getAGResult() {
        if (!(this.result instanceof AGResult)) {
            runAnalyzer();
        }
        return (AGResult) this.result;
    }

    public static SolverType solverType() {
        return SolverType.AG;
    }
}
