/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers.fluid;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.aoi.*;
import jline.api.sn.SnNonmarkovToPh;
import jline.io.Ret;
import jline.lang.FeatureSet;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.processes.Coxian;
import jline.lang.state.FromMarginal;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.*;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.fluid.analyzers.ClosingAndStateDepMethodsAnalyzer;
import jline.solvers.fluid.analyzers.DiffusionAnalyzer;
import jline.solvers.fluid.analyzers.DaeAnalyzer;
import jline.solvers.fluid.analyzers.FluidAnalyzer;
import jline.solvers.fluid.analyzers.KoPenderAnalyzer;
import jline.solvers.fluid.analyzers.MFQAnalyzer;
import jline.solvers.fluid.analyzers.QsysLimitAnalyzer;
import jline.solvers.fluid.analyzers.MatrixMethodAnalyzer;
import jline.solvers.fluid.analyzers.MinNormalAnalyzer;
import jline.solvers.fluid.analyzers.RMFAnalyzer;
import jline.solvers.fluid.analyzers.TbiAnalyzer;
import jline.solvers.fluid.moments.FluidDaeApplicable;
import jline.solvers.fluid.moments.FluidMinNormalApplicable;
import jline.solvers.fluid.moments.FluidNonHyperbolicException;
import jline.solvers.fluid.moments.MvnRectangle;
import jline.solvers.fluid.handlers.MethodStepHandler;
import jline.solvers.fluid.handlers.PassageTimeODE;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.util.FastMath;
import odesolver.LSODA;

import java.util.*;

import static java.lang.Double.*;
import static jline.api.npfqn.Npfqn_nonexp_approx.npfqn_nonexp_approx;
import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.api.sn.SnOpenProbTerms.snOpenProbTerms;
import static jline.io.InputOutput.*;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.min;
import jline.api.sym.SageRestEngine;
import jline.api.sym.SymEngine;
import jline.api.sym.SymEngines;

/**
 * FLD - Fluid/Mean-Field Approximation solver.
 * SolverFluid is based on fluid and mean-field approximation methods.
 * It provides fluid-based analysis for queueing networks, particularly
 * useful for large-scale systems where discrete-event simulation would
 * be computationally expensive.
 */
public class SolverFluid extends NetworkSolver {

    /**
     * Creates a new SolverFluid instance with default options.
     *
     * @param model The network model to analyze
     */
    public SolverFluid(Network model) {
        this(model, SolverFluid.defaultOptions());
        this.result = new FluidResult();
    }

    /**
     * Creates a new SolverFluid instance with variable arguments for options.
     *
     * @param model The network model to analyze
     * @param varargin Variable arguments for solver options
     */
    public SolverFluid(Network model, Object... varargin) {
        this(model, SolverFluid.defaultOptions());
        this.result = new FluidResult();
        this.options = Solver.parseOptions(this.options, varargin);
    }

    /**
     * Creates a new SolverFluid that warm-starts the ODE integration from the
     * steady-state solution of an auxiliary solver (see
     * {@link jline.solvers.NetworkSolver#initFromSolver}).
     *
     * @param model The network model to analyze
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @param varargin Variable arguments for solver options
     */
    public SolverFluid(Network model, NetworkSolver initSolver, Object... varargin) {
        this(model, varargin);
        this.initFromSolver(initSolver);
    }

    /**
     * Creates a new SolverFluid instance with a specific method.
     *
     * @param model The network model to analyze
     * @param method The fluid analysis method to use
     */
    public SolverFluid(Network model, String method) {
        super(model, "SolverFluid", SolverFluid.defaultOptions().method(method));
        this.result = new FluidResult();
    }

    /**
     * Creates a new SolverFluid instance with specific options.
     *
     * @param model The network model to analyze
     * @param options The solver options to use
     */
    public SolverFluid(Network model, SolverOptions options) {
        super(model, "SolverFluid", options);
        this.result = new FluidResult();
    }

    /**
     * Returns the default solver options for the Fluid solver.
     *
     * @return Default solver options with SolverType.FLUID
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.FLUID);
    }

    /**
     * Returns the feature set supported by the Fluid solver
     *
     * @return - the feature set supported by the Fluid solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "Cache", "CacheClassSwitcher",  // CacheRetrieval deliberately NOT declared: no fluid code implements delayed-hit retrieval
                "Cox2", "Coxian", "Erlang", "Exp", "HyperExp",
                // MAP and MMPP2 are accepted at their stationary rate: the fluid
                // ODE has no representation of the arrival phase process, so the
                // autocorrelation is lost while flow stays conserved, exactly as
                // in MATLAB SolverFLD.
                "APH", "Det", "MAP", "MMPP2",
                // Non-Markovian renewal distributions: converted to acyclic PH via
                // snNonmarkovToPh in runAnalyzer, so the fluid ODE can solve them.
                "Gamma", "Lognormal", "Pareto", "Uniform", "Weibull",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
                "Server", "ServiceTunnel",
                "LoadDependence", // closing family only, see getMethodFeatureSet
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS",
                "SchedStrategy_GPS", // minnormal only, see getMethodFeatureSet
                "SchedStrategy_SIRO", "SchedStrategy_LCFS", "SchedStrategy_LCFSPR",
                // Native fluid cache models: RANDOM(m)/FIFO(m) (refined mean field)
                // and strict FIFO(m) (position-resolved mean field). LRU/HLRU/CLIMB/
                // QLRU have no drift-based fluid model and are rejected at runtime.
                // Petri-net constructs, for the 'dae' Petri arm (PetriSolver). A
                // stochastic Petri net has no drift outside the DAE form: its
                // conserved quantities are P-invariants rather than chain
                // populations, and an immediate transition is an algebraic FLOW
                // rather than an event with a rate. getMethodFeatureSet takes them
                // back off every other method, which would integrate the net as an
                // empty model and report zeros without a warning.
                "Place", "Transition", "Enabling", "Inhibiting", "Timing", "Firing",
                "Storage", "Linkage",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                "ReplacementStrategy_SFIFO",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ClosedClass", "SelfLoopingClass", "Replayer",
                // Non-homogeneous Poisson arrivals: steady state uses the
                // nominal (time-average) rate, getTranAvg tracks lambda(t)
                // through the per-event rate multiplier of the closing ODE.
                "NHPP",
                // Correlated Markovian arrivals: the closing and matrix ODEs carry the
                // full (D0, D1) phase structure, so a MAP or MMPP2 keeps its exact
                // stationary rate. Omitting them here rejected models the fluid solver
                // could already solve, and matched neither MATLAB nor python.
                "MAP", "MMPP2",
                "MAPt",
                "PHt",
                // Fork-join through the MMT transformation, driven by
                // jline.solvers.fj.FJFixedPoint with fldDispatch as the inner
                // solve (as in SolverMVA and SolverNC). The transform emits only
                // Source, Delay, Queue, Router and ClassSwitch, all of which the
                // fluid drift already carries.
                "Fork", "Forker", "Join", "Joiner",
                // quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                "JoinPartial",
                "RandomSource", "Sink", "Source", "OpenClass", "JobSink",
                // c-server stations: the drifts carry min(n,c); withdrawn from
                // "diffusion" and "mfq" in getMethodFeatureSet.
                "MultiServer",
                // A binding buffer: "dae" carries it as an algebraic constraint,
                // "mol" IS the Mt/G/s/0 loss system and the AoI arm of "mfq" is a
                // bufferless or single-buffer queue. WHICH method serves one is the
                // structural rule supportsModelMethod asks and runAnalyzer stops
                // on, so no per-method delta duplicates it here.
                "FiniteCapacity"
        });
        return featSupported;
    }

    /**
     * Per-method feature envelope, mirroring the MATLAB
     * {@code SolverFLD.getMethodFeatureSet}.
     *
     * <p>Defining this is what lets the solver gate name the offending features:
     * with no method feature set the check falls back to the coarse
     * {@code supports(model)}, which returns an empty reason, so a rejection
     * could only say "features not supported" without saying which ones.</p>
     *
     * @param method the concrete method name
     * @return the features that method accepts
     */
    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        FeatureSet featSupported = SolverFluid.getFeatureSet();
        // EVERY TEST BELOW IS ON THE CANONICAL NAME. Enumerating each spelling by
        // hand is what let the four codebases drift apart over an alias: the
        // Reneging branch below listed the qualified spellings but not the bare
        // "ggisgi"/"tga" that QsysLimitAnalyzer.handles accepts, and HOL was not
        // gated here at all. Canonicalize once and an alias cannot carry a
        // different envelope than the name it resolves to.
        String m = canonicalMethod(method);
        // Limited load dependence composes with the closure as a rate multiplier
        // alpha(n_i) on the scheduling share, which only the closing family evaluates.
        // The matrix, softmin, statedep, tbi, mfq and rmf paths build their drift
        // independently and would silently ignore alpha, so they must keep rejecting it.
        if (!isClosingFamily(m)) {
            featSupported.setFalse(new String[]{"LoadDependence"});
        }
        // Scheduling disciplines with no branch in the drift. A station without a case
        // in FluidRateFactors keeps rates = x, i.e. it is integrated as an INFINITE
        // SERVER, and the answer is wrong without any warning: on Delay(Z=1) ->
        // Queue(c=1), N=4, exact Q2 = 3.0154, the fall-through returns 2.0000. The
        // closing metric reader accepts SIRO as FCFS, so the ODE integrated it as INF
        // while the metrics were read as if it shared the server. Reject at the featset
        // gate instead, where the message names the offending discipline. matrix builds
        // a PS drift for every queueing station, which is the right aggregate for any
        // work-conserving discipline, so it is unaffected.
        if (isClosingFamily(m) || "statedep".equals(m) || "softmin".equals(m)
                || "tbi".equals(m)) {
            featSupported.setFalse(new String[]{
                    "SchedStrategy_SIRO", "SchedStrategy_LCFS", "SchedStrategy_LCFSPR"});
        }
        // GPS divides the server by weight among the BACKLOGGED classes, so its share is
        // a function of the backlog INDICATOR. A first-order closure cannot express it at
        // all: with continuous x_k > 0 every class is always backlogged and the share
        // collapses to the constant w_k/sum_j w_j, the heavy-traffic limit, regardless of
        // load. Only minnormal supplies the P(X_k >= 1) that the closure needs. matrix
        // builds a PS drift, and GPS is not PS: equal-weight GPS gives each backlogged
        // CLASS an equal share, PS each JOB.
        if (!"minnormal".equals(m)) {
            featSupported.setFalse(new String[]{"SchedStrategy_GPS"});
        }
        // HOL allocates capacity in PRIORITY order, not in proportion to population,
        // and no fluid drift reads sn.classprio except the single-queue MFQ priority
        // branch. This port declared HOL for EVERY fluid method, so a priority model
        // was offered to drifts that answer it as if the classes shared the server
        // proportionally. MATLAB SolverFLD and the C++ twin both gate it here.
        if (!"mfq".equals(m)) {
            featSupported.setFalse(new String[]{"SchedStrategy_HOL"});
        }
        // MULTISERVER (registry name since 2026-09-05): the drifts carry min(n,c)
        // except the diffusion SDE, which is written for one or infinitely many
        // servers, and MFQ, a single-queue model on the same server counts (off
        // them "mfq" resolves to "matrix", so the delta binds only where it runs
        // as itself). The structural refusal keeps wording the diffusion case.
        if ("diffusion".equals(m) || "mfq".equals(m)) {
            featSupported.setFalse(new String[]{"MultiServer"});
        }
        if ("ggisgi.fluid".equals(m) || "ggingi.tga".equals(m) || "tvms".equals(m)) {
            // The only fluid methods in LINE stated for a queue customers
            // ABANDON. Reneging stays out of the base FLD envelope: the network
            // drift carries no abandonment flow, so every other method would
            // integrate the model as if nobody left.
            featSupported.setTrue(new String[]{"Reneging"});
        }
        if (QsysLimitAnalyzer.handles(m)) {
            // Every one of them is stated for a single open station; the base
            // envelope's closed classes have no meaning there.
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        }
        if ("refined".equals(m)) {
            // CLOSED MODELS ONLY, which the MATLAB runAnalyzer has always
            // enforced by name and the featset never stated: the 1/N correction
            // is solved on orth(D) over the FULL state, so on an open model it
            // adds a perturbation to the SOURCE POOL mass, a normalisation
            // constant rather than a population. Only "minnormal" was validated
            // open. Stating it here is what lets a report withdraw the pair
            // instead of offering a run that stops -- on an open fork-join model
            // the same restriction surfaced as a failure inside the MMT fixed
            // point rather than as a refusal.
            featSupported.setFalse(new String[]{"OpenClass", "Source", "Sink",
                    "RandomSource", "JobSink"});
        }
        if ("diffusion".equals(m) || "kp".equals(m)) {
            // NEITHER OF THESE TWO INTEGRATES A FORK-JOIN MODEL, and each says so
            // by answering rather than by refusing, which is the reason to state
            // it here. Measured on a SYMMETRIC closed fork-join (Delay -> Fork ->
            // two identical FCFS queues -> Join, N = 2) whose exact chain is
            // Q1 = Q2 = 0.664, J = 0.624, D = 1.024: "diffusion" returns the
            // whole population on ONE station and zero elsewhere -- a different
            // station on a rerun, so the SDE is not integrating this model at
            // all -- and "kp" returns an ALL-ZERO table on a symmetric OPEN
            // fork-join fed at rate 0.5, an empty network where jobs are
            // arriving. The C++ featset has always withheld the names.
            featSupported.setFalse(new String[]{"Fork", "Join", "Forker", "Joiner",
                    "JoinPartial"});
        }
        if ("diffusion".equals(m)) {
            // The diffusion SDE PROJECTS each class back onto its own fixed
            // population at every step, which is the closed-network constraint
            // itself: an open class has no population to project onto, and a
            // Source is not a station the SDE has a coordinate for.
            featSupported.setFalse(new String[]{"OpenClass", "Source", "Sink",
                    "RandomSource", "JobSink"});
        }
        if ("tbi".equals(m)) {
            // Trajectory-based iteration decomposes the CLOSED population into
            // cells and relaxes the waveforms between them; there is no cell for
            // an unbounded open stream. A cache model is solved by decomposition
            // rather than by one drift, so the cell partition has nothing to
            // partition -- use "rmf".
            featSupported.setFalse(new String[]{"OpenClass", "Source", "Sink",
                    "RandomSource", "JobSink", "Cache", "CacheClassSwitcher",
                    "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                    "ReplacementStrategy_SFIFO"});
        }
        if ("kp".equals(m)) {
            // The Ko-Pender limits are proved for an OPEN network of stations fed
            // by external arrival processes: a closed class has no arrival process
            // to modulate and no source phase to carry, and the cache and
            // class-switch machinery has no counterpart in the paper's event set.
            // Narrow the envelope rather than return the drift of a model the
            // limits do not describe. Mirrors native python and the MATLAB
            // reference.
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass",
                    "Cache", "CacheClassSwitcher", "ClassSwitch", "StatelessClassSwitcher",
                    "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                    "ReplacementStrategy_SFIFO"});
        }
        // A stochastic Petri net has no drift outside the DAE form: its conserved
        // quantities are P-invariants rather than chain populations, an immediate
        // transition is an algebraic FLOW rather than an event with a rate, and a
        // bounded place is a linear inequality on the marking. Every other fluid
        // method builds its drift from the station/class/phase encoding, where a
        // Place contributes no coordinate at all, so it would integrate the net as
        // an empty model and report zeros without a warning.
        if (!"dae".equals(m)) {
            featSupported.setFalse(new String[]{"Place", "Transition", "Enabling", "Inhibiting",
                    "Timing", "Firing", "Storage", "Linkage"});
        }
        if ("dae".equals(m)) {
            // DPS closes on the covariance BETWEEN the class coordinates of a
            // station, not on the station total. 'minnormal' carries those
            // blocks through its outer iteration; the DAE form has no unknown for
            // them, since a matrix block per station restores the quartic cost
            // that keeping Sigma out of the Newton vector avoids. GPS is already
            // excluded above, for 'minnormal' only.
            featSupported.setFalse(new String[]{"SchedStrategy_DPS"});
            // A finite capacity region is a linear inequality on the state, which
            // the DAE form can carry as an algebraic equation beside the drift
            // and no ODE method can carry at all. The gate is where this has to
            // be declared: the static getFeatureSet must keep Region false
            // because every other method has to go on rejecting one.
            // DaeAnalyzer still refuses, by name, the region forms that are not
            // constraints on this drift.
            featSupported.setTrue(new String[]{"Region"});
        }
        return featSupported;
    }

    /**
     * The structural finite-capacity gate {@code runAnalyzer} enforces at solve
     * time, stated here so that a CALLER can see it before running.
     *
     * <p>Nothing in the fluid tree reads sn.cap or sn.classcap, so every method
     * but two integrates a capped station as an unbounded one. "dae" carries
     * the buffer as an algebraic constraint on the drift, and "mol" is stated
     * for the Mt/G/s/0 LOSS system, where the server count IS the buffer. There is no
     * registry feature name for plain capacity, hence the structural test;
     * SolverNC and SolverMVA gate the same way. The exemption list is the one
     * runAnalyzer applies, so the two cannot disagree.
     *
     * <p>Left only in runAnalyzer the rule was invisible to every gate above it,
     * and SolverAUTO.listValidMethods offered every fluid method on the
     * BAS-blocking model of cqn_bas_blocking, each of which then threw when
     * asked to run.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        String reason = super.supportsModelMethod(method);
        if (reason != null && !reason.isEmpty()) {
            return reason;
        }
        // The time-varying single-station limits report a TRAJECTORY, so they
        // need a finite options.timespan. A horizon is an option and not a model
        // feature, hence the structural test; the predicate is the one
        // QsysLimitAnalyzer stops on, so the report and the run cannot answer
        // differently.
        if (isTimeVaryingLimit(method)) {
            String horizonReason = QsysLimitAnalyzer.horizonReason(this.options);
            if (!horizonReason.isEmpty()) {
                return "The '" + method + "' method reports a trajectory. " + horizonReason;
            }
        }
        // A fork-join model is answered by the MMT fixed point rather than by
        // one drift, and not every method can run it. Fork and OpenClass are
        // both declared names, so the featset cannot state a rule that is their
        // CONJUNCTION; it is structural, and it is the predicate runAnalyzer
        // stops on.
        if (this.model != null) {
            String fjReason = forkJoinAdmitsReason(this.model.getStruct(false), method);
            if (!fjReason.isEmpty()) {
                return fjReason;
            }
        }
        if ("dae".equals(method) || "fluid.dae".equals(method)
                || "mol".equals(method) || "fluid.mol".equals(method)) {
            return "";
        }
        if (this.model == null) {
            return "";
        }
        // 'default' IS ASKED THROUGH ITS RESOLUTION, not as a name of its own.
        // On a capped model runAnalyzer now resolves it to "dae" (see
        // blockedResolvesToDae), so gating the literal name against the
        // capacity rule would refuse the very run that goes on to succeed.
        if (("default".equals(method) || "fluid.default".equals(method))
                && blockedResolvesToDae(this.model, this.model.getStruct(false), this.options)) {
            return "";
        }
        String capReason = NetworkSolver.bindingCapacityReason(this.model,
                this.model.getStruct(false), "SolverFluid");
        if (capReason != null) {
            return capReason + " Use options.method = \"dae\", which carries the buffer as an "
                    + "algebraic constraint on the drift.";
        }
        return "";
    }

    /**
     * Does {@code options.method = "default"} stand for {@code "dae"} on this
     * model? True exactly when a buffer or a capacity region BINDS and the DAE
     * route accepts the model.
     *
     * <p>A BINDING BUFFER OR A REGION HAS ONE FLUID ROUTE, for the same reason a
     * Petri net does: nothing else in the fluid tree reads sn.cap, sn.classcap or
     * the region limit, so every other method integrates the capped station as an
     * unbounded one -- which is why runAnalyzer refuses them. Resolving "default"
     * to one of those turned a model this solver CAN answer into an error whose
     * advice was to type the very method the resolution should have picked.
     *
     * <p>The capacity test is the gate's own, so the two cannot disagree. Where
     * the DAE route declines, this returns false and the gate speaks, naming the
     * blocking feature. Mirrors the MATLAB fluid_resolve_default_method.
     *
     * @param model   the network
     * @param sn      its struct
     * @param options the solver options, for the dae applicability limits
     * @return true when "default" must resolve to "dae"
     */
    public static boolean blockedResolvesToDae(Network model, NetworkStruct sn, SolverOptions options) {
        if (model == null || sn == null) {
            return false;
        }
        boolean blocked = sn.nregions > 0
                || NetworkSolver.bindingCapacityReason(model, sn, "SolverFluid") != null;
        return blocked && FluidDaeApplicable.reasonToDecline(sn, options) == null;
    }

    /** Methods that evaluate the closing rate factors, and so the load-dependent scaling. */
    /**
     * The method name without its {@code fluid.} prefix, which names the same
     * method, with {@code butools} and {@code aoi} folded onto {@code mfq}: they
     * name the MFQ branch's backend and its age-of-information reading rather
     * than methods of their own, which is how native python spells them and how
     * getAvgAoI reaches it (it requires method='mfq').
     *
     * @param method the requested method name
     * @return the canonical spelling
     */
    public static String fluidUnqualify(String method) {
        if (method == null) {
            return null;
        }
        String m = method;
        if (m.length() > 6 && m.startsWith("fluid.")) {
            m = m.substring(6);
        }
        if ("butools".equals(m) || "aoi".equals(m)) {
            return "mfq";
        }
        return m;
    }

    /**
     * Can {@code method} run the fluid fork-join fixed point on this model?
     *
     * <p>A fork-join model is not integrated as one drift: the MMT transform
     * replaces the fork by auxiliary classes and the answer is the fixed point
     * of solving that transformed model repeatedly. On a CLOSED model the
     * transform stays closed and every fluid method takes it. On an OPEN one the
     * auxiliary classes arrive at a Source, and the DAE form has no unknowns for
     * them: the inner solve fails on the class count rather than returning a
     * drift, so the method is refused by name instead.</p>
     *
     * <p>"refined" is NOT listed here even though it fails the same way, because
     * it is already refused on every open model, fork-join or not, by its own
     * closed-model restriction (see {@link #getMethodFeatureSet}).</p>
     *
     * <p>Called by runAnalyzer, so the run stops on it, and by
     * {@link #supportsModelMethod}, so a caller sees the same verdict before
     * paying for the fixed point. One predicate, two callers.</p>
     *
     * @param sn     the network structure
     * @param method the concrete method name
     * @return empty string when the method may run the fixed point here, else the refusal
     */
    public static String forkJoinAdmitsReason(NetworkStruct sn, String method) {
        if (!("dae".equals(method) || "fluid.dae".equals(method))) {
            return "";
        }
        boolean anyFork = false;
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Fork) {
                anyFork = true;
                break;
            }
        }
        if (!anyFork) {
            return "";
        }
        boolean anyOpen = false;
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                anyOpen = true;
                break;
            }
        }
        if (!anyOpen) {
            return "";
        }
        return "The dae method has no route through the fork-join fixed point on an OPEN "
                + "model: the MMT transform hands the inner solve a mixed network whose "
                + "auxiliary open classes the DAE form carries no unknowns for. Use "
                + "options.method = \"minnormal\", which is the same closure and does run "
                + "that fixed point.";
    }

    /**
     * The three single-station limits that report a trajectory rather than a
     * stationary point, and so need a finite options.timespan; "ggisgi" and
     * "tga" are stationary and are not among them.
     *
     * @param method the concrete method name
     * @return true when the method integrates over a finite horizon
     */
    private static boolean isTimeVaryingLimit(String method) {
        return "tvms".equals(method) || "fluid.tvms".equals(method)
                || "mtginf".equals(method) || "fluid.mtginf".equals(method)
                || "mol".equals(method) || "fluid.mol".equals(method);
    }

    /**
     * The one spelling of a fluid method that every gate tests against.
     *
     * <p>Three families of alias reach this solver and they used to be expanded by
     * hand at each branch, which is why the four codebases drifted apart: a
     * {@code fluid.} qualifier the dispatch accepts on every name, the MFQ backend
     * aliases {@code butools} and {@code aoi}, and the short spellings
     * {@code ggisgi} and {@code tga} of the two single-station limits.
     * Canonicalizing once is what makes an alias carry the same feature envelope as
     * the name it resolves to; MATLAB {@code SolverFLD.canonicalMethod}, native
     * python and C++ apply the same three rules in the same order.</p>
     *
     * @param method the method name as the caller spelled it
     * @return the canonical spelling
     */
    public static String canonicalMethod(String method) {
        if (method == null) {
            return null;
        }
        String m = method.startsWith("fluid.") ? method.substring(6) : method;
        // The MFQ backend and its age-of-information reading are the same drift as
        // 'mfq': the dispatch sends all three to the MFQ analyzer.
        if ("butools".equals(m) || "aoi".equals(m)) {
            return "mfq";
        }
        // The short spellings of the two single-station limits, as
        // QsysLimitAnalyzer.canonical and the C++ fluid_qsys_canonical map them.
        if ("ggisgi".equals(m)) {
            return "ggisgi.fluid";
        }
        if ("tga".equals(m)) {
            return "ggingi.tga";
        }
        return m;
    }

    private static boolean isClosingFamily(String method) {
        String m = canonicalMethod(method);
        return "closing".equals(m)
                || "minnormal".equals(m)
                // 'refined' IS the min-normal closure with the O(1/N) correction on
                // top: same drift, same rate factors, so it evaluates the same
                // alpha(n_i) multiplier and needs the same SIRO/LCFS/LCFSPR strip.
                // Omitting it here refused a load-dependent model on 'refined' that
                // MATLAB and C++ both answer, and left 'refined' declaring three
                // disciplines its drift silently integrates as INF.
                || "refined".equals(m)
                // 'dae' IS the min-normal closure -- same drift, same rate
                // factors -- solved as one system instead of by substitution, so
                // it evaluates the same alpha(n_i) multiplier.
                || "dae".equals(m);
    }

    /**
     * Computes transient average station metrics, tracking any non-homogeneous
     * (NHPP) arrival intensity over the timespan.
     *
     * <p>The NHPP schedules of the Source stations are injected into
     * {@code options.config.nhpp_sched}, so that the closing ODE follows
     * lambda(t) rather than the baked-in time-average rate. Steady-state
     * {@code getAvg} is unaffected (no schedule is injected there), consistent
     * with defining the NHPP steady state as its time average; the injected
     * schedule is removed again on return. Mirrors MATLAB
     * {@code SolverFLD/getTranAvg.m}.</p>
     */
    @Override
    public void getTranAvg() {
        java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.NhppEntry> previous =
                (this.options.config != null) ? this.options.config.nhpp_sched : null;
        java.util.List<jline.solvers.fluid.handlers.FluidRateMultiplier.NhppEntry> sched =
                FluidNhpp.detectNhppSchedule(model.getStruct(true));
        String previousMethod = this.options.method;
        if (this.options.config != null && !sched.isEmpty()) {
            this.options.config.nhpp_sched = sched;
        }
        // THE TRANSIENT MEANS ARE READ OFF AN INTEGRATED TRAJECTORY, so the
        // method has to be one that produces one. Switching only on an NHPP
        // schedule left "default" to resolve to "minnormal", whose moment
        // closure returns a CONVERGED FIXED POINT rather than a trajectory: on
        // renv_node_breakdown's UP stage (lambda 0.8, mu 2.0, fluid mean
        // lambda/mu = 0.4) this engine reported 0.5095 where MATLAB, native
        // Python and C++ all report 0.4000, and SolverENV -- which couples its
        // stages THROUGH these trajectories -- carried that into a 6% error on
        // QLen and an 11% one on throughput, breaking flow balance against an
        // arrival rate it must reproduce. Mirrors the switch in MATLAB
        // SolverFLD/getTranAvg.m:
        //   default/matrix/closing  -> closing;
        //   tbi                     -> kept, it integrates the closing rates by
        //                              cell decomposition and is a trajectory;
        //   minnormal/refined       -> kept, they integrate the closing ODEs
        //                              with the converged variance held fixed,
        //                              unless an NHPP makes the drift
        //                              non-autonomous, which has no stationary
        //                              covariance to hold fixed.
        // `rmf` is the one arm NOT mirrored: MATLAB switches it and recovers the
        // cache transient through solver_fld_cacheqn_tran, which this port does
        // not carry, so switching here would answer a cache model's transient
        // with the non-cache closing ODE.
        // `kp` is kept for the opposite reason: KoPenderAnalyzer INTEGRATES the
        // fluid and diffusion limits over the horizon, so it already returns a
        // trajectory, and it is the only method that also carries the second
        // moment (QVart/Sigmat) getTranAvgVar reads. Switching it away discards
        // that covariance and hands a MAPt/PHt source to the closing ODE, which
        // builds its per-event multiplier from the time-averaged nominal and
        // raises on the segment matrices (1xh against hxh).
        String m = (previousMethod == null) ? "default" : previousMethod;
        boolean isTbi = "tbi".equals(m) || "fluid.tbi".equals(m);
        boolean isRmf = "rmf".equals(m) || "fluid.rmf".equals(m);
        boolean isKp = "kp".equals(m) || "fluid.kp".equals(m);
        boolean isClosure = "minnormal".equals(m) || "fluid.minnormal".equals(m)
                || "refined".equals(m) || "fluid.refined".equals(m)
                // 'dae' integrates the closure itself over the horizon, with
                // conservation as an algebraic equation and the covariance
                // advancing alongside the mean, so its trajectory is the
                // transient of the closed system rather than a first-order
                // stand-in. Same non-autonomous exclusion as the other closures:
                // a time-varying drift has no stationary covariance for the seed
                // solve to converge to.
                || "dae".equals(m) || "fluid.dae".equals(m);
        // The single-station fluid limits produce their own trajectory over the
        // horizon -- that is what they are -- so switching them to the closing
        // ODE would answer a different model with a time-averaged rate.
        boolean isQsysLimit = QsysLimitAnalyzer.handles(m);
        boolean toClosing = !isTbi && !isRmf && !isKp && !isQsysLimit
                && !(isClosure && sched.isEmpty());
        if (toClosing && !"closing".equals(m) && !"fluid.closing".equals(m)) {
            this.options.method = "closing";
            this.reset();
        }
        try {
            super.getTranAvg();
        } finally {
            if (this.options.config != null) {
                this.options.config.nhpp_sched = previous;
            }
            this.options.method = previousMethod;
        }
    }

    public DistributionResult getCdfRespT() {
        long startTime = System.nanoTime();
        this.getAvg(); // Get steady-state solution
        
        // The number of phases may have changed during fluid iterations
        if (this.result instanceof FluidResult && ((FluidResult) this.result).snFinal != null) {
            this.sn = ((FluidResult) this.result).snFinal;
	    // Also update phasessz and phaseshift to be consistent with the new phases
            this.sn.phasessz = this.sn.phases.copy();
            for (int row = 0; row < this.sn.phasessz.getNumRows(); row++) {
                for (int col = 0; col < this.sn.phasessz.getNumCols(); col++) {
                    if (this.sn.phasessz.get(row, col) < 1) {
                        this.sn.phasessz.set(row, col, 1);
                    }
                }
            }
            
            this.sn.phaseshift = new Matrix(this.sn.phases.getNumRows(), 1);
            this.sn.phaseshift = Matrix.concatColumns(this.sn.phaseshift, this.sn.phasessz.cumsumViaRow(), null);
        } else {
            this.sn = model.getStruct(true);
        }
        
        this.options.init_sol = ((FluidResult) this.result).odeStateVec;
        Matrix[][] passageTimeResults = passageTime();
        ((FluidResult) this.result).distribC = passageTimeResults;
        ((FluidResult) this.result).distribRuntime = (System.nanoTime() - startTime) / 1000000000.0;

        // Create and populate the distribution result
        DistributionResult distResult = new DistributionResult(sn.nstations, sn.nclasses, "response_time");
        
        // Populate CDF data from passage time results
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (passageTimeResults[i][k] != null && !passageTimeResults[i][k].isEmpty()) {
                    distResult.setCdf(i, k, passageTimeResults[i][k]);
                }
            }
        }
        
        distResult.runtime = ((FluidResult) this.result).distribRuntime;
        return distResult;
    }

    /**
     * Get cumulative distribution function of response times, handle form.
     *
     * <p>The handle only names which metrics the caller wants; the passage-time
     * solve computes every (station, class) anyway, exactly as the reference
     * {@code @SolverFLD/getCdfRespT.m}, so it is accepted for signature
     * compatibility and not read. Without this override a caller holding a
     * {@code NetworkSolver} reference would silently fall through to the base
     * class exponential fit.</p>
     *
     * @param R the response time handles, accepted for signature compatibility
     * @return DistributionResult containing the response time CDF data
     */
    @Override
    public DistributionResult getCdfRespT(AvgHandle R) {
        return getCdfRespT();
    }

    /** Per-station population variance as a column Matrix, null-preserving. */
    private static Matrix sigma2Matrix(double[] sigma2) {
        if (sigma2 == null) {
            return null;
        }
        Matrix out = new Matrix(sigma2.length, 1, sigma2.length);
        for (int i = 0; i < sigma2.length; i++) {
            out.set(i, 0, sigma2[i]);
        }
        return out;
    }

    /** Whether the model owns a cache node, i.e. needs the decomposition route. */
    private static boolean hasCacheNodes(NetworkStruct sn) {
        if (sn == null || sn.nodetype == null) {
            return false;
        }
        for (int ind = 0; ind < sn.nodetype.size(); ind++) {
            if (sn.nodetype.get(ind) == jline.lang.constant.NodeType.Cache) {
                return true;
            }
        }
        return false;
    }

    /**
     * Per-class job counts at a station: the caller's {@code state_a} when given,
     * else the marginal of the model's own state row.
     */
    private Matrix perClassCounts(int ist, Matrix state_a) {
        if (state_a != null && state_a.length() > 0) {
            Matrix out = new Matrix(1, sn.nclasses, sn.nclasses);
            for (int r = 0; r < sn.nclasses && r < state_a.length(); r++) {
                out.set(0, r, state_a.get(r));
            }
            return out;
        }
        State.StateMarginalStatistics stats =
                ToMarginal.toMarginal(
                        this.sn,
                        ist,
                        sn.state.get(this.model.getStations().get((int) sn.stationToStateful.get(0, ist))),
                        null,
                        null,
                        null,
                        null,
                        null);
        return stats.nir;
    }

    public ProbabilityResult getProbAggr(int ist) {
        return getProbAggr(ist, null);
    }

    /**
     * Probability of a given per-class job distribution at a station. A null
     * {@code state_a} reads the model's own state, as the one-argument form does.
     *
     * <p>The explicit form exists for the delegating callers: the model
     * interchange carries no initial state, so a caller that set one with
     * {@code initFromMarginal} has to name the cell it is asking about, or every
     * query would be answered at the default initialization.</p>
     *
     * @param ist     station index
     * @param state_a per-class job counts, or null to read the model state
     * @return scalar probability in [0,1]
     */
    @Override
    public ProbabilityResult getProbAggr(int ist, Matrix state_a) {

        if (ist > sn.nstations) {
            throw new RuntimeException("Station number exceeds the number of stations in the model.");
        }

        if (!this.hasAvgResults()) {
            this.getAvg();
        }

        // Re-read the struct AFTER the analysis: the one captured at construction
        // predates the model's default initialization, so its state rows are not
        // the ones every metric above was computed for.
        this.sn = this.model.getStruct(false);

        // The moment closure supplies the JOINT law of the per-class populations,
        // so the answer is the probability its multivariate normal assigns to the
        // unit cell around the current state. A Source is excluded because its
        // coordinate is a normalisation constant rather than a population and
        // carries no covariance (see FluidMomentTerms).
        FluidResult fres = (FluidResult) this.result;
        if (fres.momentSigma != null && fres.momentClassBlock != null
                && sn.sched.get(this.model.getStations().get(ist)) != SchedStrategy.EXT
                && !hasOpenClassAt(ist, fres.momentClassBlock)) {
            return gaussianCellProb(ist, state_a);
        }

        boolean allAreFinite = true;
        for (int i = 0; i < sn.njobs.getNumCols(); i++) {
            if (isInfinite(sn.njobs.get(0, i))) {
                allAreFinite = false;
                break;
            }
        }

        if (allAreFinite) {
            // Binomial approximation with mean fitted to queue-lengths.
            // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
            Matrix N = sn.njobs;
            Matrix Q = this.result.QN;
            Matrix nir = perClassCounts(ist, state_a);
            ((FluidResult) this.result).logPnir = 0;
            for (int r = 0; r < nir.getNumCols(); r++) {
                int Nr = (int) N.get(0, r);
                int nirVal = (int) nir.get(0, r);
                double Qir = Q.get(ist, r);
                ((FluidResult) this.result).logPnir += Maths.logBinomial(Nr, nirVal);
                ((FluidResult) this.result).logPnir += nirVal * FastMath.log(Qir / Nr);
                ((FluidResult) this.result).logPnir += (Nr - nirVal) * FastMath.log(1 - Qir / Nr);
            }
            ((FluidResult) this.result).Pnir = FastMath.exp(((FluidResult) this.result).logPnir);
            return new ProbabilityResult(((FluidResult) this.result).Pnir);
        } else {
            // Mixed or open model: the open classes' product form at this
            // station, times the closed classes' binomials.
            Matrix N = sn.njobs;
            Matrix Q = this.result.QN;
            Matrix nir = perClassCounts(ist, state_a);
            double logPnir = snOpenProbTerms(this.sn, Q, this.result.UN, nir, ist);
            for (int r = 0; r < nir.getNumCols(); r++) {
                if (isInfinite(N.get(0, r))) {
                    continue;
                }
                int Nr = (int) N.get(0, r);
                int nirVal = (int) nir.get(0, r);
                double Qir = Q.get(ist, r);
                logPnir += Maths.logBinomial(Nr, nirVal);
                logPnir += nirVal * FastMath.log(Qir / Nr);
                logPnir += (Nr - nirVal) * FastMath.log(1 - Qir / Nr);
            }
            ((FluidResult) this.result).logPnir = logPnir;
            ((FluidResult) this.result).Pnir = FastMath.exp(logPnir);
            return new ProbabilityResult(((FluidResult) this.result).Pnir);
        }
    }

    /**
     * Joint probability of the per-class populations at a station under the linear
     * noise approximation solved by the moment closure. Java twin of the local
     * {@code local_gaussian_cell} of the MATLAB {@code @SolverFLD/getProbAggr}.
     *
     * <p>The state coordinates of class r at the station are
     * {@code momentClassBlock[ist*K+r]} (one per service phase), so the class
     * population is their sum: its mean is the reported QN and the class-to-class
     * covariance is the sum of the corresponding block of momentSigma. The integer
     * count n is then read off the continuous law as the unit cell [n-1/2, n+1/2],
     * with the two ends extended to infinity at the boundaries of the state space,
     * so that the mass the normal puts on negative populations lands on the empty
     * station and the mass above a closed population lands on the full one.</p>
     */
    /**
     * Whether an OPEN class is served at a station.
     *
     * <p>The Gaussian cell is used only where it beats the alternative. For an
     * open class the first-order path is not an independence heuristic but the
     * exact product form of the underlying queue, geometric at a queue and
     * Poisson at a Delay, so replacing it by a normal approximation of the same
     * law would be a loss: on M/M/1 at rho = 0.5 the product form is exact where
     * the cell of the linear noise approximation returns 0.39 for the empty queue
     * against 0.50. The closure earns its place on the CLOSED populations, where
     * the alternative is Schmidt's binomial, itself an approximation, and where
     * correlation between the classes is real.</p>
     */
    private boolean hasOpenClassAt(int ist, int[][] classBlock) {
        int K = sn.nclasses;
        for (int r = 0; r < K; r++) {
            int[] blk = classBlock[ist * K + r];
            if (blk != null && blk.length > 0 && isInfinite(sn.njobs.get(0, r))) {
                return true;
            }
        }
        return false;
    }

    private ProbabilityResult gaussianCellProb(int ist, Matrix state_a) {
        FluidResult fres = (FluidResult) this.result;
        int K = sn.nclasses;
        Matrix nir = perClassCounts(ist, state_a);
        Matrix sigma = fres.momentSigma;
        int[][] classBlock = fres.momentClassBlock;

        List<Integer> idx = new ArrayList<Integer>();
        for (int r = 0; r < K; r++) {
            int[] blk = classBlock[ist * K + r];
            if (blk == null || blk.length == 0) {
                // the class has no service process here, so it has no coordinate:
                // any positive count is impossible rather than improbable
                if (nir.get(0, r) > 0) {
                    fres.logPnir = NEGATIVE_INFINITY;
                    fres.Pnir = 0;
                    return new ProbabilityResult(0.0);
                }
                continue;
            }
            idx.add(r);
        }

        int nr = idx.size();
        if (nr == 0) {
            fres.logPnir = 0;
            fres.Pnir = 1;
            return new ProbabilityResult(1.0);
        }

        double[] m = new double[nr];
        double[] a = new double[nr];
        double[] b = new double[nr];
        for (int u = 0; u < nr; u++) {
            int r = idx.get(u);
            m[u] = this.result.QN.get(ist, r);
            double n = nir.get(0, r);
            a[u] = (n <= 0) ? NEGATIVE_INFINITY : n - 0.5;
            double nr_jobs = sn.njobs.get(0, r);
            b[u] = (!isInfinite(nr_jobs) && n >= nr_jobs) ? POSITIVE_INFINITY : n + 0.5;
        }

        double[][] C = new double[nr][nr];
        for (int u = 0; u < nr; u++) {
            int[] bu = classBlock[ist * K + idx.get(u)];
            for (int v = u; v < nr; v++) {
                int[] bv = classBlock[ist * K + idx.get(v)];
                double acc = 0;
                for (int p = 0; p < bu.length; p++) {
                    for (int q = 0; q < bv.length; q++) {
                        acc += sigma.get(bu[p], bv[q]);
                    }
                }
                C[u][v] = acc;
                C[v][u] = acc;
            }
        }

        double p = MvnRectangle.probability(m, C, a, b);
        fres.Pnir = p;
        fres.logPnir = p > 0 ? FastMath.log(p) : NEGATIVE_INFINITY;
        return new ProbabilityResult(p);
    }

    public DistributionResult getTranCdfPassT() {
        long startTime = System.nanoTime();
        this.sn = model.getStruct(true);
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind, 0) == 1) {
                int isf = (int) sn.nodeToStateful.get(ind, 0);
                Matrix statePrior = sn.stateprior.get(this.model.getStatefulNodes().get(isf));
                if (statePrior.elementSum() != statePrior.elementMax()) {
                    throw new RuntimeException(
                            "getTranCdfPassT: multiple initial states have non-zero prior - unsupported.");
                }
                // Assign initial state to network
                sn.state.put(
                        this.model.getStatefulNodes().get(isf),
                        Matrix.extractRows(sn.state.get(this.model.getStatefulNodes().get(isf)), 0, 1, null));
            }
        }

        initSol();
        Matrix[][] passageTimeResults = passageTime();
        ((FluidResult) this.result).distribC = passageTimeResults;
        ((FluidResult) this.result).distribRuntime = (System.nanoTime() - startTime) / 1000000000.0;

        DistributionResult distResult = new DistributionResult(sn.nstations, sn.nclasses, "passage_time");
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (passageTimeResults[i][k] != null && !passageTimeResults[i][k].isEmpty()) {
                    distResult.setCdf(i, k, passageTimeResults[i][k]);
                }
            }
        }
        distResult.runtime = ((FluidResult) this.result).distribRuntime;
        return distResult;
    }

    /**
     * Initializes the solution vector for the fluid analysis.
     * This method sets up the initial state representation for the fluid solver.
     */
    public void initSol() {

        Matrix initSol = new Matrix(1, 0);
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind, 0) == 1) {
                int isf = (int) sn.nodeToStateful.get(ind);
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist < 0) {
                    // A stateful node that is not a station: a Router, which the
                    // MMT transformation puts where every fork was. It holds no
                    // jobs, so it contributes no coordinate to the initial
                    // condition. solver_fluid_initsol.m skips it the same way.
                    continue;
                }
                Matrix state_i = new Matrix(1, 0);
                // Compared to state_i, initSol_i does not track disabled classes
                // and removes Inf entries in the Sources
                Matrix initSol_i = new Matrix(1, 0);

                State.StateMarginalStatistics stats =
                        ToMarginal.toMarginal(
                                sn, ind, sn.state.get(this.model.getStatefulNodes().get(isf)), null, null, null, null, null);
                Matrix nir = stats.nir;
                List<Matrix> kir_i = stats.kir;
                
                // Declare variables outside switch to avoid scope issues
                int rMax, kMax;

                switch (sn.sched.get(this.model.getStations().get(ist))) {
                    case EXT:
                        state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
                        state_i.set(0, 0, POSITIVE_INFINITY); // Fluid does not model infinite buffer?
                        
                        // Guard against empty kir_i which occurs with unsupported node types (e.g., Fork/Join)
                        if (kir_i.isEmpty()) {
                            NodeType nodeType = sn.nodetype.get(ind);
                            throw new RuntimeException("SolverFluid does not support " + nodeType + " nodes. " +
                                    "Empty job phase data indicates incompatible node type. " +
                                    "Use SolverMVA or SolverJMT for fork-join models.");
                        }
                        
                        rMax = kir_i.get(0).getNumCols();
                        for (int r = 0; r < rMax; r++) {
                            kMax = sn.mu.get(this.model.getStations().get(ist)).get(sn.jobclasses.get(r)).length();
                            for (int k = 0; k < kMax; k++) {
                                state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
                                state_i.set(0, state_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                                if (!isNaN(sn.rates.get(ist, r))) {
                                    initSol_i.expandMatrix(
                                            1, initSol_i.getNumCols() + 1, initSol_i.getNumElements() + 1);
                                    initSol_i.set(0, initSol_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                                }
                            }
                        }
                        break;

                    case FCFS:
                    case SIRO:
                    case LCFS:
                    case LCFSPR:
                    case PS:
                    case INF:
                    case DPS:
                    case GPS:
                    case HOL:
                        // Guard against empty kir_i which occurs with unsupported node types (e.g., Fork/Join)
                        if (kir_i.isEmpty()) {
                            NodeType nodeType = sn.nodetype.get(ind);
                            throw new RuntimeException("SolverFluid does not support " + nodeType + " nodes. " +
                                    "Empty job phase data indicates incompatible node type. " +
                                    "Use SolverMVA or SolverJMT for fork-join models.");
                        }
                        rMax = kir_i.get(0).getNumCols();
                        for (int r = 0; r < rMax; r++) {
                            kMax = sn.mu.get(this.model.getStations().get(ist)).get(sn.jobclasses.get(r)).length();
                            for (int k = 0; k < kMax; k++) { // iterate over all phases
                                state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
                                if (!isNaN(sn.rates.get(ist, r))) {
                                    initSol_i.expandMatrix(
                                            1, initSol_i.getNumCols() + 1, initSol_i.getNumElements() + 1);
                                }
                                if (k == 0) {
                                    double sumKir_i = 0;
                                    for (int m = 1; m < kir_i.size(); m++) {
                                        sumKir_i += kir_i.get(m).get(0, r); // Accumulate the job counts across phases
                                    }
                                    // Jobs in waiting buffer are re-started phase 1
                                    state_i.set(0, state_i.getNumCols() - 1, nir.get(0, r) - sumKir_i);
                                    if (!isNaN(sn.rates.get(ist, r))) {
                                        initSol_i.set(0, initSol_i.getNumCols() - 1, nir.get(0, r) - sumKir_i);
                                    }
                                } else {
                                    state_i.set(0, state_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                                    if (!isNaN(sn.rates.get(ist, r))) {
                                        initSol_i.set(0, initSol_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                                    }
                                }
                            }
                        }
                        break;
                    default:
                        System.err.format("Unsupported scheduling policy at station %d.", ist);
                        return;
                }

                initSol = Matrix.concatColumns(initSol, initSol_i, null);
                sn.state.put(this.model.getStatefulNodes().get(isf), state_i);
            }
        }
        options.init_sol = initSol;
    }

    /**
     * The options the passage-time ODE runs under: the solver's own, plus the
     * closure the mean solve finished at when that solve was a moment-closure one.
     *
     * <p>The per-station variance is layout-free, so it transfers to the extended
     * (transient-class) model unchanged. The coordinate covariance block is not:
     * its rows are state coordinates, and the transient class adds one block of
     * them per station, so it is left out and the class share stays the plug-in
     * ratio. The capacity term is what separates the two closures at a busy
     * station, and that one does transfer.</p>
     *
     * @return options carrying the fixed-point closure, or the solver's own when
     *         the mean was solved by a first-order method
     */
    private SolverOptions passageTimeOptions() {
        if (!(this.result instanceof FluidResult)) {
            return options;
        }
        Matrix sigma2 = ((FluidResult) this.result).momentSigma2Drift;
        if (sigma2 == null || sigma2.isEmpty()) {
            return options;
        }
        SolverOptions ptOptions = options.copy();
        double[] s2 = new double[sigma2.getNumRows()];
        for (int i = 0; i < s2.length; i++) {
            s2[i] = sigma2.get(i, 0);
        }
        ptOptions.config.moment_sigma2 = s2;
        ptOptions.config.moment_cov = null;
        return ptOptions;
    }

    private Matrix[][] passageTime() {

        int M = sn.nstations; // Number of Stations
        int K = sn.nclasses; // Number of Classes
        double N = sn.nclosedjobs; // Population
        Matrix S = sn.nservers.copy();
        for (int i = 0; i < M; i++) {
            // Set number of servers in delay station = population
            if (Double.isInfinite(S.get(i, 0))) {
                S.set(i, 0, N);
            }
        }
        Matrix[][][] tmpRT = new Matrix[M][K][2];

        // The passage-time ODE is a SECOND solve on the fixed point the mean solve
        // reached, so it has to be driven by the same drift. Under the moment closure
        // that drift closes min(n_i,c_i) at the fixed-point variance; evaluating the
        // first-order min() there instead drains a station the mean solve holds below
        // capacity at full rate, and the response-time distribution then contradicts
        // the mean the same solver reports (measured on cdf_respt_populationsN4: CDF
        // mean 3.4504 against RN 4.5181). Only the per-STATION variance transfers:
        // it is the variance of a station population, which the transient class only
        // relabels, whereas the coordinate covariance block is indexed by a state
        // layout the transient class changes.
        SolverOptions ptOptions = passageTimeOptions();

        // Initialisation
        Matrix slowrate = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                // Service completion (exit) rates in each phase
                slowrate.set(
                        i,
                        k,
                        FastMath.min(
                                POSITIVE_INFINITY,
                                sn.mu.get(this.model.getStations().get(i)).get(sn.jobclasses.get(k)).elementMin()));
            }
        }

        // Response time analysis - starting from fixed point found
        JobClass KcClass = null;
        for (int i = 0; i < M; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) != NodeType.Source) {
                for (int k = 0; k < sn.nchains; k++) { // once for each chain
                    LinkedList<Integer> idxClassesInChain = new LinkedList<>();
                    for (int c = 0; c < sn.chains.getNumCols(); c++) {
                        if (sn.chains.get(k, c) == 1) {
                            idxClassesInChain.add(c);
                        }
                    }
                    for (Integer c : idxClassesInChain) {
                        if (sn.phases.get(i, c) > 0) {

                            // Generate ODEs Passage Time
                            Matrix phases_c;
                            int Kc = K + 1; // Add a single new class
                            // Create transient JobClass (but don't modify original sn)
                            KcClass = new JobClass(JobClassType.OPEN, "transientClass_" + c);
                            // see _kb/06-solver-catalog.md (JAR-only implementation notes: transient response-time class routing-table copy)
                            Matrix idxTranCl = new Matrix(1, K);
                            for (Integer c2 : idxClassesInChain) {
                                idxTranCl.set(0, c2, K);
                            }

                            Matrix newRT = new Matrix(M * Kc, M * Kc); // New routing table
                            Map<Station, Map<JobClass, MatrixCell>> newProc = new HashMap<>();
                            Map<Station, Map<JobClass, Matrix>> newMu = new HashMap<>();
                            Map<Station, Map<JobClass, Matrix>> newPi = new HashMap<>();

                            for (int j = 0; j < M; j++) {
                                Station station = this.model.getStations().get(j);
                                newMu.put(station, new HashMap<>());
                                newPi.put(station, new HashMap<>());
                                newProc.put(station, new HashMap<>());
                                for (int r = 0; r < K; r++) {
                                    JobClass jobClass = sn.jobclasses.get(r);
                                    newProc.get(station).put(jobClass, new MatrixCell());
                                    int procSize = sn.proc.get(station).get(jobClass).size();
                                    for (int q = 0; q < procSize; q++) {
                                        newProc
                                                .get(station)
                                                .get(jobClass)
                                                .set(q, sn.proc.get(station).get(jobClass).get(q).copy());
                                    }
                                    // Service rates
                                    newMu.get(station).put(jobClass, sn.mu.get(station).get(jobClass).copy());
                                    // Completion Probabilities
                                    newPi.get(station).put(jobClass, sn.phi.get(station).get(jobClass).copy());
                                }
                                newMu
                                        .get(station)
                                        .put(KcClass, sn.mu.get(station).get(sn.jobclasses.get(c)).copy());
                                newPi
                                        .get(station)
                                        .put(KcClass, sn.phi.get(station).get(sn.jobclasses.get(c)).copy());

                                // PHD Distribution
                                for (int r = 0; r < sn.nchains; r++) {
                                    JobClass jobClass = sn.jobclasses.get(r);
                                    int procSize = sn.proc.get(station).get(jobClass).size();
                                    for (int q = 0; q < procSize; q++) {
                                        newProc
                                                .get(station)
                                                .get(jobClass)
                                                .set(q, sn.proc.get(station).get(jobClass).get(q).copy());
                                    }
                                }
                                newProc.get(station).put(KcClass, new MatrixCell());
                                int procSize = sn.proc.get(station).get(sn.jobclasses.get(c)).size();
                                for (int q = 0; q < procSize; q++) {
                                    newProc
                                            .get(station)
                                            .get(KcClass)
                                            .set(q, sn.proc.get(station).get(sn.jobclasses.get(c)).get(q).copy());
                                }
                            }

                            // Routing/switching probabilities among basic classes
                            // MATLAB: new_rt(l:Kc:end,m:Kc:end) = rt(l:K:end,m:K:end);
                            for (int l = 0; l < K; l++) {
                                for (int m = 0; m < K; m++) {
                                    // Copy routing from class l to class m at all stations
                                    for (int si = 0; si < M; si++) {
                                        for (int sj = 0; sj < M; sj++) {
                                            // Original routing matrix index: station_i * K + class_l -> station_j * K + class_m
                                            double routingProb = sn.rt.get(si * K + l, sj * K + m);
                                            // New routing matrix index with Kc classes per station
                                            newRT.set(si * Kc + l, sj * Kc + m, routingProb);
                                        }
                                    }
                                }
                            }

                            // see _kb/06-solver-catalog.md (JAR-only implementation notes: transient response-time class routing-table copy)
                            
                            // Check if there's any routing for class c
                            double tmpSum = 0;
                            for (int si = 0; si < M; si++) {
                                for (int sj = 0; sj < M; sj++) {
                                    tmpSum += sn.rt.get(si * K + c, sj * K + c);
                                }
                            }
                            
                            if (tmpSum > 0) {
                                // Copy ALL routing for transient class from class c
                                // We'll override station i routing later in the return routing section
                                for (int si = 0; si < M; si++) {
                                    for (int sj = 0; sj < M; sj++) {
                                        // Get routing probability from station si class c to station sj class c
                                        double routingProb = sn.rt.get(si * K + c, sj * K + c);
                                        // Set same routing for transient class K to transient class K
                                        newRT.set(si * Kc + K, sj * Kc + K, routingProb);
                                    }
                                }
                            }

                            // Phases of transient classes
                            phases_c = sn.phases.copy();
                            phases_c.expandMatrix(M, Kc, sn.phases.getNumElements() + M);
                            for (int row = 0; row < M; row++) {
                                phases_c.set(row, K, sn.phases.get(row, c));
                            }

                            // see _kb/06-solver-catalog.md (JAR-only implementation notes: transient response-time class routing-table copy)
                            
                            for (Integer l : idxClassesInChain) { // For each class in chain
                                for (int j = 0; j < M; j++) {
                                    // The transient class index is K (0-indexed, so position K in the array)
                                    // At station i, transient class K completes and returns to class l at station j
                                    double routingProb = sn.rt.get(i * K + c, j * K + l);
                                    
                                    // Return fluid from transient class (at position K) to original class l
                                    newRT.set(i * Kc + K, j * Kc + l, routingProb);
                                    
                                    // Delete corresponding transition among transient classes FROM station i
                                    // (transient class K at station i should not route to transient class K at station j)
                                    newRT.set(i * Kc + K, j * Kc + K, 0);
                                }
                            }
                            
                            // Prevent other stations from routing TO the transient class at station i
                            // This prevents loops where transient class returns to station i
                            for (int j = 0; j < M; j++) {
                                if (j != i) {
                                    newRT.set(j * Kc + K, i * Kc + K, 0);
                                }
                            }

                            // Setup Initial Point and Empty Array for Solution
                            int stateLength = (int) phases_c.elementSum();
                            double[] initialState = new double[stateLength];
                            double[] nextState = new double[stateLength];
                            for (int idx = 0; idx < stateLength; idx++) {
                                initialState[idx] = 0;
                                nextState[idx] = 0;
                            }
                            double fluid_c = 0;

                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int idxNew_jl =
                                            (int) phases_c.sumSubMatrix(0, j, 0, phases_c.getNumCols())
                                                    + (int) phases_c.sumSubMatrix(j, j + 1, 0, l);
                                    int idxNew_jt =
                                            (int) phases_c.sumSubMatrix(0, j, 0, phases_c.getNumCols())
                                                    + (int) phases_c.sumSubMatrix(j, j + 1, 0, (int) idxTranCl.get(0, l));
                                    int idx_jl =
                                            (int) sn.phases.sumSubMatrix(0, j, 0, sn.phases.getNumCols())
                                                    + (int) sn.phases.sumSubMatrix(j, j + 1, 0, l);
                                    if (i == j && l == c) {
                                        // Check if init_sol is valid
                                        double mass = 0;
                                        if (options.init_sol != null && !options.init_sol.isEmpty()) {
                                            int endIdx = idx_jl + (int) sn.phases.get(j, l);
                                            if (endIdx <= options.init_sol.getNumCols()) {
                                                mass = options.init_sol.sumSubMatrix(
                                                        0, 1, idx_jl, endIdx);
                                            }
                                        }
                                        
                                        // For open networks, if there's no mass at the station, inject test fluid
                                        if (mass == 0 && sn.jobclasses.get(c).getJobClassType() == JobClassType.OPEN) {
                                            // Inject a small amount of test fluid for passage time measurement
                                            // The amount doesn't affect the CDF shape, just needs to be non-zero
                                            mass = 1.0; // Use unit mass for simplicity
                                        }
                                        
                                        initialState[idxNew_jt] = mass; // mass in phases all moved back into phase 1
                                        fluid_c += mass;
                                    } else { // Leave mass as it is
                                        int idx = idxNew_jl;
                                        for (int q = idx_jl; q < idx_jl + sn.phases.get(j, l); q++) {
                                            initialState[idx] = options.init_sol.get(0, q);
                                            idx++;
                                        }
                                    }
                                }
                            }

                            // Determine max integration time (matching MATLAB: nonZeroRates > tol)
                            double minNonZeroRate = POSITIVE_INFINITY;
                            for (int row = 0; row < M; row++) {
                                for (int col = 0; col < K; col++) {
                                    double val = slowrate.get(row, col);
                                    if (val > options.tol && val < minNonZeroRate) {
                                        minNonZeroRate = val;
                                    }
                                }
                            }
                            
                            // Solve ODE until T = 100 events with slowest exit rate
                            double T = abs(100 / minNonZeroRate);
                            double[] tRange = {0, T};
                            

                            // Indices of new classes at Station i
                            LinkedList<Integer> idxN = new LinkedList<>();
                            double end = phases_c.sumSubMatrix(i, i + 1, K, Kc);
                            for (int idx = 0; idx < end; idx++) {
                                idxN.add(
                                        idx
                                                + (int) phases_c.sumSubMatrix(0, i, 0, phases_c.getNumCols())
                                                + (int) phases_c.sumSubMatrix(i, i + 1, 0, K));
                            }

                            // Create an extended NetworkStruct with the transient class
                            NetworkStruct extendedSn = new NetworkStruct();
                            // Copy basic fields from original sn
                            extendedSn.nstations = sn.nstations;
                            extendedSn.nclasses = Kc;  // Extended number of classes
                            extendedSn.nservers = sn.nservers.copy();
                            extendedSn.nclosedjobs = sn.nclosedjobs;
                            extendedSn.stations = sn.stations;
                            extendedSn.sched = sn.sched;
                            extendedSn.schedparam = sn.schedparam.copy();
                            
                            // Create extended job classes list that includes the transient class
                            extendedSn.jobclasses = new ArrayList<>(sn.jobclasses);
                            extendedSn.jobclasses.add(KcClass);  // Add the transient class
                            
                            // Set-up the ODEs for the new QN with extended structure
                            FirstOrderDifferentialEquations ode =
                                    new PassageTimeODE(
                                            extendedSn, newMu, newPi, newProc, newRT, S, ptOptions,
                                            initialState.length);

                            // ODE analysis
                            // The window loop advances initialState, so the CDF
                            // refinement below keeps the state the trajectory
                            // STARTED from: it re-integrates the whole curve on
                            // the refined grid.
                            double[] passageStart = initialState.clone();
                            Matrix tFull = new Matrix(0, 0);
                            Matrix stateFull = new Matrix(0, 0);
                            int iter = 1;
                            boolean finished = false;
                            double tref = 0;
                            boolean stiff = options.stiff;
                            // THE CDF IS READ BACK BY RIEMANN-STIELTJES QUADRATURE, so the
                            // undepleted transient fluid the integrator is entitled to leave
                            // is weighted by t for the mean and by t^2 for the second moment.
                            // At the generic tol that residual is O(1e-4) while the horizon is
                            // O(1e2), so the SECOND moment inherits O(1) of pure integration
                            // noise: measured on cdf_respt_closed_threeclasses, 7.6e-6 of fluid
                            // still undepleted past t=20 put 0.022 into m2 and read the Exp(1)
                            // response-time SCV as 1.0416 against the reference 1.0101. MATLAB's
                            // ode15s and the C++ LSODA both settle far below their stated atol
                            // on this model; LSODAExt settles ON it, so the passage solve has to
                            // ask for more than the mean solve does. A caller asking for tighter
                            // still gets it.
                            double ptTol = Math.min(options.tol, 1e-6);

                            while (iter <= options.iter_max && !finished) {

                                // Use LSODA (stiff solver) for passage time ODE, matching MATLAB ode15s.
                                // Use tight tolerances for accurate CDF representation.
                                Matrix tmpTIter;
                                Matrix tmpStateIter;

                                try {
                                    // Use LSODA for passage time ODE
                                    // Use small min step (1e-12) to handle stiff initial conditions.
                                    // The step budget is sized to separate a converging solve from
                                    // a thrashing one: across the passage times of an LQN layer
                                    // every convergent solve took 45-530 steps, whereas a stiff one
                                    // collapses its step size and runs to the cap. At the former
                                    // 1e7 that cost minutes per solve before throwing, only for the
                                    // caller to fall back to the exponential approximation anyway;
                                    // 1e5 keeps a ~200x margin over the observed worst case and
                                    // reaches the same fallback in seconds.
                                    LSODAExt lsodaSolver = new LSODAExt(
                                            1e-12, options.odesolvers.odemaxstep,
                                            ptTol, ptTol, 12, 5, 100000);
                                    double[] tempState = new double[initialState.length];
                                    lsodaSolver.integrate(ode, 0, initialState, T, tempState);

                                    int numSteps = lsodaSolver.getStepsTaken() + 1;
                                    ArrayList<Double> tHistory = lsodaSolver.getTvec();
                                    ArrayList<Double[]> yHistory = lsodaSolver.getYvec();

                                    // Use raw LSODA step points with non-negative clamping
                                    tmpTIter = new Matrix(numSteps, 1);
                                    tmpStateIter = new Matrix(numSteps, initialState.length);
                                    for (int step = 0; step < numSteps; step++) {
                                        tmpTIter.set(step, 0, tHistory.get(step));
                                        Double[] yStep = yHistory.get(step);
                                        for (int j = 0; j < initialState.length; j++) {
                                            tmpStateIter.set(step, j, Math.max(0, yStep[j]));
                                        }
                                    }

                                    System.arraycopy(tempState, 0, nextState, 0, tempState.length);

                                } catch (RuntimeException e) {
                                    // chain the cause: dropping it here left the ODE failure with no trace
                                    throw new RuntimeException("ODE Solver Failed: " + e.getMessage(), e);
                                }

                                iter++;

                                if (tFull.isEmpty()) {
                                    tFull = tmpTIter.copy();
                                    stateFull = tmpStateIter.copy();
                                } else {
                                    int startRow = tFull.getNumRows();
                                    tFull.expandMatrix(
                                            tFull.getNumRows() + tmpTIter.getNumRows(),
                                            1,
                                            tFull.getNumElements() + tmpTIter.getNumElements());
                                    int endRow = tFull.getNumRows();
                                    for (int row = startRow; row < endRow; row++) {
                                        tFull.set(row, 0, tmpTIter.get(row - startRow, 0) + tref);
                                    }
                                    stateFull = Matrix.concatRows(stateFull, tmpStateIter, null);
                                }

                                tmpSum = 0;
                                int tmpTMax = tmpStateIter.getNumRows();
                                for (Integer idx : idxN) {
                                    tmpSum += tmpStateIter.get(tmpTMax - 1, idx);
                                }
                                if (tmpSum < 0.000000001) {
                                    finished = true;
                                }
                                tref += tmpTIter.get(tmpTMax - 1, 0);
                                for (int idx = 0; idx < initialState.length; idx++) {
                                    initialState[idx] = tmpStateIter.get(tmpTMax - 1, idx);
                                }

                            }

                            
                            // Retrieve response time CDF for class k
                            int fullTMax = tFull.getNumRows();
                            
                            tmpRT[i][c][0] = tFull;
                            tmpRT[i][c][1] = tFull.copy();
                            if (fluid_c > 0) {
                                for (int row = 0; row < fullTMax; row++) {
                                    tmpSum = 0;
                                    for (Integer idx : idxN) {
                                        tmpSum += stateFull.get(row, idx);
                                    }
                                    double cdfValue = 1 - tmpSum / fluid_c;
                                    tmpRT[i][c][1].set(row, 0, cdfValue);
                                }
                                
                                // Iterative CDF refinement - detect and refine large CDF jumps.
                                // The rule of MATLAB solver_fluid_passage_time.m, the C++
                                // fluid_passage_time and the Python passage_time: while some
                                // adjacent pair of CDF values differs by more than maxCdfJump,
                                // split EVERY offending interval and re-integrate the whole
                                // curve on the new grid. Refining ONE interval per round spends
                                // the round cap on five intervals and leaves the jump target
                                // unmet, and the curve is read back by quadrature: on
                                // cdf_respt_closed_threeclasses a 130-point grid over [0,200]
                                // made the right-endpoint mean read 1.126 for a response time
                                // that is exactly Exp(1). Both caps bound WORK, not accuracy.
                                double maxCdfJump = 0.0005;
                                int maxRefinementRounds = 5;
                                int maxPoints = 20001;
                                int numRefinedPoints = 20;

                                // The refinement re-integrates the SAME passage, so it answers to
                                // ptTol as well: leaving it on the mean solve's tol would put the
                                // residual back into the curve one round after the first pass
                                // removed it. Both helpers hand back the caller's own integrator
                                // untouched when one was supplied.
                                FirstOrderIntegrator refineOdeSolver;
                                if (options.stiff) {
                                    refineOdeSolver = options.odesolvers.stiffIntegratorFor(0, T, ptTol, false);
                                } else {
                                    refineOdeSolver = options.odesolvers.integratorFor(0, T, ptTol, false);
                                }

                                for (int round = 0; round < maxRefinementRounds; round++) {
                                    if (fullTMax >= maxPoints) {
                                        break;
                                    }
                                    // A jump at EQUAL times is an atom of the law, not a
                                    // resolution failure: its refined points would all be
                                    // the same instant.
                                    List<Double> grid = new ArrayList<Double>();
                                    boolean refined = false;
                                    for (int row = 0; row + 1 < fullTMax; row++) {
                                        double t1 = tFull.get(row, 0);
                                        double t2 = tFull.get(row + 1, 0);
                                        grid.add(t1);
                                        double jump = tmpRT[i][c][1].get(row + 1, 0) - tmpRT[i][c][1].get(row, 0);
                                        if (jump > maxCdfJump && t2 > t1) {
                                            refined = true;
                                            for (int rp = 1; rp < numRefinedPoints; rp++) {
                                                grid.add(t1 + (t2 - t1) * rp / (double) numRefinedPoints);
                                            }
                                        }
                                    }
                                    if (fullTMax > 0) {
                                        grid.add(tFull.get(fullTMax - 1, 0));
                                    }
                                    if (!refined || grid.size() > maxPoints || grid.size() < 2) {
                                        break;
                                    }

                                    Matrix newTFull = new Matrix(grid.size(), 1);
                                    Matrix newStateFull = new Matrix(grid.size(), passageStart.length);
                                    double[] walk = passageStart.clone();
                                    double tWalk = grid.get(0);
                                    boolean failed = false;
                                    for (int g = 0; g < grid.size(); g++) {
                                        double tg = grid.get(g);
                                        if (tg > tWalk) {
                                            double[] out = new double[passageStart.length];
                                            try {
                                                refineOdeSolver.integrate(ode, tWalk, walk, tg, out);
                                            } catch (Exception e) {
                                                failed = true;
                                                break;
                                            }
                                            walk = out;
                                            tWalk = tg;
                                        }
                                        newTFull.set(g, 0, tg);
                                        for (int j = 0; j < passageStart.length; j++) {
                                            newStateFull.set(g, j, Math.max(0, walk[j]));
                                        }
                                    }
                                    if (failed) {
                                        break;
                                    }

                                    tFull = newTFull;
                                    stateFull = newStateFull;
                                    fullTMax = tFull.getNumRows();
                                    tmpRT[i][c][0] = tFull;
                                    tmpRT[i][c][1] = tFull.copy();
                                    for (int cdfRow = 0; cdfRow < fullTMax; cdfRow++) {
                                        tmpSum = 0;
                                        for (Integer idx : idxN) {
                                            tmpSum += stateFull.get(cdfRow, idx);
                                        }
                                        tmpRT[i][c][1].set(cdfRow, 0, 1 - tmpSum / fluid_c);
                                    }
                                }

                                // Horizon extension - extend while the TAIL misses the law.
                                // The test is on the LAST grid point. It used to read the
                                // FIRST, tmpRT[i][c][1].get(0,0), which is the CDF at the start
                                // of the horizon: the marked class holds all of fluid_c at t=0
                                // by construction, so that value is 0 whatever the horizon is,
                                // and lengthening the horizon cannot move it. The loop it
                                // guarded was unreachable, and reachable only into harm -- its
                                // body REPLACED the refined curve with raw LSODA steps over a
                                // fresh horizon, discarding the grid the refinement rounds
                                // above had just paid for. Same fix as MATLAB
                                // solver_fluid_passage_time.m.
                                int extendIterations = 0;
                                final int maxExtendIterations = 10; // Prevent infinite loops

                                while (tmpRT[i][c][1].get(fullTMax - 1, 0) < 0.99
                                        && extendIterations < maxExtendIterations) {
                                    extendIterations++;

                                    // CONTINUE the same trajectory from where it stopped and
                                    // APPEND, as the window loop above does: initialState is
                                    // the end state and tref the elapsed time, and the ODE is
                                    // autonomous, so [0, extendedT] from it is the next stretch
                                    // of the SAME passage. Doubling each round reaches a 1024x
                                    // horizon within the cap instead of 11x.
                                    double extendedT = T * Math.pow(2, extendIterations);

                                    Matrix extTIter;
                                    Matrix extStateIter;
                                    try {
                                        LSODAExt extLsoda = new LSODAExt(
                                                1e-12, options.odesolvers.odemaxstep,
                                                ptTol, ptTol, 12, 5, 10000000);
                                        double[] extTemp = new double[initialState.length];
                                        extLsoda.integrate(ode, 0, initialState, extendedT, extTemp);

                                        int extSteps = extLsoda.getStepsTaken() + 1;
                                        ArrayList<Double> extTHist = extLsoda.getTvec();
                                        ArrayList<Double[]> extYHist = extLsoda.getYvec();

                                        // Use raw LSODA steps (no densification)
                                        extTIter = new Matrix(extSteps, 1);
                                        extStateIter = new Matrix(extSteps, initialState.length);
                                        for (int step = 0; step < extSteps; step++) {
                                            extTIter.set(step, 0, extTHist.get(step));
                                            Double[] extYStep = extYHist.get(step);
                                            for (int j = 0; j < initialState.length; j++) {
                                                extStateIter.set(step, j, Math.max(0, extYStep[j]));
                                            }
                                        }
                                    } catch (Exception e) {
                                        break;
                                    }
                                    int extRows = extTIter.getNumRows();
                                    if (extRows < 2) {
                                        break;
                                    }

                                    // drop the duplicated first row: it repeats the instant the
                                    // curve already ends on
                                    Matrix appendT = new Matrix(extRows - 1, 1);
                                    Matrix appendState = new Matrix(extRows - 1, initialState.length);
                                    for (int row = 1; row < extRows; row++) {
                                        appendT.set(row - 1, 0, extTIter.get(row, 0) + tref);
                                        for (int j = 0; j < initialState.length; j++) {
                                            appendState.set(row - 1, j, extStateIter.get(row, j));
                                        }
                                    }
                                    tFull = Matrix.concatRows(tFull, appendT, null);
                                    stateFull = Matrix.concatRows(stateFull, appendState, null);
                                    fullTMax = tFull.getNumRows();

                                    tref += extTIter.get(extRows - 1, 0);
                                    for (int idx = 0; idx < initialState.length; idx++) {
                                        initialState[idx] = extStateIter.get(extRows - 1, idx);
                                    }

                                    tmpRT[i][c][0] = tFull;
                                    tmpRT[i][c][1] = tFull.copy();
                                    for (int row = 0; row < fullTMax; row++) {
                                        tmpSum = 0;
                                        for (Integer idx : idxN) {
                                            tmpSum += stateFull.get(row, idx);
                                        }
                                        tmpRT[i][c][1].set(row, 0, 1 - tmpSum / fluid_c);
                                    }
                                }
                            } else {
                                tmpRT[i][c][1].ones();
                            }

                            if (iter > options.iter_max) {
                                line_warning("SolverFluid",
                                        "Maximum number of iterations reached when computing the response time distribution. "
                                                + "Response time distributions may be inaccurate. Increase option.iter_max (currently at %d).",
                                        options.iter_max);
                            }
                        }
                    }
                }
            }
        }

        if (KcClass != null) {
            sn.jobclasses.remove(KcClass);
        }

        Matrix[][] RTret = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            // Skip Source stations as in MATLAB (line 119)
            NodeType nodeType = sn.nodetype.get((int) sn.stationToNode.get(i));
            if (nodeType != NodeType.Source) {
                for (int c = 0; c < K; c++) {
                    if (tmpRT[i][c] != null && tmpRT[i][c][0] != null && tmpRT[i][c][1] != null) {
                        // Check if first point is (0,0) and skip it if so
                        Matrix timeMatrix = tmpRT[i][c][0];
                        Matrix cdfMatrix = tmpRT[i][c][1];
                        
                        if (timeMatrix.getNumRows() > 1 && 
                            Math.abs(timeMatrix.get(0, 0)) < GlobalConstants.FineTol && 
                            Math.abs(cdfMatrix.get(0, 0)) < GlobalConstants.FineTol) {
                            // Skip the first (0,0) point
                            int numRows = timeMatrix.getNumRows() - 1;
                            Matrix trimmedTime = new Matrix(numRows, 1);
                            Matrix trimmedCdf = new Matrix(numRows, 1);
                            
                            for (int row = 0; row < numRows; row++) {
                                trimmedTime.set(row, 0, timeMatrix.get(row + 1, 0));
                                trimmedCdf.set(row, 0, cdfMatrix.get(row + 1, 0));
                            }
                            
                            RTret[i][c] = Matrix.concatColumns(trimmedCdf, trimmedTime, null);
                        } else {
                            // Keep original behavior if not starting with (0,0)
                            RTret[i][c] = Matrix.concatColumns(cdfMatrix, timeMatrix, null);
                        }
                    }
                    // Leave RTret[i][c] as null when tmpRT[i][c] is null, matching MATLAB behavior when RT{i,c} is empty
                }
            }
        }
        return RTret;
    }

    @Override
    public boolean supportsTransientAnalysis() {
        // Transient averages are available (fluid ODE integrated over options.timespan).
        return true;
    }

    /**
     * Runs the fluid analyzer to solve the queueing network.
     * This method executes the fluid approximation algorithm and stores
     * the results in the solver's result object.
     */
    @Override
    public void runAnalyzer() {
        // A MODEL HOLDING ANY Transition NODE IS A DIFFERENT FORMALISM and goes
        // to a dedicated runner rather than through the queueing dispatch. The
        // post-processing of the queueing path rewrites UN from sn.rates (NaN at
        // a Place) and RN from the station scheduling, and would overwrite the
        // Petri conventions PetriSolver reports -- a place's utilization IS its
        // token count, and its throughput is the token departure rate, not a
        // service completion rate. The branch is taken before
        // runAnalyzerChecks, whose feature gate is written about queueing
        // stations. See PetriSolver, and matlab solver_fluid_petri.m.
        if (isPetriNet(this.model.getStruct(false))) {
            runPetriAnalyzer();
            return;
        }
        // Validate model compatibility before starting analysis
        runAnalyzerChecks(this.options);

        // Finite Capacity Region: the fluid ODEs do not enforce the aggregate
        // per-region job limit and would silently return the unconstrained answer.
        //
        // 'dae' IS THE EXCEPTION, AND THE ONLY ONE. A region cap is a linear
        // inequality on the state and blocking is a throttle on the admission
        // flow that keeps it satisfied, so the DAE form has somewhere to put it
        // -- an algebraic equation beside the drift -- where an ODE has not.
        // DaeAnalyzer refuses, by name, the forms that are NOT constraints on
        // this drift (BAS/BBS/RSRD, retrial, per-class admission weights). Every
        // other method keeps the blanket refusal, because for them it is still
        // true.
        // A BINDING BUFFER OR A CAPACITY REGION RESOLVES 'default' TO 'dae'
        // before either gate below reads the method, because 'dae' is the only
        // fluid route that carries the constraint at all: the reference resolves
        // 'default' ahead of its own gate for the same reason. Every other
        // resolution stays where it is, in the dispatch switch, which cannot run
        // before the struct is phase-type converted.
        if ("default".equals(fluidUnqualify(this.options.method))
                && blockedResolvesToDae(this.model, this.model.getStruct(false), this.options)) {
            this.options.method = "dae";
            line_debug(options.verbose, "FLD default resolved to dae: the model has a binding "
                    + "finite buffer or capacity region");
        }
        boolean daeMethod = "dae".equals(this.options.method)
                || "fluid.dae".equals(this.options.method);
        // 'mol' is exempt from the CAPACITY gate for the opposite reason to
        // 'dae': a finite capacity is not something it ignores, it IS the
        // model. The approximation is stated for the Mt/G/s/0 LOSS system, so
        // the server count is the buffer, and refusing a capped station would
        // refuse the only shape the method answers. It is NOT exempt from the
        // REGION gate above, which constrains a SET of stations and has no
        // counterpart in a single-station limit. MATLAB
        // (@@SolverFLD/runAnalyzer.m) and native python draw the line in the
        // same place; this port exempted 'dae' alone.
        boolean capExempt = daeMethod
                || "mol".equals(this.options.method) || "fluid.mol".equals(this.options.method);
        if (!daeMethod && this.model.getStruct(false).nregions > 0) {
            throw new RuntimeException("This model uses a Finite Capacity Region (addRegion), "
                    + "which is not supported by SolverFluid (the region's aggregate job limit "
                    + "is not enforced). Use options.method = \"dae\", or SolverCTMC, "
                    + "SolverJMT, SolverSSA or SolverLDES.");
        }

        // A BINDING STATION BUFFER WAS SILENTLY IGNORED, by every fluid method
        // including this one: nothing in the fluid tree reads sn.cap or sn.classcap,
        // so a capped station was integrated as an unbounded one and the table
        // reported more jobs in the buffer than the buffer holds. MVA and NC have
        // refused that model through the shared structural gate since they gained
        // one; the fluid solver now does too, except on the route that can enforce
        // it as an algebraic constraint on the drift.
        if (!capExempt) {
            String capReason = NetworkSolver.bindingCapacityReason(this.model,
                    this.model.getStruct(false), "SolverFluid");
            if (capReason != null) {
                throw new RuntimeException(capReason + " Use options.method = \"dae\", which "
                        + "carries the buffer as an algebraic constraint on the drift.");
            }
        }

        long startTime = System.nanoTime();
        line_debug(options.verbose, String.format("Fluid solver starting: method=%s, iterations=%d", 
            options.method, options.iter_max));
        
        // Store the original method before any modifications
        // Unqualify the method BEFORE the dispatch. Every method name has a
        // 'fluid.'-qualified spelling, which the reference accepts arm by arm
        // (erase(options.method,'fluid.')) and which the featset and the
        // dispatch switches here matched only for a handful; 'butools' names the
        // MFQ backend and 'aoi' its age-of-information reading, and both are
        // aliases of 'mfq'. Doing it once, here, is what keeps the two switches
        // below and getMethodFeatureSet from disagreeing on a qualified name.
        options.method = fluidUnqualify(options.method);
        String origMethod = options.method;
        
        // see _kb/06-solver-catalog.md (JAR-only implementation notes: non-Markovian to phase-type conversion before initSol)
        // the fluid ODEs read mu*phi as a flow, so the surrogate must be a genuine phase-type
        String phfit0 = this.options.config.phfit;
        this.options.config.phfit = "ph";
        sn = SnNonmarkovToPh.snNonmarkovToPh(this.model.getStruct(), this.options, false);
        this.options.config.phfit = phfit0;

        // Fork-join: the same solver-agnostic fixed point MVA and NC drive
        // (jline.solvers.fj.FJFixedPoint), with fldDispatch as the inner solve.
        // The MMT transformation emits only Source, Delay, Queue, Router and
        // ClassSwitch, all of which the fluid drift already carries.
        if (this.model.hasFork()) {
            // Not every method can run that fixed point. forkJoinAdmitsReason is
            // the same predicate supportsModelMethod asks, so the report and the
            // run cannot disagree about which forks this method serves.
            String fjReason = forkJoinAdmitsReason(sn, this.options.method);
            if (!fjReason.isEmpty()) {
                throw new RuntimeException(fjReason);
            }
            runForkJoinAnalyzer(startTime);
            return;
        }

        boolean hasOpenClasses = false;
        for (NodeType nodetype : sn.nodetype) {
            if (nodetype == NodeType.Source) {
                hasOpenClasses = true;
                break;
            }
        }
        boolean hasDPS = false;
        for (SchedStrategy sched : sn.sched.values()) {
            if (sched == SchedStrategy.DPS) {
                hasDPS = true;
                break;
            }
        }
        boolean hasCache = false;
        for (NodeType nodetype : sn.nodetype) {
            if (nodetype == NodeType.Cache) {
                hasCache = true;
                break;
            }
        }

        String actualMethod = origMethod;
        switch (origMethod) {
            case "pnorm":
                // The p-norm smoothing of Ruuskanen et al., PEVA 151 (2021): the
                // matrix drift with the hard min replaced by a p-norm, which
                // MatrixMethodODE already builds when config.pstar is non-empty.
                // pstar defaults to 20, the value solver_fluid_odes.m uses and
                // the one whose behaviour matches softmin at alpha = 20. The
                // method was neither listed nor dispatched here although the ODE
                // has always carried the smoothing.
                if (options.config.pstar == null) {
                    options.config.pstar = new java.util.ArrayList<Double>();
                }
                if (options.config.pstar.isEmpty()) {
                    options.config.pstar.add(20.0);
                }
                if (hasDPS) {
                    if (options.verbose != VerboseLevel.SILENT) {
                        line_error(mfilename(new Object(){}),
                                "The matrix solver does not support DPS scheduling. Using options.method = \"closing\" instead.");
                    }
                    actualMethod = "closing";
                    options.method = actualMethod;
                }
                break;
            case "matrix":
                if (hasDPS) {
                    if (options.verbose != VerboseLevel.SILENT) {
                        line_error(mfilename(new Object(){}),
                                "The matrix solver does not support DPS scheduling. Using options.method = \"closing\" instead.");
                    }
                    actualMethod = "closing";
                }
                break;
            case "default":
                // A blocked model never reaches here: runAnalyzer resolved it to
                // "dae" above, ahead of the finite-capacity gate that would
                // otherwise have refused whatever this arm picked.
                //
                // Preference order: rmf for cache models, then minnormal wherever it
                // applies, then the historical closing-for-DPS / matrix. The
                // second-order closure dominates the first-order methods on every
                // family measured against exact CTMC and is the only method that can
                // represent GPS at all.
                if (hasCache) {
                    actualMethod = "rmf";
                } else if (FluidMinNormalApplicable.reasonToDecline(sn, options) == null) {
                    actualMethod = "minnormal";
                } else if (hasDPS) {
                    actualMethod = "closing";
                } else {
                    actualMethod = "matrix";
                }
                options.method = actualMethod;
                break;
            case "closing":
            case "diffusion":
            case "refined":
            case "statedep":
            // 'softmin' is ode_statedep with the hard min replaced by the smooth
            // one (PassageTimeODE already builds that drift, alpha = 20). It was
            // advertised by listValidMethods and had no arm here, so it fell to
            // the default: line_error below and a listed method always threw.
            case "softmin":
            case "mfq":
            case "rmf":
            case "minnormal":
            case "ggisgi.fluid":
            case "fluid.ggisgi":
            case "ggisgi":
            case "ggingi.tga":
            case "fluid.tga":
            case "tga":
            case "tvms":
            case "fluid.tvms":
            case "mtginf":
            case "fluid.mtginf":
            case "mol":
            case "fluid.mol":
                break;
            case "dae":
                // The DAE route inherits the moment-closure envelope, with two
                // extra limits refused here so the message names the model
                // feature rather than surfacing from inside the Newton solve.
                //
                // A CACHE MODEL IS A DECOMPOSITION, not one system: the caches
                // are solved in isolation and the network with them relabeled,
                // so there is no single drift for the constraint to be attached
                // to. The 'rmf' alternation carries the closure inside its
                // network step for 'minnormal'; no such route exists for the
                // DAE form.
                if (hasCacheNodes(sn)) {
                    line_error(mfilename(new Object(){}),
                            "The dae method does not support caching stations: a cache model is solved by decomposition, so it has no single drift to constrain. Use options.method = \"minnormal\" for the same closure, or \"rmf\".");
                }
                break;
            case "kp":
                // Ko-Pender limits are proved for an OPEN network fed by external arrival
                // processes: a closed class has no arrival process to modulate.
                if (!hasOpenClasses) {
                    line_error(mfilename(new Object(){}),
                            "The kp method analyses the open (MAP_t/Ph_t/inf)^N network of Ko and Pender (2017); a closed class has no arrival process to modulate. Use options.method = \"closing\".");
                }
                break;
            case "tbi":
                // Trajectory-based iteration is a closed-network transient
                // method: the Jacobi waveform relaxation decomposes the closed
                // population into cells, which requires a bounded state and no
                // cache-driven state-dependent routing.
                if (hasOpenClasses) {
                    line_error(mfilename(new Object(){}),
                            "The tbi method does not support open classes. Use a closed model or options.method = \"closing\".");
                }
                if (hasCache) {
                    line_error(mfilename(new Object(){}),
                            "The tbi method does not support cache nodes. Use options.method = \"rmf\".");
                }
                break;
            default:
                line_error(mfilename(new Object(){}),
                        "SolverFluid does not support the specified method. Using options.method = \"default\".");
                actualMethod = "default";
        }
        result.method = actualMethod;

        // Method-aware feature gate, applied AFTER "default" has been resolved.
        // runAnalyzerChecks above can only see the coarse solver envelope, and
        // "default" is not "minnormal", so gating there would reject a GPS model
        // before the resolution ever ran. Keeping the decision and the gate in this
        // order is what stops them from disagreeing.
        // NetworkSolver.model shadows Solver.model, so the base supportsModelMethod
        // would dereference a null; resolve the feature set against this solver's model
        String methodReason =
                FeatureSet.supportsReason(getMethodFeatureSet(actualMethod), this.model.getUsedLangFeatures());
        if (methodReason != null && !methodReason.isEmpty()) {
            throw new RuntimeException("The '" + actualMethod
                    + "' method of the Fluid solver does not support this model: " + methodReason);
        }

        if (isInfinite(options.timespan[0])) {
            if (options.verbose == VerboseLevel.DEBUG) {
                line_warning(mfilename(new Object(){}),
                        "SolverFluid requires options.timespan[0] to be finite. Setting it to 0.");
            }
            options.timespan[0] = 0;
        }
        if (options.timespan[0] == options.timespan[1]) {
            line_warning(mfilename(new Object(){}),
                    "SolverFluid does not support a timespan that is a single point. Setting options.timespan[0] to 0.");
            options.timespan[0] = 0;
        }
        // THE COARSE GATE AGAIN, and deferred for 'dae' for the same reason as in
        // runAnalyzerChecks: supports() reads the STATIC getFeatureSet, while
        // getMethodFeatureSet is strictly finer and has already run. The two can
        // only disagree where a method WIDENS the set, and 'dae' is the one that
        // does -- it declares Region, which the static set must keep false so
        // that every other method goes on rejecting a finite capacity region.
        boolean daeWidens = "dae".equals(options.method) || "fluid.dae".equals(options.method)
                || QsysLimitAnalyzer.handles(options.method);
        if (this.enableChecks && !daeWidens && !supports(this.model)) {
            line_error(mfilename(new Object(){}), "This model contains features not supported by the solver.");
            return;
        }

        int M = sn.nstations;
        int K = sn.nclasses;

        // The single-station fluid limits are closed forms, not integrations of
        // the network drift: they take the whole model in one call and have no
        // initial state to average over.
        if (QsysLimitAnalyzer.handles(actualMethod)) {
            FluidAnalyzer qsysAnalyzer = new QsysLimitAnalyzer(actualMethod);
            qsysAnalyzer.analyze(sn, options.copy(), result);
            ((FluidResult) this.result).odeStateVec = qsysAnalyzer.getXVecIt();
            return;
        }

        // MFQ IS A SINGLE-QUEUE METHOD AND FALLS BACK, which is what the reference
        // does: solver_fluid_analyzer.m warns "MFQ not applicable: ... Falling
        // back to matrix method" and re-enters solver_fluid_matrix. Refusing
        // instead made 'mfq' -- and therefore its aliases 'butools' and 'aoi' --
        // reject every multi-station model that MATLAB, native python and C++ all
        // answer. The substitution has to be decided HERE and not inside
        // MFQAnalyzer: the branch below diverts around the state-space
        // preparation that the matrix method reads back (options.init_sol above
        // all), so delegating from within the analyzer indexes past it.
        if (Objects.equals(actualMethod, "mfq")) {
            String why = jline.solvers.fluid.analyzers.MFQAnalyzer.mfqNotApplicableReason(sn);
            if (why != null) {
                line_warning(mfilename(new Object() {
                }), "MFQ not applicable: %s. Falling back to matrix method.", why);
                actualMethod = "matrix";
                options.method = "matrix";
            }
        }

        // MFQ is a direct steady-state method that doesn't need state space iteration
        if (Objects.equals(actualMethod, "mfq")) {
            FluidAnalyzer analyzer = new MFQAnalyzer();
            analyzer.analyze(sn, options.copy(), result);
            ((FluidResult) this.result).odeStateVec = analyzer.getXVecIt();
            ((FluidResult) this.result).snFinal = this.sn;
            result.method = "mfq";
            // Run AoI analysis if this is a valid AoI topology
            runAoiAnalysis(sn);
            return;
        }

        // The decomposition route for cache-queueing networks. "rmf" solves the
        // network layer with the first-order matrix method; "minnormal" solves
        // the SAME decomposition with the moment closure in its place, which is
        // how the closure reaches a cache model at all (one ODE cannot express
        // it: the caches are solved in isolation and the network with them
        // relabeled as class switches).
        boolean cacheClosure = Objects.equals(actualMethod, "minnormal") && hasCacheNodes(sn);
        if (Objects.equals(actualMethod, "rmf") || cacheClosure) {
            // see _kb/06-solver-catalog.md (JAR-only implementation notes: sn.rt is a Java reference, save/restore around cache rewrite)
            Matrix rtOrig = sn.rt != null ? sn.rt.copy() : null;

            RMFAnalyzer analyzer = new RMFAnalyzer();
            SolverOptions rmfOptions = options.copy();
            rmfOptions.method = cacheClosure ? "minnormal" : "rmf";
            // THIS ROUTE HAS ITS OWN LADDER, because it returns before the one
            // below and would otherwise let a FluidNonHyperbolicException out
            // of runAnalyzer as a LineException -- reporting a fallback in the
            // message that never happened. The rungs are NOT the general
            // ladder's: dae has no decomposition arm (FluidDaeApplicable
            // declines a cache model by name), so the closure's only fallback
            // here is the first-order network step of the SAME alternation,
            // which is exactly what "rmf" is.
            try {
                analyzer.analyze(sn, rmfOptions, result);
            } catch (FluidNonHyperbolicException e) {
                if (!cacheClosure) {
                    throw e;
                }
                line_debug(options.verbose,
                        "Fluid minnormal declined inside the cache decomposition (" + e.getMessage()
                                + "); falling back to rmf");
                cacheClosure = false;
                rmfOptions.method = "rmf";
                options.method = "rmf";
                analyzer = new RMFAnalyzer();
                analyzer.analyze(sn, rmfOptions, result);
            }
            ((FluidResult) this.result).odeStateVec = analyzer.getXVecIt();
            ((FluidResult) this.result).snFinal = this.sn;
            result.method = cacheClosure ? "minnormal" : "rmf";
            if (cacheClosure && analyzer.lastMinNormal != null) {
                MinNormalAnalyzer mn = analyzer.lastMinNormal;
                FluidResult fr = (FluidResult) this.result;
                fr.momentSigma = mn.sigmaMatrix;
                fr.momentQVar = mn.qVar;
                fr.momentClassBlock = mn.classBlock;
                fr.momentStationBlock = mn.stationBlock;
                fr.momentSigma2 = sigma2Matrix(mn.sigma2);
                fr.momentSigma2Drift = sigma2Matrix(mn.sigma2Drift);
                fr.momentOuterIters = mn.outerIters;
            }

            // Store hit/miss probs on cache nodes
            FluidResult fluidResult = (FluidResult) this.result;
            if (fluidResult.hitProb != null) {
                for (int ind = 0; ind < this.model.getNodes().size(); ind++) {
                    if (this.model.getNodes().get(ind) instanceof Cache) {
                        Cache cache = (Cache) this.model.getNodes().get(ind);
                        cache.setResultHitProb(Matrix.extractRows(fluidResult.hitProb, ind, ind + 1, null));
                        cache.setResultMissProb(Matrix.extractRows(fluidResult.missProb, ind, ind + 1, null));
                    }
                }
                this.model.refreshStruct(true);
            }

            // see _kb/06-solver-catalog.md (JAR-only implementation notes: cache actual hit/miss propagation to sn.nodeparam)
            for (int ind = 0; ind < this.sn.nnodes; ind++) {
                if (this.sn.nodetype.get(ind) == jline.lang.constant.NodeType.Cache) {
                    Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                    jline.lang.nodeparam.CacheNodeParam cacheParam =
                            (jline.lang.nodeparam.CacheNodeParam) this.sn.nodeparam.get(cacheNode);
                    if (cacheParam != null) {
                        Matrix hitProb = cacheNode.getHitRatio();
                        Matrix missProb = cacheNode.getMissRatio();
                        if (hitProb != null && !hitProb.isEmpty()) {
                            cacheParam.actualhitprob = hitProb;
                        }
                        if (missProb != null && !missProb.isEmpty()) {
                            cacheParam.actualmissprob = missProb;
                        }
                    }
                }
            }

            // Compute arrival rates using the ORIGINAL rt (before cache probs were baked in).
            // snGetArvRFromTput separately applies cache probs via TN_stateful, so using
            // the refreshed rt (which already has cache probs) would double-count them.
            Matrix rtRefreshed = sn.rt;
            if (rtOrig != null) {
                sn.rt = rtOrig;
            }
            AvgHandle TH = getAvgTputHandles();
            Matrix AN = snGetArvRFromTput(sn, result.TN, TH);
            if (rtRefreshed != null) {
                sn.rt = rtRefreshed;
            }

            // Handle the "default" method case - add "default/" prefix if original method was "default"
            String finalMethod;
            if (origMethod.equals("default") && !options.method.equals("default")) {
                finalMethod = "default/" + options.method;
            } else {
                finalMethod = options.method;
            }
            result.method = finalMethod;
            this.setAvgResults(result.QN, result.UN, result.RN, result.TN, AN, new Matrix(0, 0), result.CN, result.XN, result.runtime, finalMethod, result.iter);
            return;
        }

        Matrix Q = new Matrix(M, K);
        Matrix U = Q.copy();
        Matrix R = Q.copy();
        Matrix T = Q.copy();
        Matrix C = new Matrix(1, K);
        Matrix X = C.copy();
        Matrix[][] Qt = new Matrix[M][K];
        Matrix[][] Ut = new Matrix[M][K];
        Matrix[][] Tt = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                Qt[i][k] = new Matrix(0, 0);
                Ut[i][k] = new Matrix(0, 0);
                Tt[i][k] = new Matrix(0, 0);
            }
        }

        Map<StatefulNode, Matrix> cur_state = new HashMap<StatefulNode, Matrix>(sn.nstations);
        for (int i = 0; i < sn.nnodes; i++) {
            Node node_i = this.model.getNodes().get(i);
            if (node_i.isStateful()) {
                cur_state.put((StatefulNode) node_i, sn.state.get(node_i).copy());
            }
        }

        // Initialize the state matrices
        Matrix s0_sz = new Matrix(1, sn.state.size()); // Number of possible states for each station
        Matrix s0_id = s0_sz.copy(); // Used to iterate over all possible initial states

        // s0_sz is read below at the STATEFUL index (s0_id.get(isf)), so it has to
        // be filled in stateful order, as solver_fluid_analyzer.m does with
        // cellfun over sn.space. Filling it in STATION order agreed with that
        // only while every stateful node was a station; the MMT transformation
        // puts a Router where every fork was, and the two orders then diverge.
        int i = 0;
        for (StatefulNode statefulNode : this.model.getStatefulNodes()) {
            // Use sn.space for state count (matching MATLAB sn.space), fall back to sn.state
            Matrix spaceMatrix = (sn.space != null) ? sn.space.get(statefulNode) : null;
            if (spaceMatrix != null && spaceMatrix.getNumRows() > 0) {
                s0_sz.set(0, i, spaceMatrix.getNumRows());
            } else {
                Matrix stateMatrix = sn.state.get(statefulNode);
                s0_sz.set(0, i, stateMatrix == null ? 1 : stateMatrix.getNumRows());
            }
            i++;
        }
        Matrix s0_sz_1 = s0_sz.copy();
        s0_sz_1.addEq(-1);

        s0_id = PopulationLattice.pprod(s0_sz_1); // Generates all possible combinations of initial states for the network

        while (s0_id.elementMin() >= 0) { // For all possible initial states
            double s0prior_val = 1;
            for (int ind = 0; ind < sn.nnodes; ind++) { // iterate over all nodes
                if (sn.isstateful.get(ind) == 1) { // check if node is stateful

                    int isf = (int) sn.nodeToStateful.get(ind); // get stateful index of the node
                    int istate = (int) sn.nodeToStation.get(ind);
                    if (istate < 0) {
                        // A stateful node that is not a station: a Router, which
                        // the MMT transformation puts where every fork was. It
                        // holds no jobs, so it contributes nothing to the initial
                        // condition and has no station row to write back to.
                        // solver_fluid_initsol.m skips it the same way.
                        continue;
                    }

                    // Update prior
                    s0prior_val *= sn.stateprior.get(this.model.getStatefulNodes().get(isf)).get((int) (s0_id.get(isf)));

                    // Extract row from the state matrix corresponding to the current state
                    Matrix newState =
                            Matrix.extractRows(
                                    sn.space.get(this.model.getStatefulNodes().get(isf)),
                                    (int) s0_id.get(isf),
                                    (int) s0_id.get(isf) + 1,
                                    null);

                    // Update the state of the node
                    this.model.getStations().get(istate).setState(newState);
                }
            }

            // Update sn after updating the state of the stations (use true to refresh state)
            NetworkStruct sn_cur = this.model.getStruct(true).copy();

            // see _kb/06-solver-catalog.md (JAR-only implementation notes: non-Markovian to phase-type conversion before initSol)
            // the fluid ODEs read mu*phi as a flow, so the surrogate must be a genuine phase-type
            String phfitIter0 = options.config.phfit;
            options.config.phfit = "ph";
            sn_cur = SnNonmarkovToPh.snNonmarkovToPh(sn_cur, options, false);
            options.config.phfit = phfitIter0;

            //System.out.println("Prior probability (s0prior_val): " + s0prior_val);

            if (s0prior_val > 0) {
                // Clear init_sol so initSol() is called fresh for each state iteration
                options.init_sol = new Matrix(1, 0);
                // A non-hyperbolic fluid fixed point (balanced bottlenecks, a
                // saturated multiclass station, an overloaded open station) leaves
                // the linear noise approximation with no stationary covariance. It
                // cannot be seen before the mean is solved, so
                // FluidMinNormalApplicable cannot decline it and MinNormalAnalyzer
                // throws at the Lyapunov step.
                //
                // THE LADDER HAS TWO RUNGS, AND THE FIRST ONE KEEPS THE CLOSURE.
                // Most of these failures are not a property of the model at all:
                // MinNormalAnalyzer must start its alternation at sigma2 = 0, where
                // min(n,c) has no derivative, so a saturated or balanced model's
                // first-order fixed point lands on the kink, sits on a continuum of
                // equilibria, and the Jacobian there is neutral. DaeAnalyzer seeds
                // the variance POSITIVE and never adopts sigma2 = 0 as an iterate,
                // so the smoothed E[min(X,c)] breaks the degeneracy and the fixed
                // point is isolated and hyperbolic -- it answers the same closure,
                // with a covariance, where the alternation cannot. Dropping straight
                // to first order instead is not merely a lost second moment: on a
                // balanced two-station PS cycle at N=10 it returns [9 1] against the
                // exact [5 5], because a first-order method has no reason to prefer
                // one point of the continuum over another.
                //
                // The second rung is the first-order method, taken when dae declines
                // the model in advance (FluidDaeApplicable) or throws on the same
                // exception, which is the genuinely non-hyperbolic case: an unstable
                // open station has no stationary distribution to approximate under
                // any closure. Fall back whether "minnormal" was RESOLVED from
                // "default" or REQUESTED outright: the closure has no stationary
                // covariance either way, so refusing an explicit request would only
                // deny the caller the mean that is still available. matrix is that
                // method, except under DPS where closing is the one that applies.
                try {
                    runMethodSpecificAnalyzer(sn_cur); // run analyzer
                } catch (FluidNonHyperbolicException e) {
                    if (!options.method.equals("minnormal")) {
                        throw e;
                    }
                    boolean solved = false;
                    String daeReason = FluidDaeApplicable.reasonToDecline(sn_cur, options);
                    if (daeReason == null) {
                        String requested = options.method;
                        options.method = "dae";
                        options.init_sol = new Matrix(1, 0);
                        try {
                            runMethodSpecificAnalyzer(sn_cur);
                            solved = true;
                            line_debug(options.verbose,
                                    "Fluid minnormal declined at the Lyapunov step (" + e.getMessage()
                                            + "); falling back to dae");
                        } catch (FluidNonHyperbolicException daeEx) {
                            options.method = requested;
                            line_debug(options.verbose,
                                    "Fluid dae also declined at the Lyapunov step (" + daeEx.getMessage() + ")");
                        }
                    } else {
                        line_debug(options.verbose,
                                "Fluid dae not applicable as a fallback (" + daeReason + ")");
                    }
                    if (!solved) {
                        options.method = hasDPS ? "closing" : "matrix";
                        line_debug(options.verbose,
                                "Fluid moment closure declined at the Lyapunov step (" + e.getMessage()
                                        + "); falling back to " + options.method);
                        options.init_sol = new Matrix(1, 0);
                        runMethodSpecificAnalyzer(sn_cur);
                    }
                }

                // Handles the results returned by the solver

                // Note: in LINE, only unique time-step values (and their associated metrics) are stored.
                // The time inefficiency in determining the unique values is, I believe, worse than the
                // space inefficiency in storing larger arrays than is necessary. For that reason I"ve not
                // transferred the "unique" functionality across to JLINE

                if (((FluidResult) this.result).odeStateVec.isEmpty()) { // If the solution has failed
                    for (int k = 0; k < K; k++) {
                        for (int j = 0; j < M; j++) {
                            Q.set(j, k, NaN);
                            U.set(j, k, NaN);
                            R.set(j, k, NaN);
                            T.set(j, k, NaN);
                        }
                        C.set(1, k, NaN);
                        X.set(1, k, NaN);
                    }

                    Matrix nanMatrix = new Matrix(1, 1);
                    nanMatrix.set(0, 0, NaN);
                    for (int ist = 0; ist < M; ist++) {
                        for (int r = 0; r < K; r++) {
                            Qt[ist][r] = nanMatrix.copy();
                            Ut[ist][r] = nanMatrix.copy();
                            Tt[ist][r] = nanMatrix.copy();
                        }
                    }
                } else {
                    Q = Q.add(s0prior_val, result.QN);
                    U = U.add(s0prior_val, result.UN);
                    R = R.add(s0prior_val, result.RN);
                    T = T.add(s0prior_val, result.TN);
                    C = C.add(s0prior_val, result.CN);
                    X = X.add(s0prior_val, result.XN);

                    // Check if this is the first time adding results or not
                    if (Qt[0][0].isEmpty()) {
                        // First state: create 2-column [value, time] matrices matching MATLAB format
                        // MATLAB: Qt{ist,r} = [Qfull_t{ist,r} * s0prior_val, t];
                        for (int ist = 0; ist < M; ist++) {
                            for (int r = 0; r < K; r++) {
                                int nTimePoints = result.QNt[ist][r].getNumRows();
                                Qt[ist][r] = new Matrix(nTimePoints, 2);
                                Ut[ist][r] = new Matrix(nTimePoints, 2);
                                Tt[ist][r] = new Matrix(nTimePoints, 2);
                                for (int ti = 0; ti < nTimePoints; ti++) {
                                    double tVal = result.t.get(ti, 0);
                                    Qt[ist][r].set(ti, 0, result.QNt[ist][r].get(ti, 0) * s0prior_val);
                                    Qt[ist][r].set(ti, 1, tVal);
                                    Ut[ist][r].set(ti, 0, result.UNt[ist][r].get(ti, 0) * s0prior_val);
                                    Ut[ist][r].set(ti, 1, tVal);
                                    Tt[ist][r].set(ti, 0, result.TNt[ist][r].get(ti, 0) * s0prior_val);
                                    Tt[ist][r].set(ti, 1, tVal);
                                }
                            }
                        }
                    } else {
                        for (int ist = 0; ist < M; ist++) {
                            for (int r = 0; r < K; r++) {
                                // Merge time vectors and interpolate data
                                Matrix oldTimes = Qt[ist][r].getColumn(1); // Get time column
                                Matrix newTimes = result.t.copy();
                                
                                // Create union of time points
                                Set<Double> timeUnion = new TreeSet<>();
                                for (int ti = 0; ti < oldTimes.getNumRows(); ti++) {
                                    timeUnion.add(oldTimes.get(ti, 0));
                                }
                                for (int ti = 0; ti < newTimes.getNumRows(); ti++) {
                                    timeUnion.add(newTimes.get(ti, 0));
                                }
                                
                                // Convert to arrays for interpolation
                                double[] tunion = new double[timeUnion.size()];
                                int idx = 0;
                                for (Double t : timeUnion) {
                                    tunion[idx++] = t;
                                }
                                
                                // Interpolate old data at union time points
                                double[] oldTimesArray = new double[oldTimes.getNumRows()];
                                double[] oldDataArray = new double[oldTimes.getNumRows()];
                                for (int ti = 0; ti < oldTimes.getNumRows(); ti++) {
                                    oldTimesArray[ti] = oldTimes.get(ti, 0);
                                    oldDataArray[ti] = Qt[ist][r].get(ti, 0); // Get data column
                                }
                                
                                // Interpolate new data at union time points
                                double[] newTimesArray = new double[newTimes.getNumRows()];
                                double[] newDataArrayQ = new double[newTimes.getNumRows()];
                                double[] newDataArrayU = new double[newTimes.getNumRows()];
                                double[] newDataArrayT = new double[newTimes.getNumRows()];
                                
                                for (int ti = 0; ti < newTimes.getNumRows(); ti++) {
                                    newTimesArray[ti] = newTimes.get(ti, 0);
                                    newDataArrayQ[ti] = result.QNt[ist][r].get(ti, 0);
                                    newDataArrayU[ti] = result.UNt[ist][r].get(ti, 0);
                                    newDataArrayT[ti] = result.TNt[ist][r].get(ti, 0);
                                }
                                
                                // Prepare old utilization data for interpolation
                                double[] oldDataArrayU = new double[oldTimes.getNumRows()];
                                for (int ti = 0; ti < oldTimes.getNumRows(); ti++) {
                                    oldDataArrayU[ti] = Ut[ist][r].get(ti, 0);
                                }
                                
                                // Use linear interpolation
                                org.apache.commons.math3.analysis.interpolation.LinearInterpolator interpolator = 
                                    new org.apache.commons.math3.analysis.interpolation.LinearInterpolator();
                                
                                // Create interpolation functions
                                org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction oldInterpQ = 
                                    interpolator.interpolate(oldTimesArray, oldDataArray);
                                org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction oldInterpU = 
                                    interpolator.interpolate(oldTimesArray, oldDataArrayU);
                                org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction newInterpQ = 
                                    interpolator.interpolate(newTimesArray, newDataArrayQ);
                                org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction newInterpU = 
                                    interpolator.interpolate(newTimesArray, newDataArrayU);
                                org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction newInterpT = 
                                    interpolator.interpolate(newTimesArray, newDataArrayT);
                                
                                // Create new result matrices with interpolated values
                                Matrix newQt = new Matrix(tunion.length, 2);
                                Matrix newUt = new Matrix(tunion.length, 2);
                                Matrix newTt = new Matrix(tunion.length, 2);
                                
                                for (int ti = 0; ti < tunion.length; ti++) {
                                    double t = tunion[ti];
                                    
                                    // Queue length
                                    double oldValQ = (t >= oldTimesArray[0] && t <= oldTimesArray[oldTimesArray.length-1]) 
                                        ? oldInterpQ.value(t) : 0.0;
                                    double newValQ = (t >= newTimesArray[0] && t <= newTimesArray[newTimesArray.length-1]) 
                                        ? newInterpQ.value(t) : 0.0;
                                    newQt.set(ti, 0, oldValQ + s0prior_val * newValQ);
                                    newQt.set(ti, 1, t);
                                    
                                    // Utilization
                                    double oldValU = (t >= oldTimesArray[0] && t <= oldTimesArray[oldTimesArray.length-1]) 
                                        ? oldInterpU.value(t) : 0.0;
                                    double newValU = (t >= newTimesArray[0] && t <= newTimesArray[newTimesArray.length-1]) 
                                        ? newInterpU.value(t) : 0.0;
                                    newUt.set(ti, 0, oldValU + s0prior_val * newValU);
                                    newUt.set(ti, 1, t);
                                    
                                    // Throughput (commented out in MATLAB code)
                                    // double oldValT = ...
                                    // double newValT = ...
                                    // newTt.set(i, 0, oldValT + s0prior_val * newValT);
                                    // newTt.set(i, 1, t);
                                }
                                
                                // Update the result matrices
                                Qt[ist][r] = newQt;
                                Ut[ist][r] = newUt;
                                // Tt[ist][r] = newTt; // Commented out as in MATLAB
                            }
                        }
                    }
                }
            }

            s0_id = PopulationLattice.pprod(s0_id, s0_sz_1);
        }

        for (int j = 0; j < sn.nnodes; j++) {
            Node node_j = this.model.getNodes().get(j);
            if (node_j.isStateful()) {
                ((StatefulNode) node_j).setState(cur_state.get(node_j));
            }
        }
        result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        result.QN = Q;
        result.UN = U;
        result.RN = R;
        result.TN = T;
        result.CN = C;
        result.XN = X;
        // Strip time column from Qt/Ut/Tt before storing (internally they are [value, time]
        // but result.QNt/UNt/TNt should be value-only 1-column matrices)
        Matrix[][] QtOut = new Matrix[M][K];
        Matrix[][] UtOut = new Matrix[M][K];
        Matrix[][] TtOut = new Matrix[M][K];
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                if (Qt[ist][r].getNumCols() >= 2) {
                    QtOut[ist][r] = Qt[ist][r].getColumn(0);
                    UtOut[ist][r] = Ut[ist][r].getColumn(0);
                    TtOut[ist][r] = Tt[ist][r].getColumn(0);
                } else {
                    QtOut[ist][r] = Qt[ist][r];
                    UtOut[ist][r] = Ut[ist][r];
                    TtOut[ist][r] = Tt[ist][r];
                }
            }
        }
        // Set transient results using setTranAvgResults
        Matrix[][] RNt = new Matrix[0][0];
        Matrix[][] CNt = new Matrix[0][0];
        Matrix[][] XNt = new Matrix[0][0];
        this.setTranAvgResults(QtOut, UtOut, RNt, TtOut, CNt, XNt, result.runtime);
        //sn = this.model.getStruct();
        AvgHandle TH = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(sn, result.TN, TH);
        // Handle the "default" method case - add "default/" prefix if original method was "default"
        String finalMethod;
        if (origMethod.equals("default") && !options.method.equals("default")) {
            finalMethod = "default/" + options.method;
        } else {
            finalMethod = options.method;
        }
        result.method = finalMethod;
        // WN (waiting times) is not computed by the Fluid solver (matching MATLAB which passes [])
        this.setAvgResults(result.QN, result.UN, result.RN, result.TN, AN, new Matrix(0, 0), result.CN, result.XN, result.runtime, finalMethod, result.iter);
        // The resolution of "default" is PER SOLVE, not a property of the solver:
        // the non-hyperbolic fallback above fires only for a RESOLVED "default",
        // so leaving the resolved name here would make a second solve request it
        // explicitly. Matches @SolverFLD/runAnalyzer.m, which persists origMethod.
        options.method = origMethod;
    }

    private SolverResult runMethodSpecificAnalyzer(NetworkStruct sn) {

        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix phases = sn.phases.copy();
        Matrix phasesLast = sn.phases.copy();
        Matrix rates0 = sn.rates.copy();

        Matrix SCV = sn.scv.copy();
        Matrix gamma = new Matrix(M, 1);
        Matrix S = sn.nservers.copy();

        // sn.visits is indexed by STATEFUL node, not by station, so the chain
        // visits have to be read at each station's stateful row. Adding the
        // matrices whole worked only while the two index spaces coincided, i.e.
        // on a model with no stateful non-station node; the MMT transformation
        // puts a Router where every fork was, and the shapes then disagree.
        Matrix V = new Matrix(M, K);
        for (int c = 0; c < sn.visits.size(); c++) {
            Matrix Vc = sn.visits.get(c);
            for (int i = 0; i < M; i++) {
                int isf = (int) sn.stationToStateful.get(i);
                if (isf < 0 || isf >= Vc.getNumRows()) {
                    continue;
                }
                for (int k = 0; k < K; k++) {
                    V.set(i, k, V.get(i, k) + Vc.get(isf, k));
                }
            }
        }

        FluidAnalyzer analyzer;
        switch (options.method) {
            case "statedep":
            case "closing":
            case "softmin":
                line_debug(options.verbose, "Using ClosingAndStateDepMethodsAnalyzer for fluid analysis");
                analyzer = new ClosingAndStateDepMethodsAnalyzer();
                break;
            case "tbi":
                line_debug(options.verbose, "Using TbiAnalyzer for fluid analysis");
                analyzer = new TbiAnalyzer();
                break;
            case "matrix":
            case "pnorm":
                line_debug(options.verbose, "Using MatrixMethodAnalyzer for fluid analysis");
                analyzer = new MatrixMethodAnalyzer();
                break;
            case "mfq":
                line_debug(options.verbose, "Using MFQAnalyzer for fluid analysis");
                analyzer = new MFQAnalyzer();
                break;
            case "minnormal":
            // 'refined' is the SAME closure plus the O(1/N) correction of Gast
            // (POMACS 2017), taken about the mean-field fixed point; the
            // analyzer branches on the method name, as solver_fluid_moments.m
            // does.
            case "refined":
                line_debug(options.verbose, "Using MinNormalAnalyzer for fluid analysis");
                analyzer = new MinNormalAnalyzer();
                break;
            case "diffusion":
                line_debug(options.verbose, "Using DiffusionAnalyzer for fluid analysis");
                analyzer = new DiffusionAnalyzer();
                break;
            case "dae":
                line_debug(options.verbose, "Using DaeAnalyzer for fluid analysis");
                analyzer = new DaeAnalyzer();
                break;
            case "kp":
                line_debug(options.verbose, "Using KoPenderAnalyzer for fluid analysis");
                analyzer = new KoPenderAnalyzer();
                break;
            default:
                // Use default method as fallback
                line_debug(options.verbose, "Using default method for fluid analysis");
                analyzer = new ClosingAndStateDepMethodsAnalyzer();
                break;
        }

        if (options.init_sol.isEmpty()) {
            initSol();
        }
        analyzer.analyze(sn, options.copy(), result);

        if (analyzer instanceof MinNormalAnalyzer) {
            MinNormalAnalyzer mn = (MinNormalAnalyzer) analyzer;
            FluidResult fr = (FluidResult) this.result;
            fr.momentSigma = mn.sigmaMatrix;
            fr.momentQVar = mn.qVar;
            fr.momentClassBlock = mn.classBlock;
            fr.momentStationBlock = mn.stationBlock;
            fr.momentSigma2 = sigma2Matrix(mn.sigma2);
            fr.momentSigma2Drift = sigma2Matrix(mn.sigma2Drift);
            fr.momentOuterIters = mn.outerIters;
        }

        if (analyzer instanceof DaeAnalyzer) {
            DaeAnalyzer da = (DaeAnalyzer) analyzer;
            FluidResult fr = (FluidResult) this.result;
            fr.momentSigma = da.sigmaMatrix;
            fr.momentQVar = da.qVar;
            fr.momentClassBlock = da.classBlock;
            fr.momentStationBlock = da.stationBlock;
            fr.momentSigma2 = sigma2Matrix(da.sigma2);
            fr.momentSigma2Drift = sigma2Matrix(da.sigma2Drift);
            fr.momentOuterIters = da.outerIters;
            fr.daeResidual = da.residual;
            fr.daeConverged = da.converged;
            fr.daeConservation = da.conservation;
            fr.daeCapacityLabel = da.capacityLabel;
            fr.daeCapacityB = da.capacityB;
            fr.daeCapacityValue = da.capacityValue;
            fr.daeCapacityActive = da.capacityActive;
            fr.daeStaging = da.staging;
            fr.daeStagingRegion = da.stagingRegion;
            fr.daeStagingClass = da.stagingClass;
            fr.daeBlocked = da.blocked;
            fr.daeDrain = da.drain;
            fr.daeCapacityStaged = da.capacityStaged;
            fr.daeCapacityRegion = da.capacityRegion;
            fr.daeCapacityStation = da.capacityStation;
            fr.daeCapacitySwitches = da.capacitySwitches;
            // The transient second moment, on the SAME grid the mean was
            // reported on, so the two need no interpolation onto one another.
            fr.Sigmat = da.sigmat;
            if (da.qVart != null) {
                int mm = sn.nstations;
                int kk = sn.nclasses;
                fr.QVart = new Matrix[mm][kk];
                for (int i = 0; i < mm; i++) {
                    for (int k = 0; k < kk; k++) {
                        fr.QVart[i][k] = da.qVart[i * kk + k];
                    }
                }
            }
        }

        // MFQ is handled by early return above; this is a safety check
        if (Objects.equals(options.method, "mfq")) {
            ((FluidResult) this.result).odeStateVec = analyzer.getXVecIt();
            ((FluidResult) this.result).snFinal = this.sn;
            return this.result;
        }

        // Single iteration is sufficient for statedep, so do nothing unless matrix,
        // closing or the moment closure, which shares the closing FCFS treatment
        if ((Objects.equals(options.method, "matrix")) || (Objects.equals(options.method, "closing"))
                || (Objects.equals(options.method, "minnormal"))) {
            if (sn.sched.containsValue(SchedStrategy.FCFS)) {
                int iter = 0;
                Matrix eta_1 = new Matrix(M, 1);
                Matrix eta = new Matrix(M, 1);
                for (int i = 0; i < M; i++) {
                    eta.set(i, 0, POSITIVE_INFINITY);
                }
                double tol = 0.01;

                while (iter <= options.iter_max) {
                    double maxAbsEtaEta_1 = NEGATIVE_INFINITY;
                    for (int i = 0; i < M; i++) {
                        double val = FastMath.abs(1 - (eta.get(i, 0) / eta_1.get(i, 0)));
                        if (val > maxAbsEtaEta_1) {
                            maxAbsEtaEta_1 = val;
                        }
                    }
                    if (maxAbsEtaEta_1 <= tol) {
                        break;
                    }

                    iter++;
                    eta_1 = eta.copy();
                    for (int i = 0; i < M; i++) {
                        LinkedList<Integer> sdCols = new LinkedList<>();
                        for (int k = 0; k < K; k++) {
                            if (rates0.get(i, k) > 0) {
                                sdCols.add(k);
                            }
                        }
                        for (int k : sdCols) {
                            result.UN.set(i, k, result.TN.get(i, k) / rates0.get(i, k));
                        }
                    }

                    Matrix ST0 = new Matrix(M, K);
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j < K; j++) {
                            ST0.set(i, j, 1.0 / rates0.get(i, j));
                            if (isInfinite(ST0.get(i, j))) {
                                ST0.set(i, j, 1 / GlobalConstants.Zero);
                            } else if (isNaN(ST0.get(i, j))) {
                                ST0.set(i, j, GlobalConstants.Zero);
                            }
                        }
                    }

                    Matrix XN = new Matrix(1, K);
                    for (int k = 0; k < K; k++) {
                        if (sn.refstat.get(k, 0) >= 0) { // Ignore artificial classes
                            XN.set(0, k, result.TN.get((int) sn.refstat.get(k, 0), k));
                        }
                    }

                    Ret.npfqnNonexpApprox NPFQNret = npfqn_nonexp_approx(options.config.highvar == null ? "interp" : options.config.highvar, sn, ST0, V, SCV, result.TN, result.UN, gamma, S);
                    Matrix ST = NPFQNret.ST;
                    eta = new Matrix(NPFQNret.eta);

                    Matrix rates = new Matrix(M, K);
                    for (int i = 0; i < M; i++) {
                        for (int j = 0; j < K; j++) {
                            rates.set(i, j, 1.0 / ST.get(i, j));
                            if (isInfinite(rates.get(i, j))) {
                                rates.set(i, j, 1 / GlobalConstants.Zero);
                            } else if (isNaN(rates.get(i, j))) {
                                rates.set(i, j, GlobalConstants.Zero);
                            }
                        }
                    }

                    for (int i = 0; i < M; i++) {
                        if (sn.sched.get(this.model.getStations().get(i)) == SchedStrategy.FCFS) {
                            for (int k = 0; k < K; k++) {
                                if (rates.get(i, k) > 0 && sn.scv.get(i, k) > 0) {
                                    Coxian cx = Coxian.fitMeanAndSCV(1.0 / rates.get(i, k), sn.scv.get(i, k));
                                    Matrix muik = cx.getMu();
                                    Matrix phiik = cx.getPhi();
                                    // we now handle the case that due to either numerical issues or different
                                    // relationship between scv and mean if the size of the phase-type representation
                                    // has changed
                                    phases.set(i, k, muik.length());
                                    State.StateMarginalStatistics stats = null;
                                    if (phases.get(i, k) != phasesLast.get(i, k)) { // If number of phases changed
                                        // Before we update sn we adjust the initial state
                                        int isf = (int) sn.stationToStateful.get(0, i);
                                        stats =
                                                ToMarginal.toMarginal(
                                                        sn,
                                                        i,
                                                        sn.state.get(this.model.getStatefulNodes().get(isf)),
                                                        null,
                                                        null,
                                                        null,
                                                        null,
                                                        null);
                                    }
                                    sn.proc.get(this.model.getStations().get(i)).put(sn.jobclasses.get(k), cx.getProcess());
                                    sn.mu.get(this.model.getStations().get(i)).put(sn.jobclasses.get(k), muik);
                                    sn.phi.get(this.model.getStations().get(i)).put(sn.jobclasses.get(k), phiik);
                                    sn.phases = phases.copy();

                                    sn.phasessz = sn.phases.copy();
                                    for (int row = 0; row < M; row++) {
                                        for (int col = 0; col < K; col++) {
                                            if (sn.phasessz.get(row, col) < 1) {
                                                sn.phasessz.set(row, col, 1);
                                            }
                                        }
                                    }

                                    sn.phaseshift = new Matrix(sn.phases.getNumRows(), 1);
                                    sn.phaseshift =
                                            Matrix.concatColumns(sn.phaseshift, sn.phasessz.cumsumViaRow(), null);

                                    if (phases.get(i, k) != phasesLast.get(i, k)) {
                                        int isf = (int) sn.stationToStateful.get(0, i);
                                        // We now initialise the new service process
                                        assert stats != null;
                                        sn.state.put(
                                                this.model.getStatefulNodes().get(isf),
                                                FromMarginal.fromMarginalAndStarted(sn, i, stats.nir, stats.sir));
                                        // Pick one as the marginals won"t change
                                        sn.state.put(
                                                this.model.getStatefulNodes().get(isf),
                                                Matrix.extractRows(sn.state.get(this.model.getStatefulNodes().get(isf)), 0, 1, null));
                                    }
                                }
                            }
                        }

                        options.init_sol = analyzer.getXVecIt();
                        options.init_sol.transpose();
                        // If there is a change of phases reset
                        boolean callInitSol = false;
                        for (int row = 0; row < M; row++) {
                            for (int col = 0; col < K; col++) {
                                if (phasesLast.get(row, col) - phases.get(row, col) != 0) {
                                    callInitSol = true;
                                }
                            }
                        }
                        if (callInitSol) {
                            initSol();
                        }
                    }
                    sn.phases = phases.copy();
                    analyzer.analyze(sn, options.copy(), result);
                    phasesLast = phases.copy();
                    result.iter = iter;
                } // FCFS iteration ends here

                // The FCFS iteration reinitializes at the solution of the last iterative step. We now
                // have converged in the substitution of the model parameters and we rerun everything from
                // the true initial point so that we get the correct transient.
                initSol();
                analyzer.analyze(sn, options.copy(), result);
            }
        }

        if (result.t.value() == 0) {
            result.t.set(0, 0, 0.00000001);
        }

        // A CLASS THE MODEL NEVER ROUTES HERE HAS NO RESPONSE TIME, and QN
        // alone does not say so: it holds a remnant of the initial state the
        // integrator was still draining when it stopped. See
        // fluidVisitedPairs -- the visit ratios decide it, a threshold on QN
        // or TN cannot.
        boolean[][] visited = fluidVisitedPairs(sn, M, K);
        Matrix Ufull0 = result.UN.copy();
        for (int i = 0; i < M; i++) {
            List<Integer> sdCols = new LinkedList<>();
            for (int k = 0; k < K; k++) {
                if (result.QN.get(i, k) > 0 && visited[i][k]) {
                    sdCols.add(k);
                } else {
                    result.UN.set(i, k, 0);
                    result.RN.set(i, k, 0);
                }
            }

            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                for (int k : sdCols) {
                    result.UN.set(i, k, result.QN.get(i, k));
                    result.UNt[i][k] = result.QNt[i][k].copy();
                    result.UNt[i][k].scaleEq(sn.rates.get(i, k), result.TNt[i][k]);
                }
            } else {
                double sumUfull0sd = 0;
                double sumTNDivRates0sd = 0;
                for (int k : sdCols) {
                    sumUfull0sd += Ufull0.get(i, k);
                    sumTNDivRates0sd += result.TN.get(i, k) / rates0.get(i, k);
                }
                for (int k : sdCols) {
                    // correct for the real rates, instead of the diffusion approximation rates
                    result.UN.set(
                            i,
                            k,
                            min(
                                    1,
                                    min(
                                            result.QN.get(i, k) / sn.nservers.get(i, 0),
                                            sumUfull0sd * result.TN.get(i, k) / rates0.get(i, k) / sumTNDivRates0sd)));
                    result.UNt[i][k].scaleEq(sn.rates.get(i, k) * sn.nservers.get(i, 0), result.TNt[i][k]);
                    // NO DEPARTURES MEANS NO RESIDENCE TIME TO READ, and a QN the
                    // integrator has not finished draining must not be divided by a TN of
                    // the same origin: the ratio is O(1) and the residual QN is rounded
                    // away later, leaving a queue length of 0 beside a response time of 10.
                    // The pairs the model never routes to are out of sdCols already; this
                    // guards what is left. MinNormalAnalyzer guards its own division too.
                    if (result.TN.get(i, k) > GlobalConstants.Zero) {
                        result.RN.set(i, k, (result.QN.get(i, k) / result.TN.get(i, k)));
                    } else {
                        result.RN.set(i, k, 0);
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (Double.isNaN(result.UN.get(i, j))) {
                    result.UN.set(i, j, 0);
                }
                if (Double.isNaN(result.RN.get(i, j))) {
                    result.RN.set(i, j, 0);
                }
            }
        }

        result.XN = new Matrix(1, K);
        result.CN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            if (sn.refstat.get(k, 0) >= 0) { // Ignore artificial classes
                result.XN.set(0, k, result.TN.get((int) sn.refstat.get(k, 0), k));
                result.CN.set(0, k, sn.njobs.get(0, k) / result.XN.get(0, k));
            }
        }

        ((FluidResult) this.result).odeStateVec = analyzer.getXVecIt();
        
        // IMPORTANT: Store the modified network structure for use in getCdfRespT
        // The phases may have been expanded during fluid analysis
        ((FluidResult) this.result).snFinal = this.sn;

        // Calculate arrival rates from throughputs using the routing matrix
        AvgHandle TH = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(this.sn, result.TN, TH);

        // WN (waiting times) is not computed by the Fluid solver (matching MATLAB which passes [])
        Matrix WN = new Matrix(0, 0);
        this.setAvgResults(result.QN, result.UN, result.RN, result.TN, AN, WN, result.CN, result.XN, 0.0, options.method, 0);
        // Run AoI analysis if this is a valid AoI topology
        runAoiAnalysis(sn);
        return this.result;
    }

    public SolverResult runMethodSpecificAnalyzer() {
        return runMethodSpecificAnalyzer(this.sn);
    }

    /**
     * Solves a fork-join model through the solver-agnostic fixed point.
     *
     * <p>The MMT transformation turns every fork into a router, every join into
     * a zero-service delay, and carries the parallelism on auxiliary open
     * classes; the resulting network is built from Source, Delay, Queue, Router
     * and ClassSwitch alone, every one of which the fluid drift already carries.
     * Nothing in the loop is fluid-specific, so the driver is the one MVA and NC
     * use ({@link jline.solvers.fj.FJFixedPoint}).</p>
     *
     * <p>Only the steady-state means are produced. The transient tables are
     * indexed by the ORIGINAL stations and classes, whereas each pass integrates
     * a different transformed network, so getTranAvg stays unavailable on a
     * fork-join model.</p>
     *
     * @param startTime the wall-clock marker of the enclosing analyzer
     */
    private void runForkJoinAnalyzer(long startTime) {
        jline.solvers.fj.FJFixedPoint.FJState fjState =
                new jline.solvers.fj.FJFixedPoint.FJState(null, null);
        jline.solvers.fj.FJFixedPoint.FJOutcome fjOut = jline.solvers.fj.FJFixedPoint.run(
                this.model, this.model.getStruct(true), this.options, fjState,
                new jline.solvers.fj.FJFixedPoint.InnerSolve() {
                    @Override
                    public jline.solvers.mva.MVAResult solve(Network net, NetworkStruct snIn,
                                                             SolverOptions opts) {
                        return SolverFluid.this.fldDispatch(net, snIn, opts);
                    }
                }, startTime);
        jline.solvers.mva.MVAResult fjRet = fjOut.ret;
        // The driver leaves the transformed struct behind; the tail below indexes
        // the ORIGINAL stations and classes, so recompile it.
        this.sn = this.model.getStruct(true);
        // The transform appends its own Source to carry the auxiliary open
        // classes, so the metrics come back with one row per TRANSFORMED station.
        // The original stations are the prefix of that list; trim the rest.
        int Mfj = this.sn.nstations;
        fjRet.QN = trimToStations(fjRet.QN, Mfj);
        fjRet.UN = trimToStations(fjRet.UN, Mfj);
        fjRet.RN = trimToStations(fjRet.RN, Mfj);
        fjRet.TN = trimToStations(fjRet.TN, Mfj);
        AvgHandle TH = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(this.sn, fjRet.TN, TH);
        double runtimeFj = (System.nanoTime() - startTime) / 1000000000.0;
        // Only setAvgResults writes the result matrices: it trims the transformed
        // model's extra rows (its own Source) back to the original stations, so
        // assigning result.QN directly here would leak that row to the caller.
        ((FluidResult) this.result).snFinal = this.sn;
        this.setAvgResults(fjRet.QN, fjRet.UN, fjRet.RN, fjRet.TN, AN, new Matrix(0, 0),
                fjRet.CN, fjRet.XN, runtimeFj, options.method, fjOut.iter);
    }

    /** Keeps the first {@code M} station rows of a metric matrix. */
    private static Matrix trimToStations(Matrix metric, int M) {
        if (metric == null || metric.getNumRows() <= M) {
            return metric;
        }
        return Matrix.extractRows(metric, 0, M, null);
    }

    /**
     * One inner solve of the fluid analyzer, in the contract
     * {@link jline.solvers.fj.FJFixedPoint} expects.
     *
     * <p>The transformed model is solved by a fresh SolverFluid built on it
     * rather than by reusing this one: initSol keys sn.state by the node objects
     * of the model the struct came from, so this solver, which holds the
     * ORIGINAL model, cannot integrate the transformed struct. The transformed
     * model has no fork, so this call takes the ordinary path and never
     * re-enters the fixed point.</p>
     *
     * @param net the transformed network
     * @param snIn its struct
     * @param opts the solver options
     * @return the metrics of that solve
     */
    private jline.solvers.mva.MVAResult fldDispatch(Network net, NetworkStruct snIn,
                                                    SolverOptions opts) {
        long t0 = System.nanoTime();
        SolverOptions o = opts.copy();
        o.init_sol = new Matrix(1, 0);
        // The fluid ODE integrates from an initial condition, so the transformed
        // model must carry a consistent state, stateprior and space triple. The
        // transform builds the model but not its state, and the driver rebuilds
        // the struct on every pass, so seed it here rather than in the
        // solver-agnostic driver, which MVA and NC share and which never reads
        // the state at all.
        net.initDefault();
        SolverFluid inner = new SolverFluid(net, o);
        inner.runAnalyzer();
        jline.solvers.mva.MVAResult out = new jline.solvers.mva.MVAResult();
        out.QN = inner.result.QN;
        out.UN = inner.result.UN;
        out.RN = inner.result.RN;
        out.TN = inner.result.TN;
        out.CN = inner.result.CN;
        out.XN = inner.result.XN;
        // The fluid solver produces no normalizing constant.
        out.logNormConstAggr = Double.NaN;
        out.runtime = (System.nanoTime() - t0) / 1000000000.0;
        out.iter = inner.result.iter;
        out.method = o.method;
        return out;
    }

    /**
     * Checks whether the given model is supported by the Fluid solver.
     * This method compares the features used by the model against the
     * features supported by the Fluid solver.
     *
     * @param model The network model to check
     * @return true if the model is supported, false otherwise
     */
    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverFluid.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Get model structure data structure summarizing the network
     * 
     * @return NetworkStruct containing the model structure
     */
    public NetworkStruct getStruct() {
        return this.model.getStruct(true);
    }

    /**
     * Export the mean-field ODE system to LaTeX in scalar notation.
     *
     * @return LaTeX source of the exported ODE system
     */
    public String exportODEs() {
        return exportODEs("", "scalar");
    }

    /**
     * Export the mean-field ODE system to a LaTeX file in scalar notation.
     *
     * @param fileName path of the .tex file to write; empty returns the source only
     * @return LaTeX source of the exported ODE system
     */
    public String exportODEs(String fileName) {
        return exportODEs(fileName, "scalar");
    }

    /**
     * Export the system of ODEs integrated by the mean-field methods of this
     * solver (default/matrix, pnorm, closing, statedep, softmin) as a
     * standalone LaTeX document, in a symbolic form that is both human and
     * machine readable. Mirrors the MATLAB SolverFLD.exportODEs method.
     *
     * @param fileName path of the .tex file to write; empty or null returns the source only
     * @param notation "scalar" for one expanded ODE per state variable, or
     *                 "matrix" for the compact matrix notation
     * @return LaTeX source of the exported ODE system
     */
    public String exportODEs(String fileName, String notation) {
        this.sn = this.model.getStruct(true);
        if (this.options.init_sol == null || this.options.init_sol.isEmpty()) {
            initSol();
        }
        FluidODEsExporter.SymODEs sys = FluidODEsExporter.build(this.sn, this.options);
        String tex = FluidODEsExporter.render(sys, this.options, this.model.getName(), notation);
        if (fileName != null && !fileName.isEmpty()) {
            try {
                java.io.FileWriter writer = new java.io.FileWriter(fileName);
                writer.write(tex);
                writer.close();
            } catch (java.io.IOException e) {
                throw new RuntimeException("Cannot write ODE export to file '" + fileName + "': " + e.getMessage(), e);
            }
        }
        return tex;
    }

    /**
     * Right-hand side of the mean-field ODE system as expression strings, one
     * per state variable, in the format the symbolic backend parses.
     *
     * <p>Only smooth drifts are exported; see
     * {@link FluidODEsExporter#symbolicDrift}.</p>
     *
     * @return one expression per state variable
     */
    public java.util.List<String> getSymbolicDrift() {
        this.sn = this.model.getStruct(true);
        if (this.options.init_sol == null || this.options.init_sol.isEmpty()) {
            initSol();
        }
        return FluidODEsExporter.symbolicDrift(
                FluidODEsExporter.build(this.sn, this.options));
    }

    /**
     * Jacobian of the mean-field drift, d f_i / d x_j, computed exactly by the
     * computer algebra backend.
     *
     * <p>The Jacobian is what tells a fixed point apart from a limit cycle and
     * gives the local convergence rate of the fluid approximation, neither of
     * which a numerical integration reports. Only smooth drifts have one: the
     * min-scaled methods are refused by name rather than answered with a
     * one-sided derivative, see {@link FluidODEsExporter#symbolicDrift}.</p>
     *
     * @return the Jacobian entries, row major
     * @throws RuntimeException if no symbolic backend is available or the drift
     *                          is not differentiable
     */
    public String[][] getJacobian() {
        this.sn = this.model.getStruct(true);
        if (this.options.init_sol == null || this.options.init_sol.isEmpty()) {
            initSol();
        }
        FluidODEsExporter.SymODEs sys = FluidODEsExporter.build(this.sn, this.options);
        java.util.List<String> rhs = FluidODEsExporter.symbolicDrift(sys);
        java.util.List<String> vars = FluidODEsExporter.stateVariables(sys);

        SymEngine engine = SymEngines.resolve(this.options.config.symbolic);
        if (engine == null) {
            throw new RuntimeException(
                    "No symbolic backend is available. Start one with "
                            + "'docker run -d -p 8080:8080 " + SymEngines.DOCKER_IMAGE
                            + "', point " + SymEngines.URL_ENV + " at a running service, or "
                            + "set options.config.symbolic to its URL.");
        }
        if (engine instanceof SageRestEngine) {
            ((SageRestEngine) engine).setTimeoutSeconds(this.options.config.symbolic_timeout);
        }
        try {
            return engine.fluidODEs(rhs, vars,
                    java.util.Collections.singletonList("jacobian")).jacobian;
        } catch (java.io.IOException e) {
            throw new RuntimeException("Jacobian failed on backend " + engine.name() + ": "
                    + e.getMessage(), e);
        }
    }

    /**
     * List all valid solution methods supported by this solver
     *
     * @return array of valid method names
     */
    /**
     * Transient queue-length VARIANCE per station and class, [stations][classes].
     *
     * <p>Only the "kp" method computes a second moment: it integrates the covariance of the
     * Ko-Pender diffusion limit alongside the fluid mean. The full state covariance, which
     * keeps cross-station and cross-class terms, is on {@code result.Sigmat}.
     *
     * @throws RuntimeException if the solver is not configured with options.method = "kp"
     */
    public Matrix[][] getTranAvgVar() {
        if (!"kp".equals(options.method)) {
            line_error(mfilename(new Object() {
            }), "getTranAvgVar needs options.method = \"kp\"; the other fluid methods "
                    + "integrate the mean only and carry no second moment.");
        }
        if (!(result instanceof FluidResult) || ((FluidResult) result).QVart == null) {
            getTranAvg();
        }
        return ((FluidResult) result).QVart;
    }

    public String[] listValidMethods() {
        // Every method the dispatch accepts, INCLUDING the 'fluid.'-qualified
        // spelling of each, which is the set MATLAB, native python and the C++
        // port all advertise. 'butools' (the MFQ backend) and 'aoi' (its
        // age-of-information reading) are aliases of 'mfq'. 'pnorm', 'diffusion'
        // and 'refined' were the three methods this port did not carry at all
        // and are now served by MatrixMethodODE's p-norm smoothing, by
        // DiffusionAnalyzer and by MinNormalAnalyzer's refined branch.
        return new String[]{"default",
                "matrix", "fluid.matrix", "pnorm", "fluid.pnorm",
                "softmin", "fluid.softmin",
                "statedep", "fluid.statedep",
                "closing", "fluid.closing",
                "minnormal", "fluid.minnormal",
                "refined", "fluid.refined",
                "tbi", "fluid.tbi",
                "diffusion", "fluid.diffusion",
                "mfq", "fluid.mfq", "butools",
                "rmf", "fluid.rmf",
                "aoi", "fluid.aoi",
                "kp", "fluid.kp",
                "dae", "fluid.dae",
                // The single-station fluid limits (Source -> Queue -> Sink, one
                // class). QsysLimitAnalyzer refuses any other shape by name.
                // "ggisgi" and "tga" are the SHORT spellings, mapped onto the
                // two primary names as the C++ fluid_qsys_canonical does
                "ggisgi.fluid", "fluid.ggisgi", "ggisgi",
                "ggingi.tga", "fluid.tga", "tga",
                "tvms", "fluid.tvms", "mtginf", "fluid.mtginf", "mol", "fluid.mol"};
    }

    /**
     * Validates model compatibility and method support before analysis
     * 
     * @param options solver options containing method specification
     * @throws RuntimeException if model contains unsupported features or method is invalid
     */
    /** Whether the model holds a Transition node, i.e. is a stochastic Petri net. */
    private static boolean isPetriNet(jline.lang.NetworkStruct sn) {
        if (sn == null || sn.nodetype == null) {
            return false;
        }
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == jline.lang.constant.NodeType.Transition) {
                return true;
            }
        }
        return false;
    }

    /**
     * The 'dae' method's Petri arm.
     *
     * <p>'dae' is the ONLY fluid method that can carry a net: the P-invariants, the firing flow of
     * an immediate transition and a bounded place are all EQUATIONS, and the other methods have
     * nowhere to put them. A method named explicitly is therefore checked rather than silently
     * redirected, so a caller who asked for 'closing' on a net is told why it cannot answer.
     */
    private void runPetriAnalyzer() {
        String m = fluidUnqualify(this.options.method);
        if (!("default".equals(m) || "dae".equals(m))) {
            throw new RuntimeException("SolverFLD method '" + this.options.method
                    + "' cannot solve a Petri net: its conserved quantities are P-invariants "
                    + "rather than chain populations, an immediate transition is an algebraic "
                    + "FLOW rather than an event with a rate, and a bounded place is a linear "
                    + "inequality on the marking. Only 'dae' states those as equations; every "
                    + "other fluid method builds its drift from the station/class/phase encoding, "
                    + "where a Place contributes no coordinate at all, and would integrate the "
                    + "net as an empty model and report zeros without a warning.");
        }
        jline.lang.NetworkStruct sn = this.model.getStruct(false);
        jline.solvers.fluid.petri.PetriSolver.Result pr =
                new jline.solvers.fluid.petri.PetriSolver(sn, this.options).solve();
        for (String w : pr.warnings) {
            line_warning(this.getName(), "%s", w);
        }
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix CN = new Matrix(1, K);
        Matrix XN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            double q = 0.0;
            double x = 0.0;
            for (int i = 0; i < M; i++) {
                q += pr.QN.get(i, k);
                x = Math.max(x, pr.TN.get(i, k));
            }
            XN.set(0, k, x);
            CN.set(0, k, (x > jline.GlobalConstants.Zero) ? q / x : 0.0);
        }
        this.setAvgResults(pr.QN, pr.UN, pr.RN, pr.TN, new Matrix(0, 0), new Matrix(0, 0),
                CN, XN, pr.runtime, "dae", pr.iters);
        FluidResult fr = (FluidResult) this.result;
        fr.momentSigma = pr.Sigma;
        fr.momentQVar = pr.QVar;
        fr.daeResidual = pr.resnorm;
        fr.daeConverged = pr.converged;
        fr.petri = pr.petri;
    }

    public void runAnalyzerChecks(SolverOptions options) {
        // Propagate solver verbose level to global
        if (options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        // THE COARSE GATE, and it is method-agnostic: supports() reads the STATIC
        // SolverFluid.getFeatureSet, not getMethodFeatureSet. The method-aware
        // gate is strictly finer -- it starts from the same static set and then
        // narrows it per method -- so the two can only disagree where a method
        // WIDENS the set, and 'dae' is the one that does: it declares Region,
        // which the static set must keep false because every other method has to
        // go on rejecting a finite capacity region. Deferring to the finer
        // verdict for that method is what lets the declaration stand; this check
        // runs unchanged for every other method.
        // The single-station fluid limits widen the set too: three of them
        // declare Reneging, which the static set must keep false because the
        // network drift carries no abandonment flow.
        boolean isDae = options != null
                && ("dae".equals(options.method) || "fluid.dae".equals(options.method)
                    || QsysLimitAnalyzer.handles(options.method));
        if (!isDae && !supports(this.model)) {
            throw new RuntimeException("This model contains features not supported by the Fluid solver.");
        }
        
        // Validate the selected method
        String[] validMethods = listValidMethods();
        boolean isValidMethod = false;
        for (String method : validMethods) {
            if (method.equals(options.method)) {
                isValidMethod = true;
                break;
            }
        }
        
        if (!isValidMethod) {
            throw new RuntimeException("The '" + options.method + "' method is unsupported by the Fluid solver. " +
                    "Valid methods are: " + String.join(", ", validMethods));
        }
        
        // Additional fluid-specific checks
        NetworkStruct sn = getStruct();
        
        // Check for unsupported scheduling strategies with certain methods
        if (options.method.equals("matrix")) {
            for (int i = 0; i < sn.nstations; i++) {
                if (sn.sched.get((int) sn.stationToNode.get(i)) == SchedStrategy.SJF ||
                    sn.sched.get((int) sn.stationToNode.get(i)) == SchedStrategy.LJF) {
                    throw new RuntimeException("Matrix method does not support SJF or LJF scheduling strategies.");
                }
            }
        }
        
        // Check for mixed open/closed classes with certain methods
        boolean hasOpenClasses = false;
        boolean hasClosedClasses = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (Double.isInfinite(sn.njobs.get(k))) {
                hasOpenClasses = true;
            } else {
                hasClosedClasses = true;
            }
        }
        
        if (hasOpenClasses && hasClosedClasses && options.method.equals("statedep")) {
            throw new RuntimeException("State-dependent method does not support mixed open/closed class models.");
        }
        
    }

    /**
     * Initialize performance metric handles from the model.
     * This method ensures all handles are properly configured for fluid analysis.
     */
    protected void initFluidHandles() {
        // Call parent initialization
        this.initHandles();
        
        // Ensure handles are refreshed from the model
        if (this.model != null) {
            this.avgHandles = this.model.getAvgHandles();
            this.tranHandles = this.model.getTranHandles();
        }
        
        // Validate that all required handles are available
        if (this.avgHandles == null) {
            throw new RuntimeException("Average handles not available from model");
        }
        
        // Initialize result structure if needed
        if (this.result == null) {
            this.result = new FluidResult();
        }
    }

    /**
     * Set transient average results with proper validation
     * 
     * @param Qt transient queue length matrices [time][stations x classes]
     * @param Ut transient utilization matrices [time][stations x classes]  
     * @param Rt transient response time matrices [time][stations x classes]
     * @param Tt transient throughput matrices [time][stations x classes]
     * @param Ct transient system response time matrices [time][chains]
     * @param Xt transient system throughput matrices [time][chains]
     * @param runtimet computation time for transient analysis
     */
    public void setFluidTranAvgResults(Matrix[][] Qt, Matrix[][] Ut, Matrix[][] Rt, Matrix[][] Tt, 
                                      Matrix[][] Ct, Matrix[][] Xt, double runtimet) {
        if (this.result == null) {
            this.result = new FluidResult();
        }
        
        // Store transient results using parent method
        super.setTranAvgResults(Qt, Ut, Rt, Tt, Ct, Xt, runtimet);
        
        // Update solver metadata
        this.result.solver = this.getName();
        this.result.method = this.getOptions().method;
        this.result.runtime = runtimet;
    }

    /**
     * Set distribution results with enhanced metadata for fluid solver
     * 
     * @param RD distribution data [stations x classes] containing CDFs
     * @param runtime computation time for distribution analysis
     */
    public void setFluidDistribResults(Matrix RD, double runtime) {
        if (this.result == null) {
            this.result = new FluidResult();
        }
        
        // Store in FluidResult format - note: distribC expects different type
        FluidResult fluidResult = (FluidResult) this.result;
        // distribC appears to be a different type in FluidResult, so we'll store runtime only
        fluidResult.distribRuntime = runtime;
        
        // Call parent method for standard distribution result storage
        super.setDistribResults(RD, runtime);
    }

    /**
     * Enhanced method to get average handles with validation
     * 
     * @return SolverAvgHandles with all handles properly initialized
     */
    @Override
    public SolverAvgHandles getAvgHandles() {
        if (this.avgHandles == null) {
            initFluidHandles();
        }
        return this.avgHandles;
    }

    /**
     * Enhanced method to set average handles with validation
     * 
     * @param handles the average handles to set
     */
    @Override
    public void setAvgHandles(SolverAvgHandles handles) {
        if (handles == null) {
            throw new IllegalArgumentException("Handles cannot be null");
        }
        this.avgHandles = handles;
    }

    /**
     * Enhanced method to get transient handles with validation
     * 
     * @return SolverTranHandles with all handles properly initialized
     */
    @Override
    public SolverTranHandles getTranHandles() {
        if (this.tranHandles == null) {
            initFluidHandles();
        }
        return this.tranHandles;
    }

    /**
     * Enhanced method to set transient handles with validation
     * 
     * @param handles the transient handles to set
     */
    @Override
    public void setTranHandles(SolverTranHandles handles) {
        if (handles == null) {
            throw new IllegalArgumentException("Handles cannot be null");
        }
        this.tranHandles = handles;
    }

    /**
     * Get queue length handles with fluid-specific validation
     * 
     * @return AvgHandle for queue length metrics
     */
    @Override
    public AvgHandle getAvgQLenHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.Q == null || handles.Q.isEmpty()) {
            throw new RuntimeException("Queue length handles not available");
        }
        return handles.Q;
    }

    /**
     * Get utilization handles with fluid-specific validation
     * 
     * @return AvgHandle for utilization metrics
     */
    @Override
    public AvgHandle getAvgUtilHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.U == null || handles.U.isEmpty()) {
            throw new RuntimeException("Utilization handles not available");
        }
        return handles.U;
    }

    /**
     * Get response time handles with fluid-specific validation
     * 
     * @return AvgHandle for response time metrics
     */
    @Override
    public AvgHandle getAvgRespTHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.R == null || handles.R.isEmpty()) {
            throw new RuntimeException("Response time handles not available");
        }
        return handles.R;
    }

    /**
     * Get throughput handles with fluid-specific validation
     * 
     * @return AvgHandle for throughput metrics
     */
    @Override
    public AvgHandle getAvgTputHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.T == null || handles.T.isEmpty()) {
            throw new RuntimeException("Throughput handles not available");
        }
        return handles.T;
    }

    /**
     * Get arrival rate handles with fluid-specific validation
     * 
     * @return AvgHandle for arrival rate metrics
     */
    @Override
    public AvgHandle getAvgArvRHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.A == null || handles.A.isEmpty()) {
            throw new RuntimeException("Arrival rate handles not available");
        }
        return handles.A;
    }

    /**
     * Get residence time handles with fluid-specific validation
     * 
     * @return AvgHandle for residence time metrics
     */
    @Override
    public AvgHandle getAvgResidTHandles() {
        SolverAvgHandles handles = getAvgHandles();
        if (handles.W == null || handles.W.isEmpty()) {
            throw new RuntimeException("Residence time handles not available");
        }
        return handles.W;
    }

    /**
     * Get cumulative distribution function for passage time.
     *
     * <p>Delegates verbatim to {@link #getCdfRespT()}, as the reference
     * {@code @SolverFLD/getCdfPassT.m} does: for the fluid solver the passage
     * time and the response time are the same second ODE solve from the
     * steady-state fixed point.</p>
     *
     * @return DistributionResult containing the passage time CDF data
     */
    @Override
    public DistributionResult getCdfPassT() {
        return getCdfRespT();
    }

    /**
     * Get cumulative distribution function for passage time with specific response time handles.
     *
     * <p>The handle only names which metrics the caller wants; the passage-time
     * solve computes every (station, class) anyway, exactly as the reference,
     * so it is accepted for signature compatibility and not read.</p>
     *
     * @param R the response time handles, accepted for signature compatibility
     * @return DistributionResult containing the passage time CDF data
     */
    @Override
    public DistributionResult getCdfPassT(AvgHandle R) {
        return getCdfRespT();
    }

    /**
     * Backward-compatible name for {@link #getCdfPassT()}, to which this
     * delegates -- the {@code @SolverFLD/getCdfPT.m} twin.
     *
     * @return DistributionResult containing the passage time CDF data
     */
    public DistributionResult getCdfPT() {
        return getCdfPassT();
    }

    /** Backward-compatible name for {@link #getCdfPassT(AvgHandle)}. */
    public DistributionResult getCdfPT(AvgHandle R) {
        return getCdfPassT(R);
    }

    /**
     * Run Age of Information analysis if the model is a valid AoI topology.
     * Called during runAnalyzer() when method is "mfq" or after other methods.
     * Uses the MFQ-based solvers (solveBufferless / solveSingleBuffer) to compute
     * matrix exponential representations of AoI and Peak AoI distributions.
     *
     * @param sn Network structure
     */
    private void runAoiAnalysis(NetworkStruct sn) {
        try {
            AoiValidationResult aoiInfo = Aoi_is_aoi.aoi_is_aoi(sn);
            if (!aoiInfo.isAoI()) {
                return;  // Not an AoI topology, skip silently
            }

            // Determine preemption probability from solver options
            double aoiPreemption = Double.NaN;
            if (options.config != null && !Double.isNaN(options.config.aoi_preemption)) {
                aoiPreemption = options.config.aoi_preemption;
            }

            AoiParams params = Aoi_extract_params.aoi_extract_params(sn, aoiInfo, aoiPreemption);

            AoiMfqResult aoiResult;
            if (params.getSystemType().equals("bufferless")) {
                aoiResult = Aoi_solve_bufferless.aoi_solve_bufferless(
                    params.getTau(), params.getT(), params.getSigma(), params.getS(), params.getP());
            } else {
                aoiResult = Aoi_solve_singlebuffer.aoi_solve_singlebuffer(
                    params.getLambda(), params.getSigma(), params.getS(), params.getR());
            }

            ((FluidResult) this.result).aoiResults = aoiResult;

            line_debug(options.verbose,
                String.format("AoI analysis: Mean AoI=%.4f, Mean PAoI=%.4f, system=%s",
                    aoiResult.getAoiMean(), aoiResult.getPaoiMean(), aoiResult.getSystemType()));
        } catch (Exception e) {
            line_warning("SolverFluid",
                "AoI analysis failed: %s. Standard metrics still available.", e.getMessage());
        }
    }

    /**
     * Get average Age of Information metrics.
     *
     * <p>Returns AoI and Peak AoI statistics computed by the Fluid solver.
     * Requires the model to have a valid AoI topology (single open class,
     * Source-Queue-Sink, capacity 1 or 2, single server, FCFS/LCFS/LCFSPR).</p>
     *
     * @return Map with keys "AoI" and "PAoI", each mapping to a Map with
     *         "mean", "var", "std" entries. Also includes "systemType" and
     *         "preemption" at the top level.
     *         Returns null if no AoI results are available.
     */
    public Map<String, Object> getAvgAoI() {
        if (!this.hasAvgResults()) {
            this.getAvg();
        }

        FluidResult fluidResult = (FluidResult) this.result;
        if (fluidResult.aoiResults == null) {
            // Try running AoI analysis now
            if (this.sn == null) {
                this.sn = this.model.getStruct(true);
            }
            runAoiAnalysis(this.sn);
        }

        if (fluidResult.aoiResults == null) {
            line_warning("SolverFluid",
                "No AoI results available. Ensure model has valid AoI topology and use method='mfq'.");
            return null;
        }

        AoiMfqResult aoi = fluidResult.aoiResults;

        Map<String, Object> aoiStats = new LinkedHashMap<String, Object>();
        aoiStats.put("mean", aoi.getAoiMean());
        aoiStats.put("var", aoi.getAoiVar());
        aoiStats.put("std", Math.sqrt(Math.max(0, aoi.getAoiVar())));

        Map<String, Object> paoiStats = new LinkedHashMap<String, Object>();
        paoiStats.put("mean", aoi.getPaoiMean());
        paoiStats.put("var", aoi.getPaoiVar());
        paoiStats.put("std", Math.sqrt(Math.max(0, aoi.getPaoiVar())));

        Map<String, Object> result = new LinkedHashMap<String, Object>();
        result.put("AoI", aoiStats);
        result.put("PAoI", paoiStats);
        result.put("systemType", aoi.getSystemType());
        result.put("preemption", aoi.getPreemption());

        return result;
    }

    /**
     * Get CDF of Age of Information.
     *
     * <p>Computes the cumulative distribution function of AoI and Peak AoI
     * using matrix exponential representations: F(t) = 1 + g * expm(A*t) * inv(A) * h</p>
     *
     * @param tValues Time values at which to evaluate CDF. If null, uses
     *               automatic range based on mean AoI (0 to 5*mean, 200 points).
     * @return Array of two Matrix objects: [AoI_cdf, PAoI_cdf].
     *         Each is an n x 2 matrix with columns [CDF_values, t_values].
     *         Returns null if no AoI results are available.
     */
    public Matrix[] getCdfAoI(Matrix tValues) {
        if (!this.hasAvgResults()) {
            this.getAvg();
        }

        FluidResult fluidResult = (FluidResult) this.result;
        if (fluidResult.aoiResults == null) {
            if (this.sn == null) {
                this.sn = this.model.getStruct(true);
            }
            runAoiAnalysis(this.sn);
        }

        if (fluidResult.aoiResults == null) {
            line_warning("SolverFluid",
                "No AoI results available for CDF computation.");
            return null;
        }

        AoiMfqResult aoi = fluidResult.aoiResults;

        // Generate time values if not provided
        if (tValues == null || tValues.isEmpty()) {
            double meanAoI = aoi.getAoiMean();
            if (Double.isNaN(meanAoI) || meanAoI <= 0) {
                meanAoI = 1.0;
            }
            double tMax = 5.0 * meanAoI;
            int nPoints = 200;
            tValues = new Matrix(nPoints, 1);
            for (int i = 0; i < nPoints; i++) {
                tValues.set(i, 0, tMax * i / (nPoints - 1.0));
            }
        }

        int n = tValues.getNumRows();

        // Compute AoI CDF: F(t) = 1 - S(t), S(t) = -g * expm(A*t) * inv(A) * h.
        // (g,A,h) is a DENSITY triple -- g is normalized by -g*inv(A)*h, so
        // g*expm(A*t)*h is the density and g*inv(A)^2*h the mean -- and the
        // survival function carries the extra inv(A). Subtracting the density
        // gives a curve that falls before it rises, which is not a CDF.
        Matrix aoiCdf = new Matrix(n, 2);
        Matrix aoiG = aoi.getAoiG();
        Matrix aoiA = aoi.getAoiA();
        Matrix aoiH = aoi.getAoiH();
        Matrix aoiAinvH = aoiA.inv().mult(aoiH);

        for (int i = 0; i < n; i++) {
            double t = tValues.get(i, 0);
            aoiCdf.set(i, 1, t);
            if (t <= 0) {
                aoiCdf.set(i, 0, 0.0);
            } else {
                Matrix expAt = aoiA.scale(t).expm();
                double surv = aoiG.mult(expAt).mult(aoiAinvH).get(0, 0);
                aoiCdf.set(i, 0, Math.max(0, Math.min(1, 1.0 + surv)));
            }
        }

        // Compute Peak AoI CDF: the same survival form as above
        Matrix paoiCdf = new Matrix(n, 2);
        Matrix paoiG = aoi.getPaoiG();
        Matrix paoiA = aoi.getPaoiA();
        Matrix paoiH = aoi.getPaoiH();
        Matrix paoiAinvH = paoiA.inv().mult(paoiH);

        for (int i = 0; i < n; i++) {
            double t = tValues.get(i, 0);
            paoiCdf.set(i, 1, t);
            if (t <= 0) {
                paoiCdf.set(i, 0, 0.0);
            } else {
                Matrix expAt = paoiA.scale(t).expm();
                double surv = paoiG.mult(expAt).mult(paoiAinvH).get(0, 0);
                paoiCdf.set(i, 0, Math.max(0, Math.min(1, 1.0 + surv)));
            }
        }

        return new Matrix[]{aoiCdf, paoiCdf};
    }

    /**
     * Get CDF of Age of Information with automatic time range.
     *
     * @return Array of two Matrix objects: [AoI_cdf, PAoI_cdf]
     */
    public Matrix[] getCdfAoI() {
        return getCdfAoI(null);
    }

    /**
     * Get sojourn time CDF. Alias for getCdfRespT().
     *
     * @return DistributionResult containing response time CDFs
     */
    public DistributionResult getSjrnT() {
        return this.getCdfRespT();
    }

    /**
     * Get sojourn time CDF. Lowercase alias for getSjrnT.
     *
     * @return DistributionResult containing response time CDFs
     */
    public DistributionResult sjrnT() {
        return this.getSjrnT();
    }

    /**
     * Bundled third-party libraries used by the fluid solver: rmf_tool backs
     * the refined mean-field methods. Mirrors SolverFLD.getLibrariesUsed in
     * MATLAB.
     */
    @Override
    public java.util.List<String> getLibrariesUsed(NetworkStruct sn, SolverOptions options) {
        java.util.List<String> libs = new java.util.ArrayList<String>();
        String method = (options == null || options.method == null) ? "" : options.method;
        if (method.equals("rmf") || method.equals("fluid.rmf")) {
            libs.add("rmf_tool");
        }
        return libs;
    }

    /**
     * Mark the (station, class) pairs the model actually routes a job into,
     * read off the per-chain visit ratios sn.visits.
     *
     * A fluid result cannot decide that question from the SIZE of QN or TN.
     * Both carry a decaying remnant of the initial state, spread over pairs the
     * class never reaches, and the remnant is whatever the integrator left
     * behind when it stopped: measured at QN = 1.3e-12 and TN = 1.3e-13 on
     * picard05 for test_CQN_Cox_CS_7, i.e. ABOVE GlobalConstants.Zero, so a
     * threshold on them divides one remnant by the other and reports the
     * station's own service time, 10.0000086, as a response time. The visit
     * ratios come from the routing solve instead, where an unrouted pair is
     * zero to the last bits (2.7e-17 on that pair).
     *
     * sn.visits is indexed by STATEFUL node, hence the stationToStateful
     * lookup. A struct carrying no visit information decides nothing and every
     * pair is reported visited. Mirrors fluid_visited_pairs.m.
     */
    private static boolean[][] fluidVisitedPairs(NetworkStruct sn, int M, int K) {
        boolean[][] visited = new boolean[M][K];
        boolean haveVisits = false;
        int maxCols = 0;
        if (sn.visits != null) {
            for (Matrix Vc : sn.visits.values()) {
                if (Vc == null || Vc.isEmpty()) {
                    continue;
                }
                haveVisits = true;
                maxCols = Math.max(maxCols, Vc.getNumCols());
                for (int i = 0; i < M; i++) {
                    // stationToStateful is a row vector here, read linearly as elsewhere
                    boolean haveMap = sn.stationToStateful != null
                            && i < sn.stationToStateful.length();
                    int isf = haveMap ? (int) sn.stationToStateful.get(i) : i;
                    if (isf < 0 || isf >= Vc.getNumRows()) {
                        for (int k = 0; k < K; k++) {
                            visited[i][k] = true;
                        }
                        continue;
                    }
                    for (int k = 0; k < K && k < Vc.getNumCols(); k++) {
                        if (Math.abs(Vc.get(isf, k)) > GlobalConstants.Zero) {
                            visited[i][k] = true;
                        }
                    }
                }
            }
        }
        if (!haveVisits) {
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    visited[i][k] = true;
                }
            }
        } else if (maxCols < K) {
            // A class NO visit matrix reaches is not evidence of a non-visit, only of
            // a struct whose visits were refreshed against fewer classes.
            for (int i = 0; i < M; i++) {
                for (int k = maxCols; k < K; k++) {
                    visited[i][k] = true;
                }
            }
        }
        return visited;
    }
}
