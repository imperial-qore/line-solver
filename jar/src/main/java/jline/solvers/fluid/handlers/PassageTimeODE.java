/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers.fluid.handlers;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.moments.FluidRateFactors;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.Map;
import java.util.Objects;

import static jline.api.mam.Map_pie.map_pie;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;
import static jline.lang.constant.SchedStrategy.DPS;
import static jline.lang.constant.SchedStrategy.GPS;
import static jline.util.Maths.softmin;
import static org.apache.commons.math3.util.FastMath.min;

public class PassageTimeODE implements FirstOrderDifferentialEquations {
    // Softmin method implemented - provides smooth approximation of min function for ODE stability
    private final NetworkStruct sn;
    private final Map<Station, Map<JobClass, Matrix>> mu;
    private final Map<Station, Map<JobClass, Matrix>> phi;
    private final Map<Station, Map<JobClass, MatrixCell>> proc;
    private final Matrix rt;
    private final Matrix nservers;
    private final SolverOptions options;
    private final int numDimensions;

    // Precomputed structures (computed once in constructor, reused in computeDerivatives)
    private final boolean[][] cachedEnabled;
    private final Matrix cachedQIndices;
    private final Matrix cachedKic;
    private final Matrix cachedW;
    // Precomputed for the default (closing) method only
    private final Matrix cachedAllJumps;
    /**
     * Projector onto the coordinates that survive the immediate elimination, null when nothing was
     * eliminated. The caller must apply it to the initial point: mass parked on an eliminated
     * coordinate has no event left to move it.
     */
    private final Matrix immediateAbsorb;
    /**
     * [nEventsReduced x nEventsOriginal] expected firings of each ORIGINAL event per firing of each
     * reduced one, the identity when nothing was eliminated. A caller that classifies events -- which
     * (station,class) each is a completion of -- keeps that classification on the original indexing
     * and maps it here, which is what lets the moment closure read throughputs off a reduced event
     * set. See ImmediateElimination.
     */
    private final Matrix immediateEmap;
    /** Source coordinate of each ORIGINAL event, before any elimination. */
    private final Matrix originalEventIdx;
    /** Number of leading events of the ORIGINAL list that are service completions. */
    private int originalDepartureCount;
    private final Matrix cachedRateBase;
    private final Matrix cachedEventIdx;
    /**
     * Optional time-varying per-event rate multiplier m(t). Null when no
     * time-varying channel is configured, in which case the closure is
     * numerically identical to the legacy autonomous form.
     */
    private final FluidRateMultiplier cachedRateMult;
    /**
     * Per-coordinate service shares of the closing drift, shared with the
     * moment-closure Jacobian so the two can never disagree.
     */
    private final FluidRateFactors cachedFactors;

    public PassageTimeODE(
            NetworkStruct sn,
            Map<Station, Map<JobClass, Matrix>> mu,
            Map<Station, Map<JobClass, Matrix>> phi,
            Map<Station, Map<JobClass, MatrixCell>> proc,
            Matrix rt,
            Matrix S,
            SolverOptions options,
            int numDimensions) {
        this.sn = sn;
        this.mu = mu;
        this.phi = phi;
        this.proc = proc;
        this.rt = rt;
        this.nservers = S;
        this.options = options;
        this.numDimensions = numDimensions;

        // Precompute enabled, qIndices, Kic, w (these don't change between ODE steps)
        int M = S.length();
        int K = mu.get(sn.stations.get(0)).size();

        this.cachedEnabled = new boolean[M][K];
        this.cachedQIndices = new Matrix(M, K);
        this.cachedKic = new Matrix(M, K);
        this.cachedW = new Matrix(M, K);
        int cumSum = 0;

        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int c = 0; c < K; c++) {
                JobClass jobClass = sn.jobclasses.get(c);
                cachedW.set(i, c, 1);
                cachedEnabled[i][c] = false;
                int numPhases = 0;

                Matrix muMatrix = mu.get(station).get(jobClass);
                if (muMatrix != null) {
                    int numNans = 0;
                    for (int row = 0; row < muMatrix.getNumRows(); row++) {
                        for (int col = 0; col < muMatrix.getNumCols(); col++) {
                            if (Double.isNaN(muMatrix.get(row, col))) {
                                numNans++;
                            }
                        }
                    }
                    if ((numNans != muMatrix.getNumElements()) && !muMatrix.isEmpty()) {
                        numPhases = muMatrix.length();
                        cachedEnabled[i][c] = true;
                    }
                }

                cachedQIndices.set(i, c, cumSum);
                cachedKic.set(i, c, numPhases);
                cumSum += numPhases;
            }

            if (sn.sched.get(station) == DPS || sn.sched.get(station) == GPS) {
                for (int k = 0; k < K; k++) {
                    cachedW.set(i, k, sn.schedparam.get(i, k));
                }
                double sumWI = cachedW.sumRows(i);
                if (sumWI > 0) {
                    for (int k = 0; k < K; k++) {
                        cachedW.set(i, k, cachedW.get(i, k) / sumWI);
                    }
                }
            }
        }

        SchedStrategy[] schedArray = new SchedStrategy[M];
        for (int i = 0; i < M; i++) {
            schedArray[i] = sn.sched.get(sn.stations.get(i));
        }
        this.cachedFactors = new FluidRateFactors(M, K, this.cachedEnabled, this.cachedQIndices,
                this.cachedKic, S, this.cachedW, schedArray, sn.lldscaling);

        // Precompute allJumps, rateBase, eventIdx for the default (closing) method
        if (!Objects.equals(options.method, "statedep") && !Objects.equals(options.method, "softmin")) {
            Matrix tmpAllJumps = calculateJumps(cachedEnabled, cachedQIndices, cachedKic);
            Matrix tmpRateBase = new Matrix(tmpAllJumps.getNumCols(), 1);
            Matrix tmpEventIdx = new Matrix(tmpAllJumps.getNumCols(), 1);
            calculateRateBaseAndEventIdxs(cachedEnabled, cachedQIndices, cachedKic, tmpRateBase, tmpEventIdx);

            // Stochastic-complement the instantaneous coordinates out of the event set, so no
            // integrator has to step through an InfRate mode (matching MATLAB solver_fluid_odes.m).
            if (FluidHideImmediate.resolve(sn, options)) {
                ImmediateElimination.EliminationResult result =
                        ImmediateElimination.eliminateImmediate(tmpAllJumps, tmpRateBase, tmpEventIdx, sn, options);
                this.cachedAllJumps = result.allJumpsReduced;
                this.cachedRateBase = result.rateBaseReduced;
                this.cachedEventIdx = result.eventIdxReduced;
                this.immediateAbsorb = result.absorb;
                this.immediateEmap = result.Emap;
                this.originalEventIdx = tmpEventIdx;
            } else {
                this.cachedAllJumps = tmpAllJumps;
                this.cachedRateBase = tmpRateBase;
                this.cachedEventIdx = tmpEventIdx;
                this.immediateAbsorb = null;
                this.immediateEmap = null;
                this.originalEventIdx = tmpEventIdx;
            }

            // Optional time-varying event-rate multiplier m(t), making the
            // otherwise autonomous closing ODE non-autonomous. Built from the
            // post-elimination event mapping so the two agree event by event;
            // null when no channel is configured. Mirrors the rate multiplier
            // composition of MATLAB solver_fluid_odes.m.
            this.cachedRateMult = FluidRateMultiplier.build(
                    this.cachedEventIdx.getNumRows(), this.cachedEnabled, this.cachedQIndices,
                    this.cachedKic, mu, sn.stations, sn.jobclasses, this.cachedEventIdx, options);
        } else {
            this.cachedAllJumps = null;
            this.cachedRateBase = null;
            this.cachedEventIdx = null;
            this.cachedRateMult = null;
            this.immediateAbsorb = null;
            this.immediateEmap = null;
            this.originalEventIdx = null;
        }
    }

    public PassageTimeODE(
            NetworkStruct sn,
            Map<Station, Map<JobClass, Matrix>> mu,
            Map<Station, Map<JobClass, Matrix>> phi,
            Map<Station, Map<JobClass, MatrixCell>> proc,
            Matrix rt,
            Matrix S,
            SolverOptions options) {
        this(sn, mu, phi, proc, rt, S, options, options.init_sol.length());
    }

    private Matrix calculateJumps(boolean[][] enabled, Matrix qIndices, Matrix Kic) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        int jumpsRows = (int) Kic.elementSum();
        Matrix jumps =
                new Matrix(jumpsRows, 0); // Returns state changes triggered by all the events

        for (int i = 0; i < M; i++) { // state changes from departures in service phases 2
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    int xic = (int) qIndices.get(i, c); //  index of x_ic
                    for (int j = 0; j < M; j++) {
                        for (int l = 0; l < K; l++) {
                            if (rt.get(i * K + c, j * K + l) > 0) {
                                int xjl = (int) qIndices.get(j, l); // index of x_jl
                                for (int ki = 0; ki < Kic.get(i, c); ki++) { // job can leave from any phase in i
                                    for (int kj = 0; kj < Kic.get(j, l); kj++) { // job can start from any phase in j
                                        setNextJump(jumps, xic + ki, xjl + kj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) { // state changes: "next service phase" transition
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    int xic = (int) qIndices.get(i, c);
                    // every source phase, the last included: bounding ki at Kic-1 is valid
                    // only for an acyclic PH and drops the last row of D0 for a general MAP
                    // or MMPP2, whose D0 is cyclic
                    for (int ki = 0; ki < Kic.get(i, c); ki++) {
                        for (int kip = 0; kip < Kic.get(i, c); kip++) {
                            if (ki != kip) {
                                setNextJump(jumps, xic + ki, xic + kip);
                            }
                        }
                    }
                }
            }
        }

        return jumps;
    }

    private void calculateRateBaseAndEventIdxs(
            boolean[][] enabled,
            Matrix qIndices,
            Matrix Kic,
            Matrix rateBase,
            Matrix eventIdx) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        int rateIdx = 0;

        // State changes from departures in service phases 2...
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    for (int j = 0; j < M; j++) {
                        for (int l = 0; l < K; l++) {
                            Matrix pie;
                            if (proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).isEmpty()) {
                                pie = new Matrix(1, 1, 1);
                                pie.set(0, 0, 1);
                            } else {
                                Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);

                                pie = map_pie(D0, D1);
                            }
                            if (rt.get(i * K + c, j * K + l) > 0) {
                                for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                    for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                        rateBase.set(
                                                rateIdx,
                                                0,
                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                        * mu.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                        * rt.get(i * K + c, j * K + l)
                                                        * pie.get(0, kjl));
                                        eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                                        rateIdx++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Everything emitted so far is a service COMPLETION; what follows is an
        // intra-PH phase change. The boundary is what lets a caller say which
        // (station,class) an event is a completion of.
        this.originalDepartureCount = rateIdx;

        // State changes from "next service phase" transition in phases 2...
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    // must match the ki range of calculateJumps event for event
                    for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                        for (int kicp = 0; kicp < Kic.get(i, c); kicp++) {
                            if (kicp != kicIdx) {
                                rateBase.set(
                                        rateIdx,
                                        0,
                                        proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicp));
                                eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                                rateIdx++;
                            }
                        }
                    }
                }
            }
        }
    }

    private Matrix calculatedxdtClosingMethod(
            double t,
            double[] x,
            Matrix w,
            boolean[][] enabled,
            Matrix qIndices,
            Matrix Kic,
            Matrix allJumps,
            Matrix rateBase,
            Matrix eventIdx) {

        Matrix rates = computeClosingRatesVector(x, w, enabled, qIndices, Kic);

        int numEventIndices = eventIdx.getNumRows(); // Use getNumRows(), not length() which returns max(rows,cols)
        Matrix newRates = new Matrix(numEventIndices, 1);
        for (int i = 0; i < numEventIndices; i++) {
            newRates.set(i, 0, rates.get((int) eventIdx.get(i, 0), 0));
        }
        newRates.elementMult(rateBase, newRates);
        applyRateMultiplier(t, newRates);

        return allJumps.mult(newRates, null);
    }

    /**
     * Scales the closing event rates in place by the time-varying multiplier
     * m(t), if one is configured. A no-op otherwise, which keeps the legacy
     * autonomous closure numerically unchanged.
     */
    /**
     * Whether the drift is AUTONOMOUS, i.e. depends on the state alone.
     *
     * <p>A rate schedule (NHPP, MAPt, PHt) makes the right-hand side a function of
     * t as well, and then a zero residual at one instant says nothing about the
     * next segment. The fixed-point short circuit in
     * ClosingAndStateDepMethodsAnalyzer is armed only when this is true.</p>
     */
    public boolean isAutonomous() {
        return cachedRateMult == null;
    }

    private void applyRateMultiplier(double t, Matrix newRates) {
        if (cachedRateMult == null) {
            return;
        }
        double[] mult = cachedRateMult.evalAt(t);
        for (int i = 0; i < mult.length; i++) {
            newRates.set(i, 0, newRates.get(i, 0) * mult[i]);
        }
    }

    /**
     * Raw closing-method event-rate vector, mirroring MATLAB
     * ode_rates_closing: the per-index scheduling correction gathered over the
     * event index set and scaled by the rate base. The returned vector has one
     * entry per event of the (possibly immediate-eliminated) closing machinery
     * and is left-multiplied by the jump matrix to obtain the derivative. It is
     * exposed for the trajectory-based iteration (TBI) analyzer, which applies
     * the row-restricted jump matrix of each cell to this full rate vector.
     *
     * @param x full fluid state (length equal to getDimension())
     * @return column vector of event rates (numEvents x 1)
     */
    public Matrix calculateRatesClosing(double[] x) {
        Matrix w = cachedW.copy();
        Matrix rates = computeClosingRatesVector(x, w, cachedEnabled, cachedQIndices, cachedKic);

        int numEventIndices = cachedEventIdx.getNumRows();
        Matrix newRates = new Matrix(numEventIndices, 1);
        for (int i = 0; i < numEventIndices; i++) {
            newRates.set(i, 0, rates.get((int) cachedEventIdx.get(i, 0), 0));
        }
        newRates.elementMult(cachedRateBase, newRates);
        return newRates;
    }

    /**
     * Closing-method event-rate vector at time {@code t}, that is
     * {@link #calculateRatesClosing(double[])} scaled by the time-varying rate
     * multiplier m(t) when one is configured. Identical to the autonomous form
     * otherwise.
     *
     * @param t current time
     * @param x full fluid state (length equal to getDimension())
     * @return column vector of event rates (numEvents x 1)
     */
    /**
     * Projector onto the coordinates that survive the immediate elimination, null when nothing was
     * eliminated. Applied to the initial point by the analyzer: once the instantaneous coordinates
     * are complemented away no event moves them any more, so whatever mass the initial condition
     * parked there -- a cold start puts everything in phase 1, but a warm start from an earlier LN
     * iterate does not -- would be frozen for the whole integration and lost from its chain.
     *
     * @return the [nStates x nStates] projector, or null
     */
    public Matrix getImmediateAbsorb() {
        return immediateAbsorb;
    }

    /**
     * Expected firings of each ORIGINAL event per firing of each reduced one, null when nothing was
     * eliminated (which the caller reads as the identity).
     *
     * @return the [nEventsReduced x nEventsOriginal] map, or null
     */
    public Matrix getImmediateEmap() {
        return immediateEmap;
    }

    /**
     * Source coordinate of each event BEFORE the immediate elimination. Event attributes are
     * classified on this indexing and mapped onto the reduced events through getImmediateEmap.
     *
     * @return the original [nEventsOriginal x 1] event index vector
     */
    public Matrix getOriginalEventIdx() {
        return originalEventIdx;
    }

    /**
     * Number of leading events of the ORIGINAL event list that are service completions, the rest
     * being intra-PH phase changes.
     *
     * @return the departure-event count
     */
    public int getOriginalDepartureCount() {
        return originalDepartureCount;
    }

    public Matrix calculateRatesClosing(double t, double[] x) {
        Matrix newRates = calculateRatesClosing(x);
        applyRateMultiplier(t, newRates);
        return newRates;
    }

    /** Precomputed jump matrix of the closing method (dimension x numEvents). */
    public Matrix getAllJumps() {
        return cachedAllJumps;
    }

    /** Precomputed per-(station,class) starting state index of the closing method. */
    public Matrix getQIndices() {
        return cachedQIndices;
    }

    /** Precomputed per-(station,class) phase count of the closing method. */
    public Matrix getKic() {
        return cachedKic;
    }

    /** Precomputed per-(station,class) enabled flags of the closing method. */
    public boolean[][] getEnabled() {
        return cachedEnabled;
    }

    /**
     * Per-coordinate service shares of the closing drift, delegated to
     * {@link FluidRateFactors} so that the ODE the solver integrates and the
     * Jacobian the moment-closure covariance equation reads are literally the
     * same expression. The moment closure enters through
     * {@code options.config.moment_sigma2} and {@code options.config.moment_cov};
     * with both absent the shares are the first-order (plug-in) ones and the
     * legacy code path is bit-identical.
     */
    private Matrix computeClosingRatesVector(
            double[] x,
            Matrix w,
            boolean[][] enabled,
            Matrix qIndices,
            Matrix Kic) {

        double[] sigma2 = (options.config == null) ? null : options.config.moment_sigma2;
        Matrix[] covblk = (options.config == null) ? null : options.config.moment_cov;
        return cachedFactors.factors(x, sigma2, covblk);
    }

    /** Per-coordinate service shares evaluated at an explicit closure variance. */
    public Matrix calculateFactors(double[] x, double[] sigma2, Matrix[] covblk) {
        return cachedFactors.factors(x, sigma2, covblk);
    }

    /** Analytic Jacobian of the per-coordinate service shares. */
    public Matrix calculateFactorsJacobian(double[] x, double[] sigma2, Matrix[] covblk) {
        return cachedFactors.jacobian(x, sigma2, covblk);
    }

    /** Precomputed per-event constant rate factor of the closing method. */
    public Matrix getRateBase() {
        return cachedRateBase;
    }

    /** Precomputed per-event source state coordinate of the closing method. */
    public Matrix getEventIdx() {
        return cachedEventIdx;
    }

    /** Shared per-coordinate service share evaluator of the closing method. */
    public FluidRateFactors getRateFactors() {
        return cachedFactors;
    }


    private Matrix calculatedxdtStateDepMethod(
            double[] x, boolean[][] enabled, Matrix qIndices, Matrix Kic, Matrix w) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        Matrix dxdt = new Matrix(x.length, 1);

        // Declare variables outside switch to avoid scope issues
        int idxIni, idxEnd;
        double ni;
        
        for (int i = 0; i < M; i++) {
            switch (sn.sched.get(sn.stations.get(i))) {
                case INF:
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie =
                                                map_pie(D0, D1);
                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    if (j != i) {
                                                        double rate =
                                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                        * mu.get(sn.stations.get(i))
                                                                        .get(sn.jobclasses.get(c))
                                                                        .get(kicIdx, 0)
                                                                        * rt.get(i * K + c, j * K + l)
                                                                        * pie.get(0, kjl);
                                                        dxdt.set(
                                                                xic + kicIdx,
                                                                0,
                                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                        dxdt.set(
                                                                xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case EXT:
                    // MATLAB ode_statedep.m: state-dep/softmin ODE does not support open models
                    line_error(mfilename(new Object(){}),
                            "State dependent ODE method does not support open models. Try with default method.");
                    break;

                case PS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        if (ni > sn.nservers.get(i, 0)) {
                                            dxdt.set(
                                                    xic + kicIdx,
                                                    0,
                                                    dxdt.get(xic + kicIdx, 0)
                                                            - (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                                            dxdt.set(
                                                    xic + kic_p,
                                                    0,
                                                    dxdt.get(xic + kic_p, 0)
                                                            + (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                                        } else {
                                            dxdt.set(
                                                    xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                            dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    if (ni > sn.nservers.get(i, 0)) {
                                                        rate = 1/ni * sn.nservers.get(i, 0) * rate;
                                                    }
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case FCFS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    double wni = GlobalConstants.FineTol;
                    for (int c = 0; c < K; c++) {
                        for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                            if (enabled[i][c]) {
                                int xic = (int) qIndices.get(i, c);
                                w.set(
                                        c,
                                        kicIdx,
                                        -1
                                                / proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicIdx));
                                wni += w.get(c, kicIdx) * x[xic + kicIdx];
                            }
                        }
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p)
                                                        * min(ni, sn.nservers.get(i, 0))
                                                        * w.get(c, kicIdx)
                                                        / wni;
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl)
                                                                    * min(ni, sn.nservers.get(i, 0))
                                                                    * w.get(c, kicIdx)
                                                                    / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case DPS:
                    double sumWI = w.sumRows(i);
                    for (int col = 0; col < w.getNumCols(); col++) {
                        w.set(i, col, w.get(i, col) / sumWI);
                    }
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    wni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        wni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        if (wni > sn.nservers.get(i, 0)) {
                                            dxdt.set(
                                                    xic + kicIdx,
                                                    0,
                                                    dxdt.get(xic + kicIdx, 0)
                                                            - (x[xic + kicIdx]
                                                            * rate
                                                            * sn.nservers.get(i, 0)
                                                            * w.get(c, kicIdx)
                                                            / wni));
                                            dxdt.set(
                                                    xic + kic_p,
                                                    0,
                                                    dxdt.get(xic + kic_p, 0)
                                                            + (x[xic + kicIdx]
                                                            * rate
                                                            * sn.nservers.get(i, 0)
                                                            * w.get(c, kicIdx)
                                                            / wni));
                                        } else {
                                            dxdt.set(
                                                    xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                            dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    if (wni > sn.nservers.get(i, 0)) {
                                                        rate *= sn.nservers.get(i, 0) * w.get(c, kicIdx) / wni;
                                                    }
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
            }
        }

        return dxdt;
    }

    private Matrix calculatedxdtSoftminMethod(
            double[] x, boolean[][] enabled, Matrix qIndices, Matrix Kic, Matrix w) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        Matrix dxdt = new Matrix(x.length, 1);
        double alpha = 20.0; // Softmin smoothing parameter (as used in MATLAB implementation)
        
        // Declare variables outside switch to avoid scope issues
        int idxIni, idxEnd;
        double ni;

        for (int i = 0; i < M; i++) {
            switch (sn.sched.get(sn.stations.get(i))) {
                case INF:
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);
                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    if (j != i) {
                                                        double rate =
                                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                        * mu.get(sn.stations.get(i))
                                                                        .get(sn.jobclasses.get(c))
                                                                        .get(kicIdx, 0)
                                                                        * rt.get(i * K + c, j * K + l)
                                                                        * pie.get(0, kjl);
                                                        dxdt.set(
                                                                xic + kicIdx,
                                                                0,
                                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                        dxdt.set(
                                                                xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case EXT:
                    // MATLAB ode_statedep.m: state-dep/softmin ODE does not support open models
                    line_error(mfilename(new Object(){}),
                            "Softmin ODE method does not support open models. Try with default method.");
                    break;

                case PS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        // THE PS BRANCH TAKES THE HARD SHARE, NOT THE
                                        // SMOOTH ONE. ode_softmin.m:71-77 keeps the same
                                        // `if ni > nservers` scaling ode_statedep uses here
                                        // and applies softmin only in the FCFS and DPS
                                        // branches below, where the division is by wni,
                                        // which is seeded with FineTol and can never be
                                        // zero. Dividing by ni here instead produced NaN at
                                        // t = 0 on any model whose queue starts empty, and
                                        // the integrator failed on the first step.
                                        double psFactor = (ni > sn.nservers.get(i, 0))
                                                ? sn.nservers.get(i, 0) / ni : 1.0;
                                        dxdt.set(
                                                xic + kicIdx,
                                                0,
                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate * psFactor));
                                        dxdt.set(
                                                xic + kic_p,
                                                0,
                                                dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate * psFactor));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    // The hard PS share, as ode_softmin.m:95-97
                                                    // takes it; see the phase-change site above.
                                                    if (ni > sn.nservers.get(i, 0)) {
                                                        rate = rate * sn.nservers.get(i, 0) / ni;
                                                    }
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case FCFS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    double wni = GlobalConstants.FineTol;
                    for (int c = 0; c < K; c++) {
                        for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                            if (enabled[i][c]) {
                                int xic = (int) qIndices.get(i, c);
                                w.set(
                                        c,
                                        kicIdx,
                                        -1
                                                / proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicIdx));
                                wni += w.get(c, kicIdx) * x[xic + kicIdx];
                            }
                        }
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p)
                                                        * softmin(ni, sn.nservers.get(i, 0), alpha) // Use softmin here
                                                        * w.get(c, kicIdx)
                                                        / wni;
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl)
                                                                    * softmin(ni, sn.nservers.get(i, 0), alpha) // Use softmin here
                                                                    * w.get(c, kicIdx)
                                                                    / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case DPS:
                    double sumWI = w.sumRows(i);
                    for (int col = 0; col < w.getNumCols(); col++) {
                        w.set(i, col, w.get(i, col) / sumWI);
                    }
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    wni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        wni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        // Use softmin for DPS scheduling
                                        double softminFactor = softmin(wni, sn.nservers.get(i, 0), alpha) * w.get(c, kicIdx) / wni;
                                        dxdt.set(
                                                xic + kicIdx,
                                                0,
                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate * softminFactor));
                                        dxdt.set(
                                                xic + kic_p,
                                                0,
                                                dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate * softminFactor));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    // Use softmin for capacity constraint
                                                    rate *= softmin(wni, sn.nservers.get(i, 0), alpha) * w.get(c, kicIdx) / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;
                    
                default:
                    // For other scheduling strategies, fall back to the state-dependent method
                    line_error(mfilename(new Object(){}),
                            "Softmin ODE method does not support scheduling strategy: " + sn.sched.get(sn.stations.get(i)) + 
                            ". Try with statedep method.");
                    break;
            }
        }

        return dxdt;
    }

    @Override
    public void computeDerivatives(double t, double[] x, double[] dxdt)
            throws MaxCountExceededException, DimensionMismatchException {

        // The drift is evaluated at the integrator's own state, as MATLAB
        // solver_fluid_iteration does: odeset('NonNegative') projects the ACCEPTED
        // step, it does not project the argument of the right-hand side. Clamping
        // here instead makes the drift discontinuous at x=0, which collapses the
        // step size (see _kb/06-solver-catalog.md, PassageTimeODE nonnegativity).

        // Use precomputed structures (w is copied since some methods modify it in-place)
        Matrix w = cachedW.copy();

        Matrix dxdtTmp;
        if (Objects.equals(options.method, "statedep")) {
            dxdtTmp = calculatedxdtStateDepMethod(x, cachedEnabled, cachedQIndices, cachedKic, w);
        } else if (Objects.equals(options.method, "softmin")) {
            dxdtTmp = calculatedxdtSoftminMethod(x, cachedEnabled, cachedQIndices, cachedKic, w);
        } else {
            // Use precomputed (and potentially immediate-eliminated) allJumps/rateBase/eventIdx
            dxdtTmp =
                    calculatedxdtClosingMethod(t, x, w, cachedEnabled, cachedQIndices, cachedKic,
                            cachedAllJumps, cachedRateBase, cachedEventIdx);
        }

        for (int i = 0; i < dxdt.length; i++) {
            dxdt[i] = dxdtTmp.get(i);
        }
    }

    @Override
    public int getDimension() {
        return numDimensions;
    }

    private void setNextJump(Matrix jumps, int completionIdx, int startIdx) {

        int jumpsCols = jumps.getNumCols();
        jumps.expandMatrix(jumps.getNumRows(), jumpsCols + 1, jumps.getNumElements() + 2);
        // ACCUMULATED, not assigned: a self-transition has both ends on one
        // coordinate and must net to zero, so a station with a routing self-loop
        // would otherwise create population out of nothing. Mirrors MATLAB
        // ode_jumps_new.
        jumps.set(completionIdx, jumpsCols,
                jumps.get(completionIdx, jumpsCols) - 1); // type c in stat i completes service
        jumps.set(startIdx, jumpsCols,
                jumps.get(startIdx, jumpsCols) + 1); // type c job starts in stat j
    }
}
