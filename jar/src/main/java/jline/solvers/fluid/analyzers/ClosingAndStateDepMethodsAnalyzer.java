/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.Map_mean;
import jline.api.sn.SnRtStations;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.FluidNhpp;
import jline.solvers.fluid.handlers.PassageTimeODE;
import jline.solvers.fluid.moments.FluidConservationGuard;
import jline.solvers.fluid.handlers.TransientDataHandler;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import odesolver.LSODA;
import odesolver.StepStopCondition;
import org.apache.commons.math3.ode.events.EventHandler;

public class ClosingAndStateDepMethodsAnalyzer implements FluidAnalyzer {

    public Matrix xvec_t;
    public Matrix xvec_it;
    /** The ODE built by the last solver_fluid call, carrying the immediate-elimination maps. */
    protected PassageTimeODE lastPassageOde;

    protected void solver_fluid_iteration(NetworkStruct sn,
                                        Map<Station, Map<JobClass, Matrix>> mu,
                                        Map<Station, Map<JobClass, Matrix>> phi,
                                        Matrix S,
                                        double[] yDefault,
                                        Matrix slowrate,
                                        SolverOptions options,
                                        SolverResult result) {
        int iter = 0;
        int totalSteps = 0;

        int slowrateRows = slowrate.getNumRows();
        int slowrateCols = slowrate.getNumCols();
        double minNonZeroRate = GlobalConstants.Inf;

        double maxNonZeroRate = 0.0;
        for (int i = 0; i < slowrateRows; i++) {
            for (int j = 0; j < slowrateCols; j++) {
                if ((slowrate.get(i, j) > GlobalConstants.CoarseTol) && Double.isFinite(slowrate.get(i, j))) {
                    if (slowrate.get(i, j) < minNonZeroRate) {
                        minNonZeroRate = slowrate.get(i, j);
                    }
                    if (slowrate.get(i, j) > maxNonZeroRate) {
                        maxNonZeroRate = slowrate.get(i, j);
                    }
                }
            }
        }

        // NHPP sources store proc as the {breakpoints, rates, cyclic} schedule
        // rather than a {D0, D1} MAP; substitute the equivalent one-phase
        // exponential at the nominal (time-average) rate so the closing-rate
        // builder sees a Markovian representation. The time-varying intensity
        // is applied separately by the per-event rate multiplier.
        Map<Station, Map<JobClass, MatrixCell>> proc = FluidNhpp.substituteNhppProc(sn, mu, phi);
        // sn.rt is indexed by STATEFUL NODE and the drift indexes it by STATION; the
        // two coincide only when no stateful node is a non-station (Cache, Router,
        // Logger). See SnRtStations.
        Matrix rtStations = SnRtStations.snRtStations(sn).getLeft();
        PassageTimeODE passageOde = new PassageTimeODE(sn, mu, phi, proc, rtStations, S, options);
        // Kept for ANALYZE, which reads throughputs off the state and needs the ODE's Emap to add
        // back the completions the eliminated coordinates make.
        this.lastPassageOde = passageOde;
        FirstOrderDifferentialEquations ode = passageOde;

        // THE INITIAL POINT HAS TO BE PROJECTED TOO. Once the instantaneous coordinates are
        // complemented away no event moves them any more, so whatever mass the initial condition
        // parked there would be frozen for the whole integration and lost from its chain. The
        // projector sends it where the eliminated coordinate would have sent it instantaneously.
        Matrix immediateAbsorb = passageOde.getImmediateAbsorb();
        if (immediateAbsorb != null && this.xvec_it != null
                && this.xvec_it.getNumCols() == immediateAbsorb.getNumRows()) {
            Matrix projected = this.xvec_it.mult(immediateAbsorb);
            this.xvec_it = projected;
        }

        double T0 = options.timespan[0];
        double T = 0.0;

        List<Matrix> tIterations = new LinkedList<Matrix>();
        List<Matrix> xVecIterations = new LinkedList<Matrix>();

        boolean goon = true;

        // Only a finite timespan resolves the transient between the endpoints; on an
        // infinite one the per-step trajectory costs O(iter_max * steps) memory to
        // describe a fixed point, so keep just its endpoints (see below)
        boolean keepTrajectory = Double.isFinite(options.timespan[1]);
        Matrix xvecInit = this.xvec_it.copy();

        // Wall-clock budget (options.timeout, seconds; Inf = none), as MATLAB solver_fluid.m
        double maxTime = (Double.isFinite(options.timeout) && options.timeout > 0) ? options.timeout
                : GlobalConstants.Inf;
        long startNanos = System.nanoTime();

        // ARM THE CONSERVATION GUARD ONLY WHERE THE MOMENT CLOSURE IS ACTIVE.
        // MinNormalAnalyzer runs its first pass at sigma2 = 0 -- that pass IS the
        // first-order solve -- and only the LATER passes integrate a drift that
        // can leave the simplex. Testing sigma2 keeps the guard off every
        // first-order call, so 'closing' and 'matrix' (which are also the
        // ladder's own fallback) keep identical behaviour and cannot be sent
        // down a fallback by their own watchdog. See FluidConservationGuard.
        FluidConservationGuard guard = null;
        if (FluidConservationGuard.closureActive(options)) {
            Matrix guardPhases = new Matrix(sn.nstations, sn.nclasses);
            for (int gi = 0; gi < sn.nstations; gi++) {
                Station gStation = sn.stations.get(gi);
                for (int gk = 0; gk < sn.nclasses; gk++) {
                    JobClass gClass = sn.jobclasses.get(gk);
                    guardPhases.set(gi, gk, mu.get(gStation).get(gClass).length());
                }
            }
            // ONLY WHEN THE STATE IS THE ONE THE PHASES DESCRIBE.
            // config.hide_immediate folds immediate coordinates out, which
            // SHRINKS the vector, and the block offsets read off phases would
            // then address the wrong coordinates and trip on a sum that was
            // never that class's population. The lengths agreeing is the exact
            // test for that, so the guard stays off rather than guessing.
            if ((int) guardPhases.elementSum() == yDefault.length) {
                guard = new FluidConservationGuard(sn, guardPhases);
            }
        }

        // Early stop on the GEOMETRIC TAIL of the window iteration; see below.
        boolean earlystop = options.config.fluid_earlystop == null || options.config.fluid_earlystop;
        // The residual cannot go below the error the integrator itself carries.
        double driftTol = FastMath.max(options.iter_tol, options.tol);
        double driftSafety = 0.01;   // headroom, since rho is estimated, not known
        double minHorizon = 10.0 / minNonZeroRate;
        double movedPrev = Double.POSITIVE_INFINITY;
        double[] rhoHist = new double[] { Double.NaN, Double.NaN, Double.NaN };
        int driftBelow = 0;

        jline.io.LineConsole.loop("integrating the fluid ODEs over successive time windows");
        while ((Double.isFinite(options.timespan[1]) && T < options.timespan[1]) || (goon && iter < options.iter_max)) {
            iter++;

            if ((System.nanoTime() - startNanos) / 1e9 > maxTime) {
                goon = false;
                break;
            }

            double[] initialState = new double[yDefault.length];
            double[] nextState = new double[yDefault.length];
            for (int i = 0; i < yDefault.length; i++) {
                initialState[i] = xvec_it.get(0, i);
                nextState[i] = 0.0;
            }
            if (iter == 1) {
                T = FastMath.min(options.timespan[1], FastMath.abs(10.0 / minNonZeroRate));
            } else {
                T = FastMath.min(options.timespan[1], FastMath.abs(10.0 * iter / minNonZeroRate));
            }
            double[] tRange = new double[] { T0, T };

            int Tmax = 0;
            // A FIXED POINT ENDS THE WINDOW IN CLOSED FORM, and this is what keeps a
            // window that has already converged from becoming a window that never
            // returns. Armed only for an AUTONOMOUS drift, so F(y*) = 0 means
            // y(t) = y* for every later t and the rest of the span is known exactly
            // rather than integrated.
            //
            // THE THRESHOLD IS ROUND-OFF, NOT THE SOLVER TOLERANCE. driftTol (1e-4 by
            // default) says "converged to what the caller asked for", and a state that
            // merely satisfies THAT is still moving -- cutting the window there was
            // measured to shift results by 1.8e-5 in the Python twin. A normalized
            // residual below GlobalConstants.Zero is the stronger claim that the drift
            // is zero to double precision, and that is what makes this exact.
            //
            // WHY THE WINDOW DOES NOT END ON ITS OWN. A stiff step controller handed a
            // state it is already at cannot pick a step: on the LN layer of
            // test_LQN_13 the Python twin advanced t by 0.011 in 20000 steps from a
            // state with |F| = 1.5e-16, and covered the whole 1000-unit span in 60
            // steps once that state was nudged 1e-6 off the equilibrium. That layer
            // carries an Immediate() coordinate -- an eigenvalue of exactly
            // -GlobalConstants.Immediate = -1e8 that the immediate elimination did not
            // fold out -- so the controller is pinned near 1/1e8 while the window runs
            // to 10*iter/minNonZeroRate. 272 windows took 3.0 s between them and the
            // 273rd had not returned after 143 s.
            // A FINITE timespan is a transient request, which must reach its end time
            // rather than stop at the fixed point -- the same gate the tail test carries.
            final boolean fpArmed = earlystop && !Double.isFinite(options.timespan[1])
                    && passageOde.isAutonomous();
            boolean settled = false;
            if (fpArmed) {
                settled = fluidSettled(ode, tRange[0], initialState, initialState.length,
                        0, minNonZeroRate);
            }
            // THE ENTRY TEST ABOVE IS TAKEN ONCE, at the window's first instant,
            // and a window that reaches the fixed point AFTER its first step is
            // left grinding out the rest of its span -- the same stall, entered
            // one step later. LSODA therefore carries the SAME test as a
            // per-accepted-step stop condition, which is where MATLAB's OutputFcn
            // and the native-Python step loop take it, so the four codebases stop
            // on the same condition at the same place. Not addStepHandler: that
            // method is an empty stub on LSODA and discards what it is handed.
            final int fpDim = initialState.length;
            final double fpRate = minNonZeroRate;
            final FirstOrderDifferentialEquations fpOde = ode;
            StepStopCondition fpStop = null;
            if (fpArmed) {
                fpStop = new StepStopCondition() {
                    @Override
                    public boolean stop(double t, double[] y) {
                        // y is 1-indexed here: LSODA keeps y[0] unused.
                        return fluidSettled(fpOde, t, y, fpDim, 1, fpRate);
                    }
                };
            }
            if (settled) {
                int dim = initialState.length;
                this.xvec_it = new Matrix(1, dim);
                for (int j = 0; j < dim; j++) {
                    this.xvec_it.set(0, j, FastMath.max(0.0, initialState[j]));
                }
                Tmax = 2;
                if (keepTrajectory) {
                    DMatrixRMaj denseT = new DMatrixRMaj(2, 1);
                    DMatrixRMaj denseX = new DMatrixRMaj(2, dim);
                    denseT.set(0, 0, tRange[0]);
                    denseT.set(1, 0, tRange[1]);
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < dim; j++) {
                            denseX.set(i, j, FastMath.max(0.0, initialState[j]));
                        }
                    }
                    tIterations.add(new Matrix(denseT));
                    xVecIterations.add(new Matrix(denseX));
                }
            } else if (options.stiff) {
                LSODA odeSolver = options.odesolvers.stiffIntegratorFor(tRange[0], tRange[1], options.tol,
                        options.tol > GlobalConstants.CoarseTol);
                // The caller's requested output instants that fall in THIS window;
                // the integrator answers for them from its own interpolant. Set
                // per window (and to null when there are none) because
                // stiffIntegratorFor may hand back a shared instance.
                odeSolver.setOutputTimes(windowOutputTimes(options.tranpoints, tRange[0], tRange[1]));
                // Set per window for the same reason, and to null when unarmed:
                // a condition left over from a previous window would apply here.
                odeSolver.setStopCondition(fpStop);

                try {
                    odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                } catch (RuntimeException firstFailure) {
                    // Do NOT name a cause yet. LSODA raises the same RuntimeException for a
                    // bad initial point and for repeated corrector failure, and its message
                    // ("corrector convergence failed repeatedly or with abs(h)=hmin") names
                    // both alternatives without saying which fired. Re-integrating from
                    // yDefault is the PROBE that decides: only if it succeeds was the
                    // initial point implicated.
                    try {
                        odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState);
                    } catch (RuntimeException retryFailure) {
                        // Both starting points failed, so the initial point was never the
                        // cause. Report that and let it propagate: returning here would
                        // hand back an un-integrated state as if it were a fluid solution.
                        retryFailure.addSuppressed(firstFailure);
                        throw new RuntimeException("Fluid integration failed from BOTH the supplied"
                                + " initial point and the default initialization over t in ["
                                + tRange[0] + ", " + tRange[1] + "]; the initial point is NOT"
                                + " implicated. Underlying integrator error: "
                                + retryFailure.getMessage(), retryFailure);
                    }
                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("The initial point was invalid, Fluid solver switched to default initialization.");
                    }
                }
                if (odeSolver.getStepsTaken() > 0) {
                    java.util.ArrayList<Double> tHistory = odeSolver.getTvec();
                    java.util.ArrayList<Double[]> yHistory = odeSolver.getYvec();
                    int dim = ode.getDimension();
                    // The RECORDED row count, not getStepsTaken()+1: with output
                    // instants requested the history holds the integrator's steps
                    // AND the interpolated rows, so counting steps would truncate
                    // the window and carry an early state into the next one.
                    // Without a request the two are the same number.
                    Tmax = Math.min(tHistory.size(), yHistory.size());
                    int lastIdx = Tmax - 1;
                    this.xvec_it = new Matrix(1, dim);
                    for (int j = 0; j < dim; j++) this.xvec_it.set(0, j, Math.max(0.0, yHistory.get(lastIdx)[j]));
                    if (keepTrajectory) {
                        DMatrixRMaj denseT = new DMatrixRMaj(Tmax, 1);
                        DMatrixRMaj denseX = new DMatrixRMaj(Tmax, dim);
                        for (int i = 0; i < Tmax; i++) {
                            denseT.set(i, 0, tHistory.get(i));
                            for (int j = 0; j < dim; j++) denseX.set(i, j, Math.max(0.0, yHistory.get(i)[j]));
                        }
                        tIterations.add(new Matrix(denseT));
                        xVecIterations.add(new Matrix(denseX));
                    }
                }
            } else {
                FirstOrderIntegrator odeSolver = options.odesolvers.integratorFor(tRange[0], tRange[1],
                        options.tol, options.tol > GlobalConstants.CoarseTol);
                odeSolver.clearStepHandlers();
                TransientDataHandler stepHandler = keepTrajectory
                        ? new TransientDataHandler(initialState.length) : null;
                if (stepHandler != null) {
                    odeSolver.addStepHandler(stepHandler);
                }
                // The same fixed-point stop on the non-stiff arm. Here it is an
                // EventHandler rather than the LSODA hook, because a Commons
                // Math integrator DOES honour events and Action.STOP ends the
                // window without an exception -- which matters, since the
                // catches below turn any RuntimeException into the stiff
                // fallback and would swallow a sentinel. g() is the same
                // residual, so the two arms stop on the same condition.
                odeSolver.clearEventHandlers();
                if (fpArmed) {
                    odeSolver.addEventHandler(
                            new FluidSettledEvent(ode, initialState.length, minNonZeroRate),
                            FastMath.max(tRange[1] - tRange[0], GlobalConstants.Zero),
                            FastMath.max(GlobalConstants.Zero * FastMath.max(1.0, FastMath.abs(tRange[1])),
                                    Double.MIN_NORMAL),
                            100);
                }

                boolean usedStiffFallback = false;
                try {
                    try {
                        odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                    } catch (RuntimeException e) {
                        if (e.getMessage() != null && e.getMessage().contains("step size")) {
                            usedStiffFallback = true;
                        } else {
                            odeSolver.clearStepHandlers();
                            if (stepHandler != null) {
                                odeSolver.addStepHandler(stepHandler);
                            }
                            // Message deferred until the retry succeeds: a failure here falls
                            // through to the stiff fallback below, which means the initial
                            // point was NOT what went wrong.
                            odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState);
                            if (options.verbose != VerboseLevel.SILENT) {
                                System.out.println("The initial point was invalid, Fluid solver switched to default initialization.");
                            }
                        }
                    }
                } catch (RuntimeException e) {
                    usedStiffFallback = true;
                }

                if (usedStiffFallback) {
                    LSODA stiffSolver = options.odesolvers.stiffIntegratorFor(tRange[0], tRange[1], options.tol,
                            options.tol > GlobalConstants.CoarseTol);
                    // This arm reaches LSODA too, so it carries the same
                    // per-step stop; without it a window that settles here
                    // stalls exactly as the main stiff branch would.
                    stiffSolver.setStopCondition(fpStop);
                    stiffSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                    java.util.ArrayList<Double> tHistory = stiffSolver.getTvec();
                    java.util.ArrayList<Double[]> yHistory = stiffSolver.getYvec();
                    int dim = ode.getDimension();
                    Tmax = stiffSolver.getStepsTaken() + 1;
                    int lastIdx = Tmax - 1;
                    this.xvec_it = new Matrix(1, dim);
                    for (int j = 0; j < dim; j++) this.xvec_it.set(0, j, Math.max(0.0, yHistory.get(lastIdx)[j]));

                    if (keepTrajectory) {
                        DMatrixRMaj denseT = new DMatrixRMaj(Tmax, 1);
                        DMatrixRMaj denseX = new DMatrixRMaj(Tmax, dim);
                        for (int i = 0; i < Tmax; i++) {
                            denseT.set(i, 0, tHistory.get(i));
                            for (int j = 0; j < dim; j++) denseX.set(i, j, Math.max(0.0, yHistory.get(i)[j]));
                        }
                        tIterations.add(new Matrix(denseT));
                        xVecIterations.add(new Matrix(denseX));
                    }
                } else if (stepHandler != null) {
                    tIterations.add(stepHandler.tVec);
                    xVecIterations.add(stepHandler.xVec);
                    Tmax = stepHandler.tVec.getNumRows();
                    this.xvec_it = Matrix.extractRows(stepHandler.xVec, Tmax - 1, Tmax, null);
                } else {
                    // no trajectory recorded: the integrator's terminal state is the fixed point
                    Tmax = 1;
                    this.xvec_it = new Matrix(1, nextState.length);
                    for (int j = 0; j < nextState.length; j++) {
                        this.xvec_it.set(0, j, Math.max(0.0, nextState[j]));
                    }
                }
            }
            totalSteps += Tmax;

            // The window finished; did it finish ON the model? A divergence that
            // does not trip the integrator leaves a state that is not a solution,
            // and every later window is then integrated from it.
            if (guard != null) {
                double[] endState = new double[this.xvec_it.getNumCols()];
                for (int j = 0; j < endState.length; j++) {
                    endState[j] = this.xvec_it.get(0, j);
                }
                guard.assertConserved(endState, T);
            }

            jline.io.LineConsole.iter(iter, "window %d up to t = %g: %d ODE steps", iter, T, Tmax);

            T0 = T;

            // THE TERMINATION TEST. The mass moved by this pass measures movement over
            // ONE window, and underestimates the distance still left to the fixed point
            // by exactly the geometric tail it drops: for a mode contracting by rho per
            // window that tail is r*rho/(1-rho), which is what is tested. rho is read
            // off the iteration itself, so no rate in the model has to stand in for the
            // slowest system mode -- when one does, as a bare drift norm must, the stop
            // lands 3% short. The drift, zero AT a fixed point, is an independent second
            // bound. Both must hold on two consecutive windows, past the relaxation time.
            double moved = 0.0, totalPrev = 0.0;
            for (int j = 0; j < initialState.length; j++) {
                moved += FastMath.abs(this.xvec_it.get(0, j) - initialState[j]);
                totalPrev += initialState[j];
            }
            double ratio = (totalPrev > 0.0) ? moved / 2.0 / totalPrev : 0.0;
            if (earlystop && iter > 1 && !Double.isFinite(options.timespan[1]) && T >= minHorizon) {
                rhoHist[iter % rhoHist.length] = ratio / FastMath.max(movedPrev, GlobalConstants.Zero);
                double rho = 0.0;
                for (double v : rhoHist) {
                    if (!Double.isNaN(v) && !Double.isInfinite(v) && v > rho) {
                        rho = v;
                    }
                }
                double[] yEnd = new double[initialState.length];
                double[] dy = new double[initialState.length];
                double total = 0.0;
                for (int j = 0; j < yEnd.length; j++) {
                    yEnd[j] = this.xvec_it.get(0, j);
                    total += yEnd[j];
                }
                ode.computeDerivatives(T, yEnd, dy);
                double dn = 0.0;
                for (double v : dy) {
                    dn += FastMath.abs(v);
                }
                double driftDispl = (total > 0.0) ? dn / 2.0 / total / minNonZeroRate : 0.0;
                // a non-contracting iteration has no tail to sum: it is not converging
                //
                // THE 1e-6 GATE IS NOT AN OVERSIGHT, even though it sits below the
                // integrator's own tol. Relaxing it to "the moved mass reached the
                // integrator floor, so trust the drift residual alone" was TRIED and
                // REVERTED: it stops the M/M/1 rho = 0.9 minnormal solve at 7.018088
                // against the 7.021524680 all four codebases agree on
                // (MinNormalTest.testOpenMm1), and it ends the statedep trajectory of
                // test_exampleCdfRespT2StatedepMethod at t = 60 instead of its 2000.
                // The accuracy of this loop comes from running the windows, so the stop
                // has to stay conservative. It is also not what makes a solve hang: see
                // _kb/06-solver-catalog.md, where the minnormal closure diverges
                // outright on a bounded multiserver station.
                if (rho < 1.0 && ratio * rho / (1.0 - rho) < driftSafety * driftTol
                        && driftDispl < driftTol) {
                    driftBelow++;
                    if (driftBelow >= 2) {
                        goon = false;
                    }
                } else {
                    driftBelow = 0;
                }
            }
            movedPrev = ratio;

            if (T >= options.timespan[1]) {
                goon = false;
            }
        }

        if (!keepTrajectory) {
            // QNt/UNt/TNt are read at row 0 (initial condition) and at the last row (fixed
            // point), so the two endpoints carry the whole contract of an infinite timespan
            int dim = this.xvec_it.getNumCols();
            this.xvec_t = new Matrix(2, dim);
            for (int j = 0; j < dim; j++) {
                this.xvec_t.set(0, j, FastMath.max(0.0, xvecInit.get(0, j)));
                this.xvec_t.set(1, j, this.xvec_it.get(0, j));
            }
            result.t = new Matrix(2, 1);
            result.t.set(0, 0, options.timespan[0]);
            result.t.set(1, 0, T);
        } else if (!xVecIterations.isEmpty() && totalSteps > 0) {
            int nextRow = 0;
            int cols = xVecIterations.get(0).getNumCols();
            DMatrixRMaj denseXvecT = new DMatrixRMaj(totalSteps, cols);
            DMatrixRMaj denseTvec = new DMatrixRMaj(totalSteps, 1);
            for (int i = 0; i < xVecIterations.size(); i++) {
                Matrix tIter = tIterations.get(i);
                Matrix xVecIter = xVecIterations.get(i);
                int stepsPerIter = tIter.getNumRows();
                for (int j = nextRow; j < nextRow + stepsPerIter; j++) {
                    denseTvec.set(j, 0, tIter.get(j - nextRow, 0));
                    for (int k = 0; k < cols; k++) {
                        denseXvecT.set(j, k, xVecIter.get(j - nextRow, k));
                    }
                }
                nextRow += stepsPerIter;
            }
            this.xvec_t = new Matrix(DConvertMatrixStruct.convert(denseXvecT, (DMatrixSparseCSC) null, 0.0));
            result.t = new Matrix(DConvertMatrixStruct.convert(denseTvec, (DMatrixSparseCSC) null, 0.0));
        }
    }

    /**
     * The caller's requested output instants that fall strictly inside a window.
     *
     * @param pts requested instants, increasing, or null
     * @param t0  window start
     * @param t1  window end
     * @return the instants in (t0, t1), or null when there are none
     */
    /**
     * The mid-window fixed-point stop for the Commons Math arm.
     *
     * <p>{@link #g(double, double[])} is the SAME normalized residual the LSODA
     * hook and the entry test read, offset by {@code GlobalConstants.Zero} so
     * the root the integrator brackets is the instant the drift becomes zero to
     * double precision. It starts positive -- a window whose entry state was
     * already settled never reaches an integrator -- and crossing DOWN through
     * it ends the window.</p>
     */
    private static final class FluidSettledEvent implements EventHandler {
        private final FirstOrderDifferentialEquations ode;
        private final int dim;
        private final double minNonZeroRate;

        FluidSettledEvent(FirstOrderDifferentialEquations ode, int dim, double minNonZeroRate) {
            this.ode = ode;
            this.dim = dim;
            this.minNonZeroRate = minNonZeroRate;
        }

        @Override
        public void init(double t0, double[] y0, double t) {
        }

        @Override
        public double g(double t, double[] y) {
            if (dim <= 0 || !(minNonZeroRate > 0.0)) {
                return 1.0;
            }
            double[] dy = new double[dim];
            ode.computeDerivatives(t, y, dy);
            double residual = 0.0, total = 0.0;
            for (int j = 0; j < dim; j++) {
                residual += FastMath.abs(dy[j]);
                total += y[j];
            }
            if (!(total > 0.0)) {
                return 1.0;
            }
            return residual / 2.0 / total / minNonZeroRate - GlobalConstants.Zero;
        }

        @Override
        public Action eventOccurred(double t, double[] y, boolean increasing) {
            return Action.STOP;
        }

        @Override
        public void resetState(double t, double[] y) {
        }
    }

    /**
     * True when the drift at {@code y} is zero to double precision, i.e. the
     * state IS a fixed point rather than merely near one.
     *
     * <p>THE THRESHOLD IS ROUND-OFF, NOT THE SOLVER TOLERANCE, and that is the
     * whole safety of it. {@code driftTol} (1e-4 by default) says "converged to
     * what the caller asked for", and a state that merely satisfies THAT is
     * still moving -- cutting a window there was measured to shift results by
     * 1.8e-5 in the Python twin. A normalized residual below
     * {@code GlobalConstants.Zero} is the stronger claim, which is what makes
     * holding the state for the rest of the span exact rather than
     * approximate.</p>
     *
     * <p>The residual is normalized by the population AND by the slowest exit
     * rate, so the test is scale-free in time as well as in mass, and is the
     * same quantity the window loop's drift displacement reads.</p>
     *
     * @param ode   the drift
     * @param t     the instant to evaluate at
     * @param y     the state, laid out from {@code offset}
     * @param dim   number of coordinates
     * @param offset 0 for a plain array, 1 for LSODA's 1-indexed state
     * @param minNonZeroRate the slowest exit rate
     * @return true when the window can be ended in closed form
     */
    private static boolean fluidSettled(FirstOrderDifferentialEquations ode, double t,
                                        double[] y, int dim, int offset,
                                        double minNonZeroRate) {
        if (dim <= 0 || !(minNonZeroRate > 0.0)) {
            return false;
        }
        double[] state = new double[dim];
        System.arraycopy(y, offset, state, 0, dim);
        double[] dy = new double[dim];
        ode.computeDerivatives(t, state, dy);
        double residual = 0.0, total = 0.0;
        for (int j = 0; j < dim; j++) {
            residual += FastMath.abs(dy[j]);
            total += state[j];
        }
        return total > 0.0 && residual / 2.0 / total / minNonZeroRate < GlobalConstants.Zero;
    }

    private static double[] windowOutputTimes(double[] pts, double t0, double t1) {
        if (pts == null || pts.length == 0) {
            return null;
        }
        int count = 0;
        for (int i = 0; i < pts.length; i++) {
            if (pts[i] > t0 && pts[i] < t1) {
                count++;
            }
        }
        if (count == 0) {
            return null;
        }
        double[] out = new double[count];
        int j = 0;
        for (int i = 0; i < pts.length; i++) {
            if (pts[i] > t0 && pts[i] < t1) {
                out[j++] = pts[i];
            }
        }
        return out;
    }

    /**
     * Integrates the closing ODE to its fixed point and fills the queue-length
     * measures, leaving {@link #xvec_it} and {@link #xvec_t} at the mean
     * trajectory. Protected so the moment-closure driver can drive the same mean
     * solve under a supplied closure variance, exactly as MATLAB
     * {@code solver_fluid_moments} calls {@code solver_fluid}.
     */
    protected void solver_fluid(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix S = sn.nservers.copy();

        Map<Station, Map<JobClass, Matrix>> mu = new HashMap<Station, Map<JobClass, Matrix>>();
        Map<Station, Map<JobClass, Matrix>> phi = new HashMap<Station, Map<JobClass, Matrix>>();
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            Map<JobClass, Matrix> muCopy = new HashMap<JobClass, Matrix>();
            Map<JobClass, Matrix> phiCopy = new HashMap<JobClass, Matrix>();
            for (int k = 0; k < K; k++) {
                JobClass jobClass = sn.jobclasses.get(k);
                muCopy.put(jobClass, sn.mu.get(station).get(jobClass).copy());
                phiCopy.put(jobClass, sn.phi.get(station).get(jobClass).copy());
            }
            mu.put(station, muCopy);
            phi.put(station, phiCopy);
        }

        Matrix match = new Matrix(M, K);
        Matrix phases = new Matrix(M, K);
        Matrix slowrate = new Matrix(M, K);
        // station-major routing, hoisted: the reduction inverts a matrix
        Matrix rtSt = SnRtStations.snRtStations(sn).getLeft();
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int k = 0; k < K; k++) {
                JobClass jobClass = sn.jobclasses.get(k);

                if (mu.get(station).get(jobClass).hasNaN()) {
                    mu.get(station).put(jobClass, new Matrix(0, 0));
                    phi.get(station).put(jobClass, new Matrix(0, 0));
                }

                if (rtSt.sumCols((i * K) + k) > 0) {
                    match.set(i, k, 1);
                }

                if (Double.isInfinite(S.get(i, 0))) {
                    S.set(i, 0, sn.nclosedjobs);
                }

                phases.set(i, k, mu.get(station).get(jobClass).length());

                if (mu.get(station).get(jobClass).isEmpty()) {
                    slowrate.set(i, k, Double.POSITIVE_INFINITY);
                } else {
                    slowrate.set(i, k, mu.get(station).get(jobClass).elementMin());
                }
            }
        }

        Matrix y0 = new Matrix(1, 0);
        Matrix assigned = new Matrix(1, K);
        double toAssign;
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (match.get(i, k) > 0 && phases.get(i, k) > 0) {
                    if (Double.isInfinite(sn.njobs.get(0, k))) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                            toAssign = 1.0;
                        } else {
                            toAssign = 0.0;
                        }
                    } else {
                        toAssign = FastMath.floor(sn.njobs.get(0, k) / match.sumCols(k));
                        if (match.sumSubMatrix(i + 1, match.getNumRows(), k, k + 1) == 0.0) {
                            toAssign = sn.njobs.get(0, k) - assigned.get(0, k);
                        }
                    }

                    int originalY0Length = y0.getNumCols();
                    int newY0Length = originalY0Length + 1 + (int) (phases.get(i, k) - 1);
                    y0.expandMatrix(1, newY0Length, newY0Length);
                    y0.set(0, originalY0Length, toAssign);
                    int col = 0;
                    while (col < phases.get(i, k) - 1) {
                        y0.set(0, originalY0Length + 1 + col, 0);
                        col++;
                    }
                    assigned.set(0, k, assigned.get(0, k) + toAssign);
                } else {
                    int originalY0Length = y0.getNumCols();
                    int newY0Length = originalY0Length + (int) phases.get(i, k);
                    y0.expandMatrix(1, newY0Length, newY0Length);
                    int col = 0;
                    while (col < phases.get(i, k)) {
                        y0.set(0, originalY0Length + col, 0);
                        col++;
                    }
                }
            }
        }

        int y0cols = y0.getNumCols();
        double[] yDefault = new double[y0cols];
        if (options.init_sol.isEmpty()) {
            xvec_it = y0;
            for (int i = 0; i < y0cols; i++) {
                yDefault[i] = y0.get(0, i);
            }
        } else {
            xvec_it = options.init_sol;
            for (int i = 0; i < y0cols; i++) {
                yDefault[i] = options.init_sol.get(0, i);
            }
        }

        solver_fluid_iteration(sn, mu, phi, S, yDefault, slowrate, options, result);

        // PS, DPS, GPS scheduling
        result.QN = new Matrix(M, K);
        result.QNt = new Matrix[M][K];
        Matrix[] Qt = new Matrix[M];
        result.UNt = new Matrix[M][K];
        int Tmax = xvec_t.getNumRows();

        for (int i = 0; i < M; i++) {
            Qt[i] = new Matrix(Tmax, 1);
            for (int k = 0; k < K; k++) {
                int shift = (int) phases.sumSubMatrix(0, i, 0, phases.getNumCols())
                        + (int) phases.sumSubMatrix(i, i + 1, 0, k);
                result.QN.set(i, k, xvec_it.sumSubMatrix(0, 1, shift, shift + (int) phases.get(i, k)));

                result.QNt[i][k] = new Matrix(Tmax, 1);
                for (int step = 0; step < Tmax; step++) {
                    result.QNt[i][k].set(step, 0,
                            xvec_t.sumSubMatrix(step, step + 1, shift, shift + (int) phases.get(i, k)));
                }

                Qt[i] = Qt[i].add(1.0, result.QNt[i][k]);
            }
            // A Source holds no jobs: its coordinate is the arrival process's
            // normalisation constant, not a population (the same reason
            // SolverFluid.getProbAggr excludes it), and the station-level
            // routing folds every departure of an open class back onto it, so
            // from an initial state above the fixed point it accumulates the
            // mass that LEFT the system. Reporting that as a queue length also
            // feeds it back as jobs at the Source whenever a caller carries the
            // marginal across (SolverENV).
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                for (int k = 0; k < K; k++) {
                    result.QN.set(i, k, 0.0);
                    result.QNt[i][k] = new Matrix(Tmax, 1);
                }
                Qt[i] = new Matrix(Tmax, 1);
            }
        }

        for (int i = 0; i < M; i++) {
            if (sn.nservers.get(i, 0) > 0) {
                for (int k = 0; k < K; k++) {
                    result.UNt[i][k] = new Matrix(Tmax, 1);
                    for (int step = 0; step < Tmax; step++) {
                        result.UNt[i][k].set(step, 0,
                                FastMath.min(result.QNt[i][k].get(step, 0) / S.get(i, 0),
                                        result.QNt[i][k].get(step, 0) / Qt[i].get(step, 0)));
                        if (Double.isNaN(result.UNt[i][k].get(step, 0))) {
                            result.UNt[i][k].set(step, 0, 0);
                        }
                    }
                }
            } else {
                System.arraycopy(result.QNt[i], 0, result.UNt[i], 0, K);
            }
        }
    }

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Map<Station, Map<JobClass, Matrix>> lambda = sn.mu;

        solver_fluid(sn, options, result);

        Matrix delayNodes = new Matrix(1, M);
        for (int i = 0; i < M; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                delayNodes.set(0, i, 1);
            }
        }

        result.TN = new Matrix(M, K);
        result.TNt = new Matrix[M][K];
        int Tmax = result.t.getNumRows();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                result.TNt[i][k] = new Matrix(Tmax, 1);
            }
        }

        Matrix[][] Xservice = new Matrix[M][K];

        for (int i = 0; i < M; i++) {
            if (delayNodes.get(0, i) == 1.0) {
                for (int k = 0; k < K; k++) {
                    int idx = (int) sn.phases.sumSubMatrix(0, i, 0, K)
                            + (int) sn.phases.sumSubMatrix(i, i + 1, 0, k);
                    Xservice[i][k] = new Matrix((int) sn.phases.get(i, k), 1);
                    for (int f = 0; f < (int) sn.phases.get(i, k); f++) {
                        double lambdaIKF = lambda.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(f, 0);
                        double phiIKF = sn.phi.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(f, 0);

                        result.TN.set(i, k, result.TN.get(i, k) + (xvec_it.get(0, idx + f) * lambdaIKF * phiIKF));

                        Matrix tmpTNt = result.QNt[i][k].copy();
                        tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt);
                        result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt);

                        Xservice[i][k].set(f, 0, (xvec_it.get(0, idx + f) * lambdaIKF));
                    }
                }
            } else {
                double xi = result.QN.sumRows(i);
                Matrix xi_t = result.QNt[i][0].copy();
                for (int r = 1; r < K; r++) {
                    xi_t = xi_t.add(1.0, result.QNt[i][r]);
                }

                double wni = GlobalConstants.FineTol;
                Matrix wi = new Matrix(1, K);
                if (xi > 0 || sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                        for (int k = 0; k < K; k++) {
                            int idx = (int) sn.phases.sumSubMatrix(0, i, 0, K)
                                    + (int) sn.phases.sumSubMatrix(i, i + 1, 0, k);
                            Matrix D0 = sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(0);
                            Matrix D1 = sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(1);
                            wi.set(0, k, Map_mean.map_mean(D0, D1));

                            wni += wi.get(0, k) * xvec_it.sumSubMatrix(0, 1, idx, idx + (int) sn.phases.get(i, k));
                        }
                    }

                    for (int k = 0; k < K; k++) {
                        int idx = (int) sn.phases.sumSubMatrix(0, i, 0, K)
                                + (int) sn.phases.sumSubMatrix(i, i + 1, 0, k);
                        Xservice[i][k] = new Matrix((int) sn.phases.get(i, k), 1);

                        for (int f = 0; f < (int) sn.phases.get(i, k); f++) {
                            double lambdaIKF = lambda.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(f, 0);
                            double phiIKF = sn.phi.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(f, 0);

                            SchedStrategy ss = sn.sched.get(sn.stations.get(i));
                            if (ss == SchedStrategy.EXT) {
                                if (f == 0) {
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + (lambdaIKF * phiIKF * (1
                                                    - xvec_it.sumSubMatrix(0, 1, idx + 1,
                                                            idx + (int) sn.phases.get(i, k)))));

                                    Matrix tmpTNt = xvec_t.sumRows(idx + 1, idx + (int) sn.phases.get(i, k));
                                    Matrix ones = new Matrix(tmpTNt.getNumRows(), tmpTNt.getNumCols());
                                    ones.ones();
                                    tmpTNt = ones.sub(tmpTNt);
                                    tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt);
                                    result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt);

                                    Xservice[i][k].set(f, 0,
                                            lambdaIKF * (1 - xvec_it.sumSubMatrix(0, 1, idx + 1,
                                                    idx + (int) sn.phases.get(i, k))));
                                } else {
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + (lambdaIKF * phiIKF * xvec_it.get(0, idx + f)));
                                    Matrix tmpTNt = xvec_t.sumRows(idx + f, idx + f + 1);
                                    tmpTNt.scaleEq(lambdaIKF * phiIKF, tmpTNt);
                                    result.TNt[i][k] = result.TNt[i][k].add(1.0, tmpTNt);
                                    Xservice[i][k].set(f, 0, lambdaIKF * xvec_it.get(0, idx + f));
                                }
                            } else if (ss == SchedStrategy.INF || ss == SchedStrategy.PS) {
                                result.TN.set(i, k, result.TN.get(i, k)
                                        + ((lambdaIKF * phiIKF / xi) * xvec_it.get(0, idx + f)
                                                * FastMath.min(xi, sn.nservers.get(i, 0))));

                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / xi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                Xservice[i][k].set(f, 0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f));
                            } else if (ss == SchedStrategy.DPS) {
                                Matrix w = new Matrix(1, K);
                                int p = 0;
                                while (p < K) {
                                    w.set(0, p, sn.schedparam.get(i, p));
                                    p++;
                                }

                                Matrix tmpQ = new Matrix(1, result.QN.getNumCols());
                                Matrix.extractRows(result.QN, i, i + 1, tmpQ);
                                tmpQ = tmpQ.transpose();
                                Matrix wxi = new Matrix(1, 1);
                                wxi = w.mult(tmpQ, wxi);

                                Matrix wxi_t = result.QNt[i][0].copy();
                                wxi_t.scaleEq(w.value(), wxi_t);
                                int r = 1;
                                while (r < result.QNt[i][0].getNumCols()) {
                                    Matrix tmpWxi_t = result.QNt[i][r].copy();
                                    tmpWxi_t.scaleEq(w.get(0, r), tmpWxi_t);
                                    wxi_t = wxi_t.add(1.0, tmpWxi_t);
                                    r++;
                                }

                                result.TN.set(i, k, result.TN.get(i, k)
                                        + (((lambdaIKF * phiIKF * w.get(0, k)) / wxi.value())
                                                * xvec_it.get(0, idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))));

                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF * w.get(0, k), tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / wxi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                Xservice[i][k].set(f, 0,
                                        ((lambdaIKF * w.get(0, k) / wxi.value())
                                                * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f)));
                            } else if (ss == SchedStrategy.FCFS || ss == SchedStrategy.SIRO || ss == SchedStrategy.LCFS) {
                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / xi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                if ("default".equals(options.method) || "closing".equals(options.method)
                                        || "softmin".equals(options.method)) {
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + ((lambdaIKF * phiIKF / xi) * xvec_it.get(0, idx + f)
                                                    * FastMath.min(xi, sn.nservers.get(i, 0))));
                                    Xservice[i][k].set(f, 0,
                                            lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f));
                                    break;
                                } else if ("statedep".equals(options.method)) {
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + (((lambdaIKF * phiIKF * wi.get(0, k)) / wni)
                                                    * xvec_it.get(0, idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))));

                                    Xservice[i][k].set(f, 0,
                                            ((lambdaIKF * wi.get(0, k) / wni)
                                                    * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f)));
                                }
                            } else if (ss == SchedStrategy.HOL || ss == SchedStrategy.FCFSPRIO || ss == SchedStrategy.SJF
                                    || ss == SchedStrategy.LJF || ss == SchedStrategy.SEPT || ss == SchedStrategy.LEPT) {
                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / xi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                if ("default".equals(options.method) || "closing".equals(options.method)
                                        || "softmin".equals(options.method)) {
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + ((lambdaIKF * phiIKF / xi) * xvec_it.get(0, idx + f)
                                                    * FastMath.min(xi, sn.nservers.get(i, 0))));
                                    Xservice[i][k].set(f, 0,
                                            lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f));
                                } else if ("statedep".equals(options.method)) {
                                    double priority = sn.schedparam.get(i, k);
                                    result.TN.set(i, k, result.TN.get(i, k)
                                            + (((lambdaIKF * phiIKF * priority) / wni)
                                                    * xvec_it.get(0, idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))));

                                    Xservice[i][k].set(f, 0,
                                            ((lambdaIKF * priority / wni)
                                                    * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f)));
                                }
                            } else if (ss == SchedStrategy.GPS || ss == SchedStrategy.PSPRIO
                                    || ss == SchedStrategy.DPSPRIO || ss == SchedStrategy.GPSPRIO) {
                                Matrix w = new Matrix(1, K);
                                int p = 0;
                                while (p < K) {
                                    w.set(0, k, sn.schedparam.get(i, k));
                                    p++;
                                }

                                Matrix tmpQ = new Matrix(1, result.QN.getNumCols());
                                Matrix.extractRows(result.QN, i, i + 1, tmpQ);
                                tmpQ = tmpQ.transpose();
                                Matrix wxi = new Matrix(1, 1);
                                wxi = w.mult(tmpQ, wxi);

                                result.TN.set(i, k, result.TN.get(i, k)
                                        + (((lambdaIKF * phiIKF * w.get(0, k)) / wxi.value())
                                                * xvec_it.get(0, idx + f) * FastMath.min(xi, sn.nservers.get(i, 0))));

                                Xservice[i][k].set(f, 0,
                                        ((lambdaIKF * w.get(0, k) / wxi.value())
                                                * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f)));
                            } else if (ss == SchedStrategy.LCFSPR) {
                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / xi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                result.TN.set(i, k, result.TN.get(i, k)
                                        + ((lambdaIKF * phiIKF / xi) * xvec_it.get(0, idx + f)
                                                * FastMath.min(xi, sn.nservers.get(i, 0))));
                                Xservice[i][k].set(f, 0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f));
                            } else if (ss == SchedStrategy.POLLING || ss == SchedStrategy.FORK || ss == SchedStrategy.REF) {
                                Matrix tmpTNT = new Matrix(xvec_t.getNumRows(), 1);
                                Matrix.extractColumn(xvec_t, idx + f, tmpTNT);
                                tmpTNT.scaleEq(lambdaIKF * phiIKF, tmpTNT);
                                int rows = tmpTNT.getNumRows();
                                int row = 0;
                                while (row < rows) {
                                    result.TNt[i][k].set(row, 0,
                                            result.TNt[i][k].get(row, 0)
                                                    + (tmpTNT.get(row, 0) / xi_t.get(row, 0)
                                                            * FastMath.min(xi_t.get(row, 0), sn.nservers.get(i, 0))));
                                    row++;
                                }

                                result.TN.set(i, k, result.TN.get(i, k)
                                        + ((lambdaIKF * phiIKF / xi) * xvec_it.get(0, idx + f)
                                                * FastMath.min(xi, sn.nservers.get(i, 0))));
                                Xservice[i][k].set(f, 0,
                                        lambdaIKF / xi * FastMath.min(xi, sn.nservers.get(i, 0)) * xvec_it.get(0, idx + f));
                            } else {
                                throw new RuntimeException("Unsupported scheduling policy: " + sn.sched.get(sn.stations.get(i)));
                            }
                        }
                    }
                }
            }
        }

        // THE COMPLETIONS AN ELIMINATED COORDINATE MAKES ARE NOT LOST. The loop above
        // reads throughputs off the STATE, as x_f * mu_f * phi_f summed over phases, and an
        // eliminated coordinate holds no mass there -- so its completions, which are FINITE
        // because mu_f is InfRate, would silently vanish and the station would stop balancing
        // against its neighbours. Their total rate is exactly what Emap carries: the composed
        // event that replaced the inflow stands for the original completion too.
        addImmediateThroughput(this.lastPassageOde, sn, result, M, K);

        // Response Times
        result.RN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                // See SolverFluid: TN is zero only to the integrator's accuracy.
                if (result.TN.get(i, k) > GlobalConstants.Zero) {
                    result.RN.set(i, k, result.QN.get(i, k) / result.TN.get(i, k));
                }
            }
        }

        // Utilisation
        result.UN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                double sumForUN = 0.0;
                double sumForTN = 0.0;
                if (Xservice[i] != null && Xservice[i][k] != null) {
                    int rows = Xservice[i][k].getNumRows();
                    for (int row = 0; row < rows; row++) {
                        if (Xservice[i][k].get(row, 0) > 0) {
                            sumForUN += (Xservice[i][k].get(row, 0)
                                    / lambda.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(row, 0));
                            sumForTN += Xservice[i][k].get(row, 0);
                        }
                    }
                }
                result.UN.set(i, k, sumForUN);
                if ((sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) && ("statedep".equals(options.method))) {
                    result.TN.set(i, k, sumForTN);
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (delayNodes.get(0, i) == 0.0) {
                    result.UN.set(i, k, result.UN.get(i, k) / sn.nservers.get(i, 0));
                }
            }
        }
        result.WN = new Matrix(0, 0);
        result.AN = new Matrix(0, 0);
    }

    @Override
    public Matrix getXVecIt() {
        return this.xvec_it;
    }


    /**
     * Add to {@code result.TN} the service completions made by coordinates that the immediate
     * elimination removed from the event set. A no-op when nothing was eliminated.
     *
     * @param ode    the ODE object carrying Emap and the original event identity
     * @param sn     the network structure
     * @param result the result whose TN is corrected in place
     * @param M      station count
     * @param K      class count
     */
    private void addImmediateThroughput(PassageTimeODE ode, NetworkStruct sn,
                                        SolverResult result, int M, int K) {
        if (ode == null) {
            return;
        }
        Matrix emap = ode.getImmediateEmap();
        Matrix absorb = ode.getImmediateAbsorb();
        Matrix origEventIdx = ode.getOriginalEventIdx();
        if (emap == null || absorb == null || origEventIdx == null) {
            return;
        }
        Matrix qIndices = ode.getQIndices();
        Matrix kic = ode.getKic();
        int nstate = (int) absorb.getNumRows();
        int[] coordStation = new int[nstate];
        int[] coordClass = new int[nstate];
        for (int i = 0; i < nstate; i++) {
            coordStation[i] = -1;
            coordClass[i] = -1;
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int lo = (int) qIndices.get(i, k);
                for (int f = 0; f < (int) kic.get(i, k); f++) {
                    if (lo + f < nstate) {
                        coordStation[lo + f] = i;
                        coordClass[lo + f] = k;
                    }
                }
            }
        }
        // A coordinate was eliminated exactly when its row of the projector is no longer the
        // identity row: the reduction replaced it by its absorption distribution.
        boolean[] eliminated = new boolean[nstate];
        boolean any = false;
        for (int f = 0; f < nstate; f++) {
            if (absorb.get(f, f) != 1.0) {
                eliminated[f] = true;
                any = true;
            }
        }
        if (!any) {
            return;
        }
        double[] x = new double[(int) this.xvec_it.getNumCols()];
        for (int j = 0; j < x.length; j++) {
            x[j] = this.xvec_it.get(0, j);
        }
        Matrix r = ode.calculateRatesClosing(x);
        int nDep = Math.min(ode.getOriginalDepartureCount(), (int) origEventIdx.getNumRows());
        for (int o = 0; o < nDep; o++) {
            int c = (int) origEventIdx.get(o, 0);
            if (c < 0 || c >= nstate || !eliminated[c]) {
                continue;
            }
            int i = coordStation[c];
            int k = coordClass[c];
            if (i < 0 || k < 0) {
                continue;
            }
            double extra = 0;
            for (int e = 0; e < emap.getNumRows(); e++) {
                double w = emap.get(e, o);
                if (w != 0) {
                    extra += w * r.get(e, 0);
                }
            }
            result.TN.set(i, k, result.TN.get(i, k) + extra);
        }
    }

}
