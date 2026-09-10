/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mc.Dtmc_stochcomp;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.handlers.FluidHideImmediate;
import jline.solvers.fluid.handlers.MatrixMethodODE;
import jline.solvers.fluid.handlers.TransientDataHandler;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.matrix.MatrixEquation;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrixRMaj;

import static jline.io.InputOutput.line_warning;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;

public class MatrixMethodAnalyzer implements FluidAnalyzer {
    public Matrix xvec_t;
    public Matrix xvec_it;

    /**
     * Whether (station ist, class r) carries a non-homogeneous Poisson process.
     * An NHPP slot of sn.proc holds the rate schedule {breakpoints, rates,
     * cyclic}, not a {D0, D1} MAP, so it must never be read as phase-type
     * matrices. Mirrors the MATLAB solver_fluid_matrix substitution.
     */
    private static boolean isNhpp(NetworkStruct sn, Station station, JobClass jobClass) {
        if (sn.procid == null) {
            return false;
        }
        java.util.Map<JobClass, ProcessType> stProcid = sn.procid.get(station);
        if (stProcid == null) {
            return false;
        }
        return stProcid.get(jobClass) == ProcessType.NHPP;
    }

    /**
     * The process representation the matrix method integrates. The matrix
     * method is a steady-state method, so an NHPP enters at its nominal
     * (time-average) rate as a single-phase exponential MAP {-lam, lam}; the
     * time-varying intensity is a transient property carried by the closing
     * method's per-event rate multiplier.
     */
    private static MatrixCell effProc(NetworkStruct sn, Station station, JobClass jobClass,
                                      int ist, int r) {
        if (isNhpp(sn, station, jobClass)) {
            double lam = sn.rates.get(ist, r);
            if (!Double.isNaN(lam) && !Double.isInfinite(lam) && lam > 0) {
                Matrix d0 = new Matrix(1, 1, 1);
                d0.set(0, 0, -lam);
                Matrix d1 = new Matrix(1, 1, 1);
                d1.set(0, 0, lam);
                MatrixCell cell = new MatrixCell(2);
                cell.set(0, d0);
                cell.set(1, d1);
                return cell;
            }
        }
        return sn.proc.get(station).get(jobClass);
    }

    /** The entry-phase distribution matching {@link #effProc}. */
    private static Matrix effPie(NetworkStruct sn, Station station, JobClass jobClass,
                                 int ist, int r) {
        if (isNhpp(sn, station, jobClass)) {
            double lam = sn.rates.get(ist, r);
            if (!Double.isNaN(lam) && !Double.isInfinite(lam) && lam > 0) {
                Matrix pie = new Matrix(1, 1, 1);
                pie.set(0, 0, 1.0);
                return pie;
            }
        }
        return sn.pie.get(station).get(jobClass);
    }

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;

        Matrix S = sn.nservers.copy();
        double initialPopulation = sn.njobs.elementSum();
        int SRows = S.getNumRows();
        for (int i = 0; i < SRows; i++) {
            if (Double.isInfinite(S.get(i, 0))) {
                S.set(i, 0, initialPopulation);
            }
        }

        List<Integer> stationIndices = new ArrayList<Integer>();
        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            for (int r = 0; r < K; r++) {
                stationIndices.add(isf * K + r);
            }
        }
        Matrix P = Dtmc_stochcomp.dtmc_stochcomp(sn.rt, stationIndices);

        for (int srcIst = 0; srcIst < M; srcIst++) {
            if (sn.sched.get(sn.stations.get(srcIst)) == SchedStrategy.EXT) {
                for (int r = 0; r < K; r++) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        int srcCol = srcIst * K + r;
                        for (int fromIst = 0; fromIst < M; fromIst++) {
                            if (fromIst != srcIst) {
                                for (int fromR = 0; fromR < K; fromR++) {
                                    P.set(fromIst * K + fromR, srcCol, 0.0);
                                }
                            }
                        }
                    }
                }
            }
        }

        Matrix psi = new Matrix(0, 0);
        Matrix A = new Matrix(0, 0);
        Matrix B = new Matrix(0, 0);

        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int r = 0; r < K; r++) {
                JobClass jobClass = sn.jobclasses.get(r);
                if (sn.phases.get(i, r) == 0.0) {
                    Matrix zeroMatrix = new Matrix(1, 1, 1);
                    Matrix nanMatrix = new Matrix(1, 1, 1);
                    nanMatrix.set(0, 0, Double.NaN);
                    psi = psi.createBlockDiagonal(zeroMatrix);
                    A = A.createBlockDiagonal(nanMatrix);
                    B = B.createBlockDiagonal(zeroMatrix);
                } else {
                    psi = psi.createBlockDiagonal(effProc(sn, station, jobClass, i, r).get(0));
                    A = A.createBlockDiagonal(effPie(sn, station, jobClass, i, r).transpose());
                    B = B.createBlockDiagonal(effProc(sn, station, jobClass, i, r).get(1).sumRows());
                }
            }
        }

        Matrix W = calculateW(psi, A, B, P);

        if (W.getNumRows() == 0 || W.getNumCols() == 0) {
            result.QN = new Matrix(M, K);
            result.UN = new Matrix(M, K);
            result.RN = new Matrix(M, K);
            result.TN = new Matrix(M, K);
            result.WN = new Matrix(0, 0);
            result.AN = new Matrix(0, 0);
            result.t = new Matrix(1, 1);
            result.QNt = new Matrix[M][K];
            result.UNt = new Matrix[M][K];
            result.TNt = new Matrix[M][K];
            for (int i = 0; i < M; i++) {
                for (int kk = 0; kk < K; kk++) {
                    result.QNt[i][kk] = new Matrix(1, 1);
                    result.UNt[i][kk] = new Matrix(1, 1);
                    result.TNt[i][kk] = new Matrix(1, 1);
                }
            }
            this.xvec_t = new Matrix(1, options.init_sol.length());
            this.xvec_it = new Matrix(1, options.init_sol.length());
            return;
        }

        Matrix sourceArrivals = new Matrix(M, K);
        for (int srcIst = 0; srcIst < M; srcIst++) {
            if (sn.sched.get(sn.stations.get(srcIst)) == SchedStrategy.EXT) {
                for (int r = 0; r < K; r++) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        sourceArrivals.set(srcIst, r, sn.rates.get(srcIst, r));
                    }
                }
            }
        }

        int totalStates = 0;
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                int np = (int) sn.phases.get(i, r);
                totalStates += (np == 0) ? 1 : np;
            }
        }

        Matrix ALambdaFull = new Matrix(totalStates, 1);
        int state = 0;
        for (int ist = 0; ist < M; ist++) {
            Station station = sn.stations.get(ist);
            for (int r = 0; r < K; r++) {
                JobClass jobClass = sn.jobclasses.get(r);
                int nPhases = (int) sn.phases.get(ist, r);
                if (nPhases > 0) {
                    if (sn.sched.get(station) == SchedStrategy.EXT) {
                        state += nPhases;
                    } else {
                        double arrivalRateToQueue = 0.0;
                        for (int srcIst = 0; srcIst < M; srcIst++) {
                            if (sourceArrivals.get(srcIst, r) > 0) {
                                arrivalRateToQueue += sourceArrivals.get(srcIst, r) * P.get(srcIst * K + r, ist * K + r);
                            }
                        }
                        if (arrivalRateToQueue > 0) {
                            Matrix pieVec = effPie(sn, station, jobClass, ist, r);
                            for (int kk = 0; kk < nPhases; kk++) {
                                ALambdaFull.set(state, 0, pieVec.get(kk) * arrivalRateToQueue);
                                state++;
                            }
                        } else {
                            state += nPhases;
                        }
                    }
                } else {
                    state++;
                }
            }
        }

        List<Integer> keep = new ArrayList<Integer>();
        Matrix sumWCols = W.sumCols();
        for (int col = 0; col < W.getNumCols(); col++) {
            if (!Double.isNaN(sumWCols.get(0, col))) {
                keep.add(col);
            }
        }

        Matrix WFiltered = new Matrix(keep.size(), keep.size());
        for (int i = 0; i < keep.size(); i++) {
            for (int j = 0; j < keep.size(); j++) {
                WFiltered.set(i, j, W.get(keep.get(i), keep.get(j)));
            }
        }

        Matrix ALambda = new Matrix(keep.size(), 1);
        for (int i = 0; i < keep.size(); i++) {
            ALambda.set(i, 0, ALambdaFull.get(keep.get(i), 0));
        }

        Matrix QaFull = new Matrix(1, totalStates);
        // (station, class) owning each ODE state, used to expand a per-(station,
        // class) rate schedule onto the state vector.
        int[] stationOfStateFull = new int[totalStates];
        int[] classOfStateFull = new int[totalStates];
        Matrix SQCFull = new Matrix(M * K, totalStates);
        Matrix SUCFull = new Matrix(M * K, totalStates);
        Matrix STCFull = new Matrix(M * K, totalStates);
        Matrix x0Build = new Matrix(totalStates, 1);

        state = 0;
        int initSolIdx = 0;
        for (int ist = 0; ist < M; ist++) {
            Station station = sn.stations.get(ist);
            for (int r = 0; r < K; r++) {
                JobClass jobClass = sn.jobclasses.get(r);
                int nPhases = (int) sn.phases.get(ist, r);
                if (nPhases == 0) {
                    QaFull.set(0, state, (double) ist);
                    stationOfStateFull[state] = ist;
                    classOfStateFull[state] = r;
                    state++;
                } else {
                    for (int kk = 0; kk < nPhases; kk++) {
                        QaFull.set(0, state, (double) ist);
                        stationOfStateFull[state] = ist;
                        classOfStateFull[state] = r;
                        SQCFull.set(ist * K + r, state, 1.0);
                        SUCFull.set(ist * K + r, state, 1.0 / S.get(ist, 0));
                        STCFull.set(ist * K + r, state, effProc(sn, station, jobClass, ist, r).get(1).sumRows(kk));
                        x0Build.set(state, 0, options.init_sol.get(0, initSolIdx));
                        initSolIdx++;
                        state++;
                    }
                }
            }
        }

        int nStatesFiltered = keep.size();
        Matrix Qa = new Matrix(1, nStatesFiltered);
        Matrix SQC = new Matrix(M * K, nStatesFiltered);
        Matrix SUC = new Matrix(M * K, nStatesFiltered);
        Matrix STC = new Matrix(M * K, nStatesFiltered);
        double[] x0 = new double[nStatesFiltered];

        for (int i = 0; i < nStatesFiltered; i++) {
            int ki = keep.get(i);
            Qa.set(0, i, QaFull.get(0, ki));
            for (int row = 0; row < M * K; row++) {
                SQC.set(row, i, SQCFull.get(row, ki));
                SUC.set(row, i, SUCFull.get(row, ki));
                STC.set(row, i, STCFull.get(row, ki));
            }
            x0[i] = x0Build.get(ki, 0);
        }

        // STOCHASTIC-COMPLEMENT THE INSTANTANEOUS STATES OUT OF THE LINEAR GENERATOR, the same
        // reduction the event-set routes take, so this route does not integrate an InfRate mode
        // either. Mirrors MATLAB solver_fluid_matrix.m / eliminate_immediate_matrix.m: W is reduced
        // by the complement, and every per-state quantity is projected onto the states that survive.
        // Recorded before the reduction: XVEC_IT is handed back to the caller as options.init_sol
        // for the next FCFS iterate, so it has to stay in the layout that iterate rebuilds, with
        // the eliminated coordinates present and empty. Shrinking it there is what made the next
        // pass index past the end of init_sol.
        int nStatesPreElim = nStatesFiltered;
        int[] keptAfterElim = null;
        if (FluidHideImmediate.resolve(sn, options)) {
            int nAll = nStatesFiltered;
            double immTol = GlobalConstants.Immediate * (1 - 1e-2);
            if (options.config != null && options.config.containsKey("immediate_tol")) {
                immTol = (Double) options.config.get("immediate_tol");
            }
            List<Integer> immList = new ArrayList<Integer>();
            List<Integer> timedList = new ArrayList<Integer>();
            for (int i = 0; i < nAll; i++) {
                double maxAbs = 0;
                for (int j = 0; j < nAll; j++) {
                    maxAbs = Math.max(maxAbs, Math.abs(WFiltered.get(i, j)));
                }
                if (maxAbs >= immTol) {
                    immList.add(i);
                } else {
                    timedList.add(i);
                }
            }
            if (!immList.isEmpty() && timedList.size() > 1) {
                int nT = timedList.size();
                int nI = immList.size();
                Matrix QTT = new Matrix(nT, nT);
                Matrix QTI = new Matrix(nT, nI);
                Matrix QIT = new Matrix(nI, nT);
                Matrix QII = new Matrix(nI, nI);
                for (int a = 0; a < nT; a++) {
                    for (int b = 0; b < nT; b++) {
                        QTT.set(a, b, WFiltered.get(timedList.get(a), timedList.get(b)));
                    }
                    for (int b = 0; b < nI; b++) {
                        QTI.set(a, b, WFiltered.get(timedList.get(a), immList.get(b)));
                    }
                }
                for (int a = 0; a < nI; a++) {
                    for (int b = 0; b < nT; b++) {
                        QIT.set(a, b, WFiltered.get(immList.get(a), timedList.get(b)));
                    }
                    for (int b = 0; b < nI; b++) {
                        QII.set(a, b, WFiltered.get(immList.get(a), immList.get(b)));
                    }
                }
                Matrix reduced = null;
                Matrix absorbII = null;
                Matrix sojournTI = null;
                try {
                    Matrix negQIIinv = QII.scale(-1.0).inv();
                    absorbII = negQIIinv.mult(QIT);         // absorption distribution of the block
                    // Expected mass held by each immediate state per unit mass in each timed one:
                    // x_I = x_T * Q_TI * (-Q_II)^-1. It is O(1/InfRate) and yet carries a FINITE
                    // throughput, because the read-off rate on those states is InfRate itself. That
                    // is why SQC/SUC/STC are corrected below rather than truncated: dropping the
                    // columns would silently delete every completion the instantaneous phase makes,
                    // and the station's throughput would stop balancing against its neighbours.
                    sojournTI = QTI.mult(negQIIinv);
                    reduced = QTT.add(1.0, QTI.mult(absorbII));
                } catch (Exception ex) {
                    reduced = null;
                }
                boolean finite = reduced != null;
                if (finite) {
                    for (int a = 0; a < nT && finite; a++) {
                        for (int b = 0; b < nT; b++) {
                            if (!Double.isFinite(reduced.get(a, b))) {
                                finite = false;
                                break;
                            }
                        }
                    }
                }
                if (!finite) {
                    line_warning("MatrixMethodAnalyzer", "stochastic complementation of the "
                            + "immediate states produced non-finite values; integrating the "
                            + "unreduced system instead");
                } else {
                    Matrix QaR = new Matrix(1, nT);
                    Matrix SQCR = new Matrix(M * K, nT);
                    Matrix SUCR = new Matrix(M * K, nT);
                    Matrix STCR = new Matrix(M * K, nT);
                    double[] x0R = new double[nT];
                    List<Integer> keepR = new ArrayList<Integer>();
                    for (int a = 0; a < nT; a++) {
                        int t = timedList.get(a);
                        QaR.set(0, a, Qa.get(0, t));
                        for (int row = 0; row < M * K; row++) {
                            double sq = SQC.get(row, t);
                            double su = SUC.get(row, t);
                            double st = STC.get(row, t);
                            for (int b = 0; b < nI; b++) {
                                double g = sojournTI.get(a, b);
                                if (g != 0) {
                                    int ii = immList.get(b);
                                    sq += g * SQC.get(row, ii);
                                    su += g * SUC.get(row, ii);
                                    st += g * STC.get(row, ii);
                                }
                            }
                            SQCR.set(row, a, sq);
                            SUCR.set(row, a, su);
                            STCR.set(row, a, st);
                        }
                        // The initial point is PROJECTED, not truncated: mass parked on an
                        // eliminated coordinate is jobs, and dropping it would lose population
                        // before the first step.
                        x0R[a] = x0[t];
                        // KEEP still indexes the FULL state space, so the rate-schedule expansion
                        // below keeps reading stationOfStateFull/classOfStateFull unchanged.
                        keepR.add(keep.get(t));
                    }
                    for (int b = 0; b < nI; b++) {
                        double mass = x0[immList.get(b)];
                        if (mass == 0) {
                            continue;
                        }
                        for (int a = 0; a < nT; a++) {
                            x0R[a] += mass * absorbII.get(b, a);
                        }
                    }
                    // The arrival vector follows the same complement: an arrival landing on an
                    // eliminated state arrives, instantaneously, where that state sends it.
                    Matrix ALambdaR = new Matrix(nT, 1);
                    for (int a = 0; a < nT; a++) {
                        ALambdaR.set(a, 0, ALambda.get(timedList.get(a), 0));
                    }
                    for (int b = 0; b < nI; b++) {
                        double li = ALambda.get(immList.get(b), 0);
                        if (li == 0) {
                            continue;
                        }
                        for (int a = 0; a < nT; a++) {
                            ALambdaR.set(a, 0, ALambdaR.get(a, 0) + li * absorbII.get(b, a));
                        }
                    }
                    WFiltered = reduced;
                    Qa = QaR;
                    SQC = SQCR;
                    SUC = SUCR;
                    STC = STCR;
                    x0 = x0R;
                    ALambda = ALambdaR;
                    keep = keepR;
                    nStatesFiltered = nT;
                    keptAfterElim = new int[nT];
                    for (int a = 0; a < nT; a++) {
                        keptAfterElim[a] = timedList.get(a);
                    }
                }
            }
        }

        // Expand any per-(station,class) rate schedule (options.config.rate_sched)
        // onto the kept ODE states. The fluid drift out of a state is linear in
        // that (station,class) service rate, so a time-varying rate is exactly a
        // per-state multiplicative factor m(t) on the effective service
        // population theta. Mirrors MATLAB solver_fluid_ratemult, which applies
        // the same factor per event of the closing ODE.
        FluidStateRateMultiplier stateMult = FluidStateRateMultiplier.build(
                options, keep, stationOfStateFull, classOfStateFull, sn);

        boolean[] isSourceState = new boolean[nStatesFiltered];
        for (int s = 0; s < nStatesFiltered; s++) {
            int ist = (int) Qa.get(0, s);
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                isSourceState[s] = true;
                x0[s] = 0.0;
            }
        }

        Matrix SQ = new Matrix(nStatesFiltered, nStatesFiltered);
        for (int s = 0; s < nStatesFiltered; s++) {
            int ist = (int) Qa.get(0, s);
            for (int col = 0; col < nStatesFiltered; col++) {
                if (Qa.get(0, col) == (double) ist) {
                    SQ.set(s, col, 1);
                }
            }
        }

        double minNonZeroRate = GlobalConstants.Inf;
        for (int i = 0; i < WFiltered.getNumRows(); i++) {
            for (int j = 0; j < WFiltered.getNumCols(); j++) {
                double tmpRate = FastMath.abs(WFiltered.get(i, j));
                if (tmpRate < minNonZeroRate && tmpRate > 0) {
                    minNonZeroRate = tmpRate;
                }
            }
        }
        double[] tRange = new double[]{
                options.timespan[0],
                FastMath.min(options.timespan[1], FastMath.abs(10 * options.iter_max / minNonZeroRate))
        };

        double[] initialState = x0.clone();
        double[] nextState = new double[nStatesFiltered];

        MatrixMethodODE rawOde;
        if (options.config.pstar.size() == 0) {
            rawOde = new MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState, stateMult);
        } else {
            rawOde = new MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState, sn, options.config.pstar, stateMult);
        }
        FirstOrderDifferentialEquations ode = rawOde;

        int Tmax = integrateDrift(ode, tRange, initialState, nextState, nStatesFiltered, options, result);

        // DEGENERATE DRIFT: re-integrate with a closed saturation term, do not
        // touch the answer that came back. min(E[n], c) is FLAT above the server
        // count, so a network of saturated stations has a CONTINUUM of fixed
        // points and this method returns whichever one the integrator stopped at
        // -- [9 1] against an exact [5 5] on two identical saturated stations in
        // a closed cycle, and [8 2] with two servers each. The repair is applied
        // to the DRIFT, not to the point: the same trajectory is integrated again
        // with E[min(n, c)] in place of min(E[n], c), which is strictly increasing
        // and so isolates one fixed point.
        //
        // WHY A CLOSURE AND NOT A SMOOTHED min: any smoothing sharp enough to stay
        // faithful to min away from the kink is numerically FLAT far from it. The
        // Boltzmann softmin at alpha = 20 carries a restoring force of exp(-160)
        // at the [9 1] point, and the p-norm trades the two off directly
        // (pstar = 2 recovers [5 5], pstar = 128 gives [8.94 1.06]). The closure
        // escapes the trade-off because its slope comes from the VARIANCE of the
        // marginal rather than from a smoothing width.
        //
        // Only a model that is ACTUALLY degenerate pays for it: the test is a
        // null-direction probe at the returned point, so a well-posed model
        // integrates once and is unchanged. See BUGS.md, _kb/06-solver-catalog.md.
        if (options.config.pstar.size() == 0 && Tmax > 0 && this.xvec_t != null) {
            boolean[] isInfState = new boolean[nStatesFiltered];
            for (int i = 0; i < nStatesFiltered; i++) {
                isInfState[i] = Double.isInfinite(sn.nservers.get((int) Qa.get(0, i), 0));
            }
            // The floor is 1e-12, NOT 1: flooring at 1 would make the null
            // threshold absolute for a slow model and read a live direction as
            // null. Same constant in all four codebases.
            double rateScale = 0.0;
            for (int r = 0; r < WFiltered.getNumRows(); r++) {
                for (int c = 0; c < WFiltered.getNumCols(); c++) {
                    rateScale = FastMath.max(rateScale, FastMath.abs(WFiltered.get(r, c)));
                }
            }
            Matrix xFinal = Matrix.extractRows(this.xvec_t, Tmax - 1, Tmax, null).transpose();
            if (fixedPointIsDegenerate(rawOde, xFinal, SQC, K, isSourceState, isInfState, rateScale)) {
                Matrix savedT = result.t;
                Matrix savedX = this.xvec_t;
                MatrixMethodODE closedOde = new MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda,
                        nStatesFiltered, isSourceState, stateMult);
                closedOde.setVarClosure(isInfState);
                boolean repaired = false;
                try {
                    int TmaxC = integrateDrift(closedOde, tRange, x0.clone(),
                            new double[nStatesFiltered], nStatesFiltered, options, result);
                    if (TmaxC > 0 && this.xvec_t != null) {
                        boolean finite = true;
                        for (int j = 0; j < this.xvec_t.getNumCols(); j++) {
                            double v = this.xvec_t.get(TmaxC - 1, j);
                            if (Double.isNaN(v) || Double.isInfinite(v)) {
                                finite = false;
                                break;
                            }
                        }
                        if (finite) {
                            Tmax = TmaxC;
                            rawOde = closedOde;
                            repaired = true;
                        }
                    }
                } catch (RuntimeException e) {
                    // fall through and restore
                }
                if (!repaired) {
                    // A failed repair leaves the unrepaired answer standing rather
                    // than turning a wrong number into no number.
                    result.t = savedT;
                    this.xvec_t = savedX;
                }
            }
        }

        this.xvec_it = Matrix.extractRows(xvec_t, Tmax - 1, Tmax, null);
        if (keptAfterElim != null && this.xvec_it.getNumCols() == keptAfterElim.length) {
            Matrix expanded = new Matrix(1, nStatesPreElim);
            for (int a = 0; a < keptAfterElim.length; a++) {
                expanded.set(0, keptAfterElim[a], this.xvec_it.get(0, a));
            }
            this.xvec_it = expanded;
        }

        result.QNt = new Matrix[M][K];
        result.UNt = new Matrix[M][K];
        result.TNt = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            for (int kk = 0; kk < K; kk++) {
                result.QNt[i][kk] = new Matrix(Tmax, 1);
                result.UNt[i][kk] = new Matrix(Tmax, 1);
                result.TNt[i][kk] = new Matrix(Tmax, 1);
            }
        }

        for (int step = 0; step < Tmax; step++) {
            Matrix x = Matrix.extractRows(xvec_t, step, step + 1, null);
            x = x.transpose();
            // The same share the drift used, smoothed or not
            Matrix theta = rawOde.theta(x);

            Matrix QNtmp = SQC.mult(x, null);
            Matrix UNtmp = SUC.mult(theta, null);
            // Under a time-varying rate schedule the completion rate baked into
            // STC is the nominal one, so the throughput is scaled by the state
            // multiplier at this time point. The queue length is rate-free and
            // the utilization is an occupancy ratio, so neither is scaled.
            Matrix thetaT = theta;
            if (stateMult != null) {
                double[] m = stateMult.multAt(result.t.get(step, 0));
                thetaT = theta.copy();
                for (int phase = 0; phase < nStatesFiltered; phase++) {
                    thetaT.set(phase, 0, theta.get(phase, 0) * m[phase]);
                }
            }
            Matrix TNtmp = STC.mult(thetaT, null);

            for (int i = 0; i < M; i++) {
                for (int kk = 0; kk < K; kk++) {
                    result.QNt[i][kk].set(step, 0, QNtmp.get(i * K + kk, 0));
                    result.UNt[i][kk].set(step, 0, UNtmp.get(i * K + kk, 0));
                    result.TNt[i][kk].set(step, 0, TNtmp.get(i * K + kk, 0));
                }
            }
        }

        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.WN = new Matrix(0, 0);
        result.AN = new Matrix(0, 0);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                result.QN.set(i, j, result.QNt[i][j].get(Tmax - 1, 0));
                result.UN.set(i, j, result.UNt[i][j].get(Tmax - 1, 0));
                result.TN.set(i, j, result.TNt[i][j].get(Tmax - 1, 0));
                result.RN.set(i, j, result.QN.get(i, j) / result.TN.get(i, j));
            }
        }

        for (int ist = 0; ist < M; ist++) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                for (int r = 0; r < K; r++) {
                    if (!Double.isNaN(sn.rates.get(ist, r)) && sn.rates.get(ist, r) > 0) {
                        result.TN.set(ist, r, sn.rates.get(ist, r));
                    }
                }
            }
        }
    }

    @Override
    public Matrix getXVecIt() {
        return this.xvec_it;
    }

    private Matrix calculateW(Matrix psi, Matrix A, Matrix B, Matrix P) {
        MatrixEquation calc = new MatrixEquation();
        calc.alias(psi, "psi", A, "A", P, "P", B, "B");
        calc.process("W = psi + B*P*A'");
        // see _kb/06-solver-catalog.md (JAR-only implementation notes: MatrixMethodAnalyzer full-size W / NaN-column strip fix)
        return new Matrix(calc.lookupSimple("W"));
    }

    /**
     * Runs one integration of the drift, filling result.t and this.xvec_t, and
     * returns the number of time points. Factored out so the degeneracy repair
     * can integrate the SAME trajectory a second time with a closed saturation
     * term rather than adjusting the point that came back.
     */
    private int integrateDrift(FirstOrderDifferentialEquations ode, double[] tRange,
                               double[] initialState, double[] nextState,
                               int nStatesFiltered, SolverOptions options, SolverResult result) {
        int Tmax;
        double[] y0 = initialState.clone();
        if (options.stiff) {
            LSODA odeSolver = options.odesolvers.stiffIntegratorFor(tRange[0], tRange[1], options.tol,
                    options.tol > GlobalConstants.CoarseTol);
            try {
                odeSolver.integrate(ode, tRange[0], y0, tRange[1], nextState);
            } catch (RuntimeException e) {
                if (odeSolver.getStepsTaken() > 0 && options.verbose != VerboseLevel.SILENT) {
                    System.out.println("WARNING: LSODA step limit reached, using partial results ("
                            + odeSolver.getStepsTaken() + " steps)");
                }
            }
            Tmax = odeSolver.getStepsTaken() + 1;
            if (odeSolver.getStepsTaken() == 0) {
                result.t = new Matrix(1, 1);
                result.t.set(0, 0, tRange[1] > 0 ? tRange[1] : 1.0);
                this.xvec_t = new Matrix(1, nStatesFiltered);
                for (int i = 0; i < nStatesFiltered; i++) {
                    this.xvec_t.set(0, i, Math.max(0.0, y0[i]));
                }
            } else {
                List<Double> tHistory = odeSolver.getTvec();
                java.util.ArrayList<Double[]> yHistory = odeSolver.getYvec();
                int dim = ode.getDimension();
                DMatrixRMaj denseT = new DMatrixRMaj(Tmax, 1);
                DMatrixRMaj denseX = new DMatrixRMaj(Tmax, dim);
                for (int i = 0; i < Tmax; i++) {
                    denseT.set(i, 0, tHistory.get(i));
                    Double[] row = yHistory.get(i);
                    for (int j = 0; j < dim; j++) {
                        denseX.set(i, j, Math.max(0.0, row[j]));
                    }
                }
                result.t = new Matrix(DConvertMatrixStruct.convert(denseT, (DMatrixSparseCSC) null, 0.0));
                this.xvec_t = new Matrix(DConvertMatrixStruct.convert(denseX, (DMatrixSparseCSC) null, 0.0));
            }
        } else {
            FirstOrderIntegrator odeSolver = options.odesolvers.integratorFor(tRange[0], tRange[1],
                    options.tol, options.tol > GlobalConstants.CoarseTol);
            odeSolver.clearStepHandlers();
            TransientDataHandler stepHandler = new TransientDataHandler(nStatesFiltered);
            odeSolver.addStepHandler(stepHandler);
            odeSolver.integrate(ode, tRange[0], y0, tRange[1], nextState);
            result.t = stepHandler.tVec;
            Tmax = result.t.getNumRows();
            this.xvec_t = stepHandler.xVec;
        }
        return Tmax;
    }

    /**
     * Is the returned point one of a CONTINUUM of fixed points?
     *
     * <p>A station whose queue exceeds its server count has theta pinned at the
     * server count, so the drift cannot tell one split of the mass between two
     * such stations from another. The test is direct rather than structural: move
     * a little mass from one station to another along a POPULATION-CONSERVING
     * direction and see whether the drift moves at all. Both directions are
     * tried, because the integrator typically stops on the BOUNDARY of the
     * degenerate set, where one of the two does change the drift.</p>
     *
     * <p>The DIRECTIONAL DERIVATIVE is the scale-free quantity to threshold: a
     * live direction moves the drift at the station's own service rate (measured
     * 9.99e-01) and a null one only by the FineTol the share carries (measured
     * 1.00e-08), four orders apart. Source and INF states are excluded: neither
     * can be the pinned coordinate. Returns false whenever the point is not a
     * fixed point at all, so a transient run is never repaired.</p>
     */
    /**
     * Is the returned point one of a CONTINUUM of fixed points?
     *
     * <p>A station whose queue exceeds its server count has theta pinned at the
     * server count, so the drift cannot tell one split of the mass between two
     * such stations from another. The test is direct: move a little mass of ONE
     * CLASS from one station to another along a population-conserving direction
     * and see whether the drift moves at all. Both directions are tried, because
     * the integrator typically stops on the BOUNDARY of the degenerate set.</p>
     *
     * <p>PER CLASS, NOT PER STATION. A direction that moves a station's mass
     * across ALL its classes is not one the model can take: a SelfLoopingClass is
     * pinned at one station and can never leave, so the direction is infeasible,
     * the drift is trivially unchanged along it, and a well-posed model reads as
     * degenerate -- which is what it did to
     * sanity_CQN_2q_psfcfs_1class_1slcateachqueue, RespT 1.4336 against a
     * baseline of 0.726303. Moving ONE class between two stations it occupies IS
     * feasible, and a self-looping class occupies exactly one station.</p>
     *
     * <p>The DIRECTIONAL DERIVATIVE is the scale-free quantity to threshold: a
     * live direction moves the drift at the station's own service rate and a null
     * one only by the FineTol the share carries, four orders apart. Source and
     * INF states are excluded: neither can be the pinned coordinate. Returns
     * false whenever the point is not a fixed point at all.</p>
     */
    private static boolean fixedPointIsDegenerate(FirstOrderDifferentialEquations ode, Matrix x,
                                                  Matrix SQC, int K, boolean[] isSourceState,
                                                  boolean[] isInfState, double rateScale) {
        int n = x.getNumRows();
        double[] xv = new double[n];
        double maxAbsX = 1.0;
        for (int i = 0; i < n; i++) {
            xv[i] = x.get(i, 0);
            maxAbsX = FastMath.max(maxAbsX, FastMath.abs(xv[i]));
        }
        double[] d0 = new double[n];
        ode.computeDerivatives(0.0, xv, d0);
        double maxD0 = 0.0;
        for (int i = 0; i < n; i++) {
            maxD0 = FastMath.max(maxD0, FastMath.abs(d0[i]));
        }
        if (maxD0 > 1e-6 * maxAbsX || K <= 0) {
            return false;
        }
        int M = SQC.getNumRows() / K;
        double step = 1e-3 * maxAbsX;
        double scale = FastMath.max(rateScale, 1e-12);
        double[] xp = new double[n];
        double[] dp = new double[n];
        for (int r = 0; r < K; r++) {
            List<int[]> groups = new ArrayList<int[]>();
            for (int i = 0; i < M; i++) {
                List<Integer> members = new ArrayList<Integer>();
                for (int a = 0; a < n; a++) {
                    boolean src = isSourceState != null && a < isSourceState.length && isSourceState[a];
                    boolean inf = isInfState != null && a < isInfState.length && isInfState[a];
                    if (!src && !inf && SQC.get(i * K + r, a) > 0) {
                        members.add(a);
                    }
                }
                if (!members.isEmpty()) {
                    int[] g = new int[members.size()];
                    for (int t = 0; t < g.length; t++) {
                        g[t] = members.get(t);
                    }
                    groups.add(g);
                }
            }
            if (groups.size() < 2) {
                continue;
            }
            double[] mass = new double[groups.size()];
            for (int a = 0; a < groups.size(); a++) {
                for (int idx : groups.get(a)) {
                    mass[a] += xv[idx];
                }
            }
            for (int a = 0; a < groups.size(); a++) {
                if (mass[a] <= step) {
                    continue;
                }
                for (int b = 0; b < groups.size(); b++) {
                    if (a == b) {
                        continue;
                    }
                    System.arraycopy(xv, 0, xp, 0, n);
                    for (int idx : groups.get(a)) {
                        xp[idx] -= step * xv[idx] / mass[a];      // take, proportionally
                    }
                    int[] gb = groups.get(b);
                    for (int idx : gb) {
                        if (mass[b] > 0) {
                            xp[idx] += step * xv[idx] / mass[b];  // give, proportionally
                        } else {
                            xp[idx] += step / gb.length;
                        }
                    }
                    ode.computeDerivatives(0.0, xp, dp);
                    double maxDd = 0.0;
                    for (int i = 0; i < n; i++) {
                        maxDd = FastMath.max(maxDd, FastMath.abs(dp[i] - d0[i]));
                    }
                    if (maxDd / step <= 1e-4 * scale) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
