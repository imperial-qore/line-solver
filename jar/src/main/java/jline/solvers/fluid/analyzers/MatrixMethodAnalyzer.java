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

        Matrix Sa = new Matrix(nStatesFiltered, 1);
        for (int i = 0; i < nStatesFiltered; i++) {
            Sa.set(i, 0, S.get((int) Qa.get(0, i), 0));
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

        FirstOrderDifferentialEquations rawOde;
        if (options.config.pstar.size() == 0) {
            rawOde = new MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState, stateMult);
        } else {
            rawOde = new MatrixMethodODE(WFiltered, SQ, S, Qa, ALambda, nStatesFiltered, isSourceState, sn, options.config.pstar, stateMult);
        }
        FirstOrderDifferentialEquations ode = rawOde;

        int Tmax;
        if (options.stiff) {
            LSODA odeSolver;
            if (options.tol > GlobalConstants.CoarseTol) {
                odeSolver = options.odesolvers.fastStiffODESolver;
            } else {
                odeSolver = options.odesolvers.accurateStiffODESolver;
            }
            try {
                odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
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
                    this.xvec_t.set(0, i, Math.max(0.0, initialState[i]));
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
            FirstOrderIntegrator odeSolver = (options.tol > GlobalConstants.CoarseTol)
                    ? options.odesolvers.fastODESolver
                    : options.odesolvers.accurateODESolver;
            odeSolver.clearStepHandlers();
            TransientDataHandler stepHandler = new TransientDataHandler(nStatesFiltered);
            odeSolver.addStepHandler(stepHandler);
            odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
            result.t = stepHandler.tVec;
            Tmax = result.t.getNumRows();
            this.xvec_t = stepHandler.xVec;
        }
        this.xvec_it = Matrix.extractRows(xvec_t, Tmax - 1, Tmax, null);

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
            Matrix theta = x.copy();
            Matrix SQx = SQ.mult(x, null);
            for (int phase = 0; phase < nStatesFiltered; phase++) {
                if (isSourceState[phase]) {
                    theta.set(phase, 0, 0.0);
                } else {
                    double valSQx = SQx.get(phase, 0) + GlobalConstants.FineTol;
                    double valSa = Sa.get(phase, 0);
                    theta.set(phase, 0, x.get(phase, 0) / valSQx * FastMath.min(valSa, valSQx));
                }
            }

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
}
