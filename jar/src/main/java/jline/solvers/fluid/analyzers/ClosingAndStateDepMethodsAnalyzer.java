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
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.FluidNhpp;
import jline.solvers.fluid.handlers.PassageTimeODE;
import jline.solvers.fluid.handlers.TransientDataHandler;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import odesolver.LSODA;

public class ClosingAndStateDepMethodsAnalyzer implements FluidAnalyzer {

    public Matrix xvec_t;
    public Matrix xvec_it;

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
        FirstOrderDifferentialEquations ode = new PassageTimeODE(sn, mu, phi, proc, sn.rt, S, options);

        options.stiff = detectStiffnessUsingOstrowski(sn, slowrate);

        double T0 = options.timespan[0];
        double T = 0.0;

        List<Matrix> tIterations = new LinkedList<Matrix>();
        List<Matrix> xVecIterations = new LinkedList<Matrix>();

        boolean goon = true;

        while ((Double.isFinite(options.timespan[1]) && T < options.timespan[1]) || (goon && iter < options.iter_max)) {
            iter++;

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
                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("The initial point is invalid, Fluid solver switching to default initialization.");
                    }
                    odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState);
                }
                if (odeSolver.getStepsTaken() > 0) {
                    java.util.ArrayList<Double> tHistory = odeSolver.getTvec();
                    java.util.ArrayList<Double[]> yHistory = odeSolver.getYvec();
                    int dim = ode.getDimension();
                    Tmax = odeSolver.getStepsTaken() + 1;
                    int lastIdx = Tmax - 1;
                    this.xvec_it = new Matrix(1, dim);
                    for (int j = 0; j < dim; j++) this.xvec_it.set(0, j, Math.max(0.0, yHistory.get(lastIdx)[j]));
                    DMatrixRMaj denseT = new DMatrixRMaj(Tmax, 1);
                    DMatrixRMaj denseX = new DMatrixRMaj(Tmax, dim);
                    for (int i = 0; i < Tmax; i++) {
                        denseT.set(i, 0, tHistory.get(i));
                        for (int j = 0; j < dim; j++) denseX.set(i, j, Math.max(0.0, yHistory.get(i)[j]));
                    }
                    tIterations.add(new Matrix(denseT));
                    xVecIterations.add(new Matrix(denseX));
                }
            } else {
                FirstOrderIntegrator odeSolver;
                if (options.tol > GlobalConstants.CoarseTol) {
                    odeSolver = options.odesolvers.fastODESolver;
                } else {
                    odeSolver = options.odesolvers.accurateODESolver;
                }
                odeSolver.clearStepHandlers();
                TransientDataHandler stepHandler = new TransientDataHandler(initialState.length);
                odeSolver.addStepHandler(stepHandler);

                boolean usedStiffFallback = false;
                try {
                    try {
                        odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                    } catch (RuntimeException e) {
                        if (e.getMessage() != null && e.getMessage().contains("step size")) {
                            usedStiffFallback = true;
                        } else {
                            if (options.verbose != VerboseLevel.SILENT) {
                                System.out.println("The initial point is invalid, Fluid solver switching to default initialization.");
                            }
                            odeSolver.clearStepHandlers();
                            odeSolver.addStepHandler(stepHandler);
                            odeSolver.integrate(ode, tRange[0], yDefault, tRange[1], nextState);
                        }
                    }
                } catch (RuntimeException e) {
                    usedStiffFallback = true;
                }

                if (usedStiffFallback) {
                    LSODA stiffSolver;
                    if (options.tol > GlobalConstants.CoarseTol) {
                        stiffSolver = options.odesolvers.fastStiffODESolver;
                    } else {
                        stiffSolver = options.odesolvers.accurateStiffODESolver;
                    }
                    stiffSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                    java.util.ArrayList<Double> tHistory = stiffSolver.getTvec();
                    java.util.ArrayList<Double[]> yHistory = stiffSolver.getYvec();
                    int dim = ode.getDimension();
                    Tmax = stiffSolver.getStepsTaken() + 1;
                    int lastIdx = Tmax - 1;
                    this.xvec_it = new Matrix(1, dim);
                    for (int j = 0; j < dim; j++) this.xvec_it.set(0, j, Math.max(0.0, yHistory.get(lastIdx)[j]));

                    DMatrixRMaj denseT = new DMatrixRMaj(Tmax, 1);
                    DMatrixRMaj denseX = new DMatrixRMaj(Tmax, dim);
                    for (int i = 0; i < Tmax; i++) {
                        denseT.set(i, 0, tHistory.get(i));
                        for (int j = 0; j < dim; j++) denseX.set(i, j, Math.max(0.0, yHistory.get(i)[j]));
                    }
                    tIterations.add(new Matrix(denseT));
                    xVecIterations.add(new Matrix(denseX));
                } else {
                    tIterations.add(stepHandler.tVec);
                    xVecIterations.add(stepHandler.xVec);
                    Tmax = stepHandler.tVec.getNumRows();
                    this.xvec_it = Matrix.extractRows(stepHandler.xVec, Tmax - 1, Tmax, null);
                }
            }
            totalSteps += Tmax;

            T0 = T;

            if (T >= options.timespan[1]) {
                goon = false;
            }
        }

        if (!xVecIterations.isEmpty() && totalSteps > 0) {
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

    private void solver_fluid(NetworkStruct sn, SolverOptions options, SolverResult result) {
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
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int k = 0; k < K; k++) {
                JobClass jobClass = sn.jobclasses.get(k);

                if (mu.get(station).get(jobClass).hasNaN()) {
                    mu.get(station).put(jobClass, new Matrix(0, 0));
                    phi.get(station).put(jobClass, new Matrix(0, 0));
                }

                if (sn.rt.sumCols((i * K) + k) > 0) {
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

                                if ("default".equals(options.method) || "closing".equals(options.method)) {
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

                                if ("default".equals(options.method) || "closing".equals(options.method)) {
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

        // Response Times
        result.RN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (result.TN.get(i, k) > 0) {
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

    public boolean detectStiffnessUsingOstrowski(NetworkStruct sn, Matrix rate) {
        Matrix transitionMatrix = sn.rt.copy();
        for (int i = 0; i < transitionMatrix.getNumRows(); i++) {
            double p = 0.0;
            for (int j = 0; j < transitionMatrix.getNumCols(); j++) {
                if (i != j) p -= transitionMatrix.get(i, j);
            }
            transitionMatrix.set(i, i, p);
        }

        Matrix expandRate = new Matrix(1, rate.getNumCols() * rate.getNumRows());
        for (int i = 0; i < rate.getNumCols() * rate.getNumRows(); i++) {
            int r = i / rate.getNumCols();
            int c = i % rate.getNumCols();
            expandRate.set(0, i, rate.get(r, c));
        }
        transitionMatrix = transitionMatrix.elementMultWithVector(expandRate);
        double[][] bound = new double[transitionMatrix.getNumCols()][2];
        int n = transitionMatrix.getNumCols();
        double alpha = 0.5;
        boolean stiff = false;

        for (int i = 0; i < n; i++) {
            bound[i][0] = transitionMatrix.get(i, i);
            double rSum = transitionMatrix.sumAbsRows(i) - FastMath.abs(transitionMatrix.get(i, i));
            double cSum = transitionMatrix.sumAbsCols(i) - FastMath.abs(transitionMatrix.get(i, i));
            bound[i][1] = FastMath.pow(rSum, alpha) * FastMath.pow(cSum, 1 - alpha);
            stiff = stiff || (bound[i][0] < 0 && FastMath.abs(bound[i][1]) < FastMath.abs(bound[i][0]));
        }
        return stiff;
    }
}
