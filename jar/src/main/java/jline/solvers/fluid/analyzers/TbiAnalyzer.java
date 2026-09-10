/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers;

import java.util.ArrayList;
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
import jline.api.sn.SnRtStations;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.FluidInterp;
import jline.solvers.fluid.FluidNhpp;
import jline.solvers.fluid.handlers.PassageTimeODE;
import jline.solvers.fluid.handlers.TransientDataHandler;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import odesolver.LSODA;

/**
 * Trajectory-based iteration (TBI) analyzer for the transient fluid solution.
 *
 * <p>The station set is partitioned into cells (see {@link #tbiPartition}). On
 * each growing-horizon time segment, the initial value problem of every cell is
 * solved with the state of the other cells frozen at the trajectory computed in
 * the previous sweep (Jacobi waveform relaxation); sweeps repeat until the
 * trajectory sup-norm gap falls below {@code options.config.tbi_tol}. Cross-cell
 * inflows are therefore evaluated on frozen trajectories, interior flows on the
 * live cell state, matching the decomposed closing-method ODEs of TBI.</p>
 *
 * <p>Only the ODE integration loop is overridden; the state initialization and
 * the QN/QNt/UNt/TN/UN postprocessing are inherited from
 * {@link ClosingAndStateDepMethodsAnalyzer}, run under the closing method
 * convention.</p>
 *
 * <p>Reference: Sheldon, Tuncer, Casale, "TBI: Transient Hierarchical Modeling
 * of Large-Scale Vehicle Sharing Systems", IEEE T-ITS.</p>
 */
public class TbiAnalyzer extends ClosingAndStateDepMethodsAnalyzer {

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        // TBI reuses the closing-method postprocessing (QN/TN/UN); the FCFS and
        // priority throughput branches key on the "closing"/"default" method
        // string, so present the closing convention downstream.
        options.method = "closing";
        super.analyze(sn, options, result);
        result.method = "tbi";
    }

    @Override
    protected void solver_fluid_iteration(NetworkStruct sn,
                                          Map<Station, Map<JobClass, Matrix>> mu,
                                          Map<Station, Map<JobClass, Matrix>> phi,
                                          Matrix S,
                                          double[] yDefault,
                                          Matrix slowrate,
                                          SolverOptions options,
                                          SolverResult result) {

        double tbiTol = 1e-3;
        int tbiIterMax = 50;
        if (options.config != null) {
            tbiTol = options.config.tbi_tol;
            tbiIterMax = options.config.tbi_iter_max;
        }

        // NHPP sources carry a rate schedule rather than a {D0, D1} MAP; see
        // ClosingAndStateDepMethodsAnalyzer.solver_fluid_iteration.
        Map<Station, Map<JobClass, MatrixCell>> proc = FluidNhpp.substituteNhppProc(sn, mu, phi);
        // station-major routing; sn.rt is indexed by stateful node (see SnRtStations)
        PassageTimeODE pt = new PassageTimeODE(sn, mu, phi, proc, SnRtStations.snRtStations(sn).getLeft(), S, options);
        int ndim = pt.getDimension();
        Matrix allJumps = pt.getAllJumps();
        Matrix qIndices = pt.getQIndices();
        Matrix Kic = pt.getKic();

        // Partition stations into cells and build per-cell state masks. The
        // restriction of the jump matrix to a cell's rows retains interior and
        // boundary events only through their effect on the cell state, so
        // inbound events read the frozen complement while outbound mass leaves.
        List<int[]> cells = tbiPartition(sn, options);
        int ncells = cells.size();
        int[][] cellmask = new int[ncells][];
        for (int kc = 0; kc < ncells; kc++) {
            List<Integer> mask = new ArrayList<Integer>();
            int[] cell = cells.get(kc);
            for (int ci = 0; ci < cell.length; ci++) {
                int i = cell[ci];
                for (int c = 0; c < sn.nclasses; c++) {
                    int kicIC = (int) Kic.get(i, c);
                    if (kicIC > 0) {
                        int start = (int) qIndices.get(i, c);
                        for (int p = 0; p < kicIC; p++) {
                            mask.add(start + p);
                        }
                    }
                }
            }
            int[] idx = new int[mask.size()];
            for (int p = 0; p < idx.length; p++) {
                idx[p] = mask.get(p);
            }
            cellmask[kc] = idx;
        }

        // Horizon growth heuristic, mirroring solver_fluid_iteration.
        double minNonZeroRate = GlobalConstants.Inf;
        for (int i = 0; i < slowrate.getNumRows(); i++) {
            for (int j = 0; j < slowrate.getNumCols(); j++) {
                double r = slowrate.get(i, j);
                if (r > GlobalConstants.CoarseTol && Double.isFinite(r) && r < minNonZeroRate) {
                    minNonZeroRate = r;
                }
            }
        }
        if (!Double.isFinite(minNonZeroRate)) {
            minNonZeroRate = 1.0; // fallback when all rates are zero or infinite
        }

        List<Matrix> tIterations = new LinkedList<Matrix>();
        List<Matrix> xVecIterations = new LinkedList<Matrix>();
        int totalSteps = 0;

        double T0 = options.timespan[0];
        double T = 0.0;
        int iter = 0;
        boolean goon = true;

        // Only a finite timespan resolves the transient between the endpoints; the
        // per-segment waveform still has to be assembled, but retaining every segment
        // costs O(iter_max * steps) to describe a fixed point (see the closing analyzer)
        boolean keepTrajectory = Double.isFinite(options.timespan[1]);
        Matrix xvecInit = this.xvec_it.copy();

        // Wall-clock budget (options.timeout, seconds; Inf = none), as MATLAB
        // solver_fluid_tbi_iteration.m guards both the segment loop and the sweep loop
        double maxTime = (Double.isFinite(options.timeout) && options.timeout > 0) ? options.timeout
                : GlobalConstants.Inf;
        long startNanos = System.nanoTime();

        while ((Double.isFinite(options.timespan[1]) && T < options.timespan[1])
                || (goon && iter < options.iter_max)) {
            iter++;

            if ((System.nanoTime() - startNanos) / 1e9 > maxTime) {
                goon = false;
                break;
            }

            double[] y0 = new double[ndim];
            for (int i = 0; i < ndim; i++) {
                y0[i] = xvec_it.get(0, i);
            }

            if (iter == 1) {
                T = FastMath.min(options.timespan[1], FastMath.abs(10.0 / minNonZeroRate));
            } else {
                T = FastMath.min(options.timespan[1], FastMath.abs(10.0 * iter / minNonZeroRate));
            }

            // Frozen trajectory on this segment, initialized constant at the
            // segment entry state (warm start of the waveform relaxation).
            double[] frozenT = new double[] { T0, T };
            double[][] frozenY = new double[][] { y0.clone(), y0.clone() };

            double delta = GlobalConstants.Inf;
            for (int sweep = 0; sweep < tbiIterMax; sweep++) {
                double[][] cellT = new double[ncells][];
                double[][][] cellY = new double[ncells][][];
                double[] tgrid = frozenT.clone();

                for (int kc = 0; kc < ncells; kc++) {
                    int[] idx = cellmask[kc];
                    if (idx.length == 0) {
                        continue;
                    }
                    double[] y0c = gather(y0, idx);
                    CellODE odeC = new CellODE(pt, allJumps, idx, ndim, frozenT, frozenY);
                    CellTraj traj = integrateCell(odeC, T0, y0c, T, options, yDefault, idx);
                    cellT[kc] = traj.t;
                    cellY[kc] = traj.y;
                    tgrid = unionSorted(tgrid, traj.t);
                }

                // Assemble the new full trajectory on the union time grid.
                double[][] Ynew = new double[tgrid.length][ndim];
                for (int kc = 0; kc < ncells; kc++) {
                    int[] idx = cellmask[kc];
                    if (idx.length == 0) {
                        continue;
                    }
                    double[] tc = cellT[kc];
                    double[][] yc = cellY[kc];
                    for (int r = 0; r < tgrid.length; r++) {
                        double tq = FastMath.min(FastMath.max(tgrid[r], tc[0]), tc[tc.length - 1]);
                        double[] yi = FluidInterp.interp(tc, yc, tq, idx.length);
                        for (int p = 0; p < idx.length; p++) {
                            Ynew[r][idx[p]] = yi[p];
                        }
                    }
                }

                // Sup-norm gap against the previous sweep trajectory.
                delta = 0.0;
                for (int r = 0; r < tgrid.length; r++) {
                    double[] yold = FluidInterp.interp(frozenT, frozenY, tgrid[r], ndim);
                    for (int j = 0; j < ndim; j++) {
                        double d = FastMath.abs(Ynew[r][j] - yold[j]);
                        if (d > delta) {
                            delta = d;
                        }
                    }
                }

                frozenT = tgrid;
                frozenY = Ynew;
                if (delta < tbiTol || (System.nanoTime() - startNanos) / 1e9 > maxTime) {
                    break;
                }
            }

            if (delta >= tbiTol && options.verbose != VerboseLevel.SILENT) {
                System.out.println("TBI sweeps did not converge within tbi_iter_max=" + tbiIterMax
                        + " on segment [" + T0 + "," + T + "], residual gap " + delta + ".");
            }

            int steps = frozenT.length;
            if (keepTrajectory) {
                DMatrixRMaj denseT = new DMatrixRMaj(steps, 1);
                DMatrixRMaj denseX = new DMatrixRMaj(steps, ndim);
                for (int r = 0; r < steps; r++) {
                    denseT.set(r, 0, frozenT[r]);
                    for (int j = 0; j < ndim; j++) {
                        denseX.set(r, j, FastMath.max(0.0, frozenY[r][j]));
                    }
                }
                tIterations.add(new Matrix(denseT));
                xVecIterations.add(new Matrix(denseX));
                totalSteps += steps;
            }

            this.xvec_it = new Matrix(1, ndim);
            for (int j = 0; j < ndim; j++) {
                this.xvec_it.set(0, j, FastMath.max(0.0, frozenY[steps - 1][j]));
            }

            T0 = T;
            if (T >= options.timespan[1]) {
                goon = false;
            }
        }

        if (!keepTrajectory) {
            // QNt/UNt/TNt are read at row 0 (initial condition) and at the last row (fixed
            // point), so the two endpoints carry the whole contract of an infinite timespan
            this.xvec_t = new Matrix(2, ndim);
            for (int j = 0; j < ndim; j++) {
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
     * Partition the station set into cells for trajectory-based iteration.
     * Honors {@code options.config.tbi_cells}, a list of disjoint zero-based
     * station-index vectors covering 0..nstations-1. Otherwise stations are
     * agglomerated greedily on the symmetrized station-level routing weights,
     * targeting {@code options.config.tbi_cellsize} stations per cell.
     *
     * @param sn      network struct
     * @param options solver options
     * @return list of cells, each an array of zero-based station indices
     */
    public static List<int[]> tbiPartition(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        if (options.config != null && options.config.tbi_cells != null && !options.config.tbi_cells.isEmpty()) {
            List<int[]> explicit = options.config.tbi_cells;
            boolean[] covered = new boolean[M];
            int count = 0;
            for (int c = 0; c < explicit.size(); c++) {
                int[] cell = explicit.get(c);
                for (int p = 0; p < cell.length; p++) {
                    int i = cell[p];
                    if (i < 0 || i >= M || covered[i]) {
                        throw new RuntimeException("options.config.tbi_cells must be a partition of the station set 0.."
                                + (M - 1) + ".");
                    }
                    covered[i] = true;
                    count++;
                }
            }
            if (count != M) {
                throw new RuntimeException("options.config.tbi_cells must be a partition of the station set 0.."
                        + (M - 1) + ".");
            }
            return explicit;
        }

        int cellsize = 5;
        if (options.config != null && options.config.tbi_cellsize > 0) {
            cellsize = options.config.tbi_cellsize;
        }

        // Station-level coupling weights, aggregated over classes and symmetrized.
        // station-major routing, hoisted: the reduction inverts a matrix
        Matrix rtSt = SnRtStations.snRtStations(sn).getLeft();
        double[][] A = new double[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                double w = 0.0;
                for (int c = 0; c < K; c++) {
                    for (int l = 0; l < K; l++) {
                        w += rtSt.get(i * K + c, j * K + l);
                    }
                }
                A[i][j] = w;
            }
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (i != j) {
                    A[i][j] = A[i][j] + A[j][i];
                }
            }
        }
        for (int i = 0; i < M; i++) {
            A[i][i] = 0.0;
        }

        int ncellsTarget = FastMath.max(1, (int) FastMath.ceil((double) M / cellsize));
        List<List<Integer>> cells = new ArrayList<List<Integer>>();
        for (int i = 0; i < M; i++) {
            List<Integer> singleton = new ArrayList<Integer>();
            singleton.add(i);
            cells.add(singleton);
        }
        double[][] C = new double[M][M];
        for (int i = 0; i < M; i++) {
            System.arraycopy(A[i], 0, C[i], 0, M);
        }

        while (cells.size() > ncellsTarget) {
            int n = cells.size();
            double maxc = -1.0;
            int besta = -1, bestb = -1;
            for (int a = 0; a < n; a++) {
                for (int b = a + 1; b < n; b++) {
                    if (cells.get(a).size() + cells.get(b).size() <= 2 * cellsize && C[a][b] > maxc) {
                        maxc = C[a][b];
                        besta = a;
                        bestb = b;
                    }
                }
            }
            if (besta == -1) {
                // Every merge exceeds the size cap: merge the two smallest cells.
                int s1 = -1, s2 = -1;
                for (int a = 0; a < n; a++) {
                    if (s1 == -1 || cells.get(a).size() < cells.get(s1).size()) {
                        s2 = s1;
                        s1 = a;
                    } else if (s2 == -1 || cells.get(a).size() < cells.get(s2).size()) {
                        s2 = a;
                    }
                }
                besta = FastMath.min(s1, s2);
                bestb = FastMath.max(s1, s2);
            }

            cells.get(besta).addAll(cells.get(bestb));
            for (int j = 0; j < n; j++) {
                C[besta][j] = C[besta][j] + C[bestb][j];
                C[j][besta] = C[j][besta] + C[j][bestb];
            }
            C[besta][besta] = 0.0;

            cells.remove(bestb);
            int m = cells.size();
            double[][] Cnew = new double[m][m];
            int ri = 0;
            for (int i = 0; i < n; i++) {
                if (i == bestb) {
                    continue;
                }
                int cj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == bestb) {
                        continue;
                    }
                    Cnew[ri][cj] = C[i][j];
                    cj++;
                }
                ri++;
            }
            C = Cnew;
        }

        List<int[]> out = new ArrayList<int[]>();
        for (int c = 0; c < cells.size(); c++) {
            List<Integer> cell = cells.get(c);
            int[] arr = new int[cell.size()];
            for (int p = 0; p < arr.length; p++) {
                arr[p] = cell.get(p);
            }
            out.add(arr);
        }
        return out;
    }

    private CellTraj integrateCell(FirstOrderDifferentialEquations ode, double t0, double[] y0c, double t1,
                                   SolverOptions options, double[] yDefault, int[] idx) {
        int d = y0c.length;
        double[] next = new double[d];

        if (options.stiff) {
            LSODA solver = options.odesolvers.stiffIntegratorFor(t0, t1, options.tol,
                    options.tol > GlobalConstants.CoarseTol);
            try {
                solver.integrate(ode, t0, y0c, t1, next);
            } catch (RuntimeException firstFailure) {
                // See ClosingAndStateDepMethodsAnalyzer: the retry from yDefault is the
                // probe that decides whether the initial point was the cause, so the
                // message belongs after it succeeds, not before it runs.
                try {
                    solver.integrate(ode, t0, gather(yDefault, idx), t1, next);
                } catch (RuntimeException retryFailure) {
                    retryFailure.addSuppressed(firstFailure);
                    throw new RuntimeException("TBI fluid integration failed from BOTH the supplied"
                            + " initial point and the default initialization over t in ["
                            + t0 + ", " + t1 + "]; the initial point is NOT implicated."
                            + " Underlying integrator error: " + retryFailure.getMessage(),
                            retryFailure);
                }
                if (options.verbose != VerboseLevel.SILENT) {
                    System.out.println("The initial point was invalid, Fluid solver switched to default initialization.");
                }
            }
            return fromLsoda(solver, d);
        }

        FirstOrderIntegrator solver = options.odesolvers.integratorFor(t0, t1, options.tol,
                options.tol > GlobalConstants.CoarseTol);
        solver.clearStepHandlers();
        TransientDataHandler handler = new TransientDataHandler(d);
        solver.addStepHandler(handler);

        boolean usedStiffFallback = false;
        try {
            try {
                solver.integrate(ode, t0, y0c, t1, next);
            } catch (RuntimeException e) {
                if (e.getMessage() != null && e.getMessage().contains("step size")) {
                    usedStiffFallback = true;
                } else {
                    solver.clearStepHandlers();
                    handler = new TransientDataHandler(d);
                    solver.addStepHandler(handler);
                    // Message deferred until the retry succeeds: a failure here falls through
                    // to the stiff fallback below, which means the initial point was NOT
                    // what went wrong.
                    solver.integrate(ode, t0, gather(yDefault, idx), t1, next);
                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("The initial point was invalid, Fluid solver switched to default initialization.");
                    }
                }
            }
        } catch (RuntimeException e) {
            usedStiffFallback = true;
        }

        if (usedStiffFallback) {
            LSODA stiffSolver = options.odesolvers.stiffIntegratorFor(t0, t1, options.tol,
                    options.tol > GlobalConstants.CoarseTol);
            stiffSolver.integrate(ode, t0, y0c, t1, next);
            return fromLsoda(stiffSolver, d);
        }

        int n = handler.tVec.getNumRows();
        double[] t = new double[n];
        double[][] y = new double[n][d];
        for (int r = 0; r < n; r++) {
            t[r] = handler.tVec.get(r, 0);
            for (int j = 0; j < d; j++) {
                y[r][j] = FastMath.max(0.0, handler.xVec.get(r, j));
            }
        }
        return new CellTraj(t, y);
    }

    private static CellTraj fromLsoda(LSODA solver, int d) {
        java.util.ArrayList<Double> tHistory = solver.getTvec();
        java.util.ArrayList<Double[]> yHistory = solver.getYvec();
        int n = solver.getStepsTaken() + 1;
        double[] t = new double[n];
        double[][] y = new double[n][d];
        for (int r = 0; r < n; r++) {
            t[r] = tHistory.get(r);
            for (int j = 0; j < d; j++) {
                y[r][j] = FastMath.max(0.0, yHistory.get(r)[j]);
            }
        }
        return new CellTraj(t, y);
    }

    private static double[] gather(double[] x, int[] idx) {
        double[] out = new double[idx.length];
        for (int p = 0; p < idx.length; p++) {
            out[p] = x[idx[p]];
        }
        return out;
    }

    /** Merge two ascending, strictly-increasing time grids into their sorted union. */
    private static double[] unionSorted(double[] a, double[] b) {
        double[] merged = new double[a.length + b.length];
        int i = 0, j = 0, k = 0;
        while (i < a.length && j < b.length) {
            if (a[i] < b[j]) {
                merged[k++] = a[i++];
            } else if (a[i] > b[j]) {
                merged[k++] = b[j++];
            } else {
                merged[k++] = a[i++];
                j++;
            }
        }
        while (i < a.length) {
            merged[k++] = a[i++];
        }
        while (j < b.length) {
            merged[k++] = b[j++];
        }
        double[] out = new double[k];
        System.arraycopy(merged, 0, out, 0, k);
        return out;
    }

    /** Cell trajectory: ascending time grid and per-time state rows. */
    private static class CellTraj {
        final double[] t;
        final double[][] y;

        CellTraj(double[] t, double[][] y) {
            this.t = t;
            this.y = y;
        }
    }

    /**
     * ODE of a single cell: the full state is the frozen-complement trajectory
     * with the cell entries overridden by the live cell state; the cell
     * derivative is the restriction of the closing-method derivative to the
     * cell's state indices.
     */
    private static class CellODE implements FirstOrderDifferentialEquations {
        private final PassageTimeODE pt;
        private final Matrix allJumps;
        private final int[] idx;
        private final int ndim;
        private final double[] frozenT;
        private final double[][] frozenY;

        CellODE(PassageTimeODE pt, Matrix allJumps, int[] idx, int ndim,
                double[] frozenT, double[][] frozenY) {
            this.pt = pt;
            this.allJumps = allJumps;
            this.idx = idx;
            this.ndim = ndim;
            this.frozenT = frozenT;
            this.frozenY = frozenY;
        }

        @Override
        public int getDimension() {
            return idx.length;
        }

        @Override
        public void computeDerivatives(double t, double[] xc, double[] dxc) {
            double[] xfull = FluidInterp.interp(frozenT, frozenY, t, ndim);
            for (int p = 0; p < idx.length; p++) {
                double v = xc[p];
                xfull[idx[p]] = (v < 0.0) ? 0.0 : v;
            }
            Matrix rates = pt.calculateRatesClosing(t, xfull);
            Matrix full = allJumps.mult(rates, null);
            for (int p = 0; p < idx.length; p++) {
                dxc[p] = full.get(idx[p], 0);
            }
        }
    }
}
