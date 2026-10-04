package jline.api.infer;

import jline.lang.layered.LayeredNetwork;
import jline.solvers.LayeredNetworkAvgTable;
import jline.util.matrix.Matrix;

/**
 * Options for LQN parameter identification via the Extended Kalman Filter
 * ({@link InferLqn}).
 *
 * <p>The tuning factors follow Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking
 * Time-Varying Parameters in Software Systems with Extended Kalman Filters",
 * CASCON 2005: the drift covariance is Q_ii = (QFac*a0_i*cvA)^2 (eq 9a) and the
 * measurement covariance is R_ii = ((RFac*zbar_i)/1.96)^2/gammaT (eq 9b), with
 * gammaT = T/Tstar. Explicit Q, R, P0, a0 override the constructed values.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London. All rights reserved.
 */
public class InferLqnOptions {

    /** Observation model: solve the LQN and return its per-element average table. */
    public interface ObservationSolver {
        LayeredNetworkAvgTable evaluate(LayeredNetwork model);
    }

    /** Observation-model solver; when null a silent SolverLN is used. */
    public ObservationSolver solver = null;

    public double QFac = 0.1;    // drift-noise factor (eq 9a)
    public double RFac = 0.2;    // measurement-noise factor (eq 9b)
    public double cvA = 1.0;     // parameter drift coefficient of variation

    public Double gammaT = null; // T/Tstar ratio for R; else derived from T,Tstar; else 1
    public Double T = null;
    public Double Tstar = null;

    public Matrix Q = null;      // explicit covariances (override the tuning factors)
    public Matrix R = null;
    public Matrix P0 = null;
    public Matrix a0 = null;     // initial estimate (default: current model values)
    public Matrix aTrue = null;  // ground truth, enables the Ea RMS metric

    public double fdStep = 1e-3;     // finite-difference relative step
    public double fdFloor = 1e-6;    // finite-difference / positivity floor
    public boolean clampPositive = true;
    public boolean verbose = false;
}
