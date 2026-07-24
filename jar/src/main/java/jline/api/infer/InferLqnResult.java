package jline.api.infer;

import jline.util.matrix.Matrix;

import java.util.List;

/**
 * Result of an LQN parameter identification run ({@link InferLqn}).
 *
 * <p>{@link #ahat} holds the parameter estimate trajectory (np x nsteps).
 * {@link #Er} is the RMS prediction error; {@link #Ea} the RMS parameter
 * tracking error against the supplied ground truth (null when none is given),
 * as defined in Section 4 of the CASCON 2005 reference.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London. All rights reserved.
 */
public class InferLqnResult {

    public Matrix ahat;        // np x nsteps parameter estimate trajectory
    public Matrix e;           // no x nsteps prediction errors
    public Matrix zpred;       // no x nsteps predicted measurements
    public Matrix P;           // final estimation-error covariance (np x np)
    public List<Matrix> Phist; // per-step covariances

    public double Er;          // RMS prediction error
    public Double Ea;          // RMS parameter tracking error (null if no aTrue)

    public Matrix a0;          // initial estimate used
    public Matrix Q;           // drift covariance used
    public Matrix R;           // measurement covariance used
    public Matrix P0;          // initial covariance used
}
