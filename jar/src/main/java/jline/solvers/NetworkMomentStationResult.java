/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Raw moment results behind {@link NetworkMomentStationTable}, mirroring the
 * {@code mom} struct that MATLAB's {@code getMomentStationTable} returns as its
 * second output.
 *
 * <p>Carries the per-station-total moments together with the cross-station
 * covariance matrix {@code Cov} that the table does not show, and the underlying
 * {@link jline.api.pfqn.sens.Pfqn_sens_mom} (exact) or
 * {@link jline.api.pfqn.sens.Pfqn_sens_linearizer} (approximate) result.</p>
 *
 * @see NetworkMomentStationTable
 * @see NetworkSolver#getMomentStationTable(int[])
 */
public class NetworkMomentStationResult {

    /**
     * The exact result, non-null when the method was {@code "exact"}.
     */
    public Ret.pfqnSensMom exact;

    /**
     * The Linearizer result, non-null when the method was {@code "lin"}.
     */
    public Ret.pfqnSensLinearizer linearizer;

    /**
     * M x 1, m(i) = E[Q_i], the mean total queue length at station i.
     */
    public Matrix m;

    /**
     * M x 1, Var[Q_i].
     */
    public Matrix Var;

    /**
     * M x M, Cov(i,j) = Cov[Q_i,Q_j]; the diagonal is {@code Var}. Non-null only in
     * the single-group case, i.e. always for the station table and for the chain
     * table of a single-chain model, mirroring the collapse of the MATLAB reference.
     */
    public Matrix Cov;

    /**
     * The general (M x G x M x G) covariance, {@code CovG[i][g]} being an M x G
     * matrix whose (j,g2) entry is Cov[Q_(i,g),Q_(j,g2)]. Populated whenever the
     * exact path ran, and the only form available when there is more than one group
     * (more than one chain), where cross-chain and cross-station covariances both
     * exist. Null on the Linearizer path, which approximates station totals only.
     */
    public Matrix[][] CovG;

    /**
     * M x 1, E[Q_i^2].
     */
    public Matrix M2;

    /**
     * M x 1, E[Q_i^3].
     */
    public Matrix M3;

    /**
     * M x 1, skewness of Q_i; NaN where the variance vanishes.
     */
    public Matrix Skew;

    /**
     * Asymmetry residual of the underlying recursion. Roundoff for the exact
     * method; a genuine measure of approximation error for the Linearizer, whose
     * fixed point does not enforce the symmetry the product form guarantees.
     */
    public double CovAsym;
}
