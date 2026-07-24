/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Raw moment results behind {@link NetworkMomentTable}, mirroring the {@code mom}
 * struct that MATLAB's {@code getMomentTable} returns as its second output.
 *
 * <p>The table shows only the diagonal of the queue-length covariance; this holder
 * carries the full per-station covariance blocks, the underlying
 * {@link jline.api.pfqn.sens.Pfqn_sens_mva} / {@link jline.api.pfqn.sens.Pfqn_sens_mvaldmx}
 * result, and the {@link jline.api.pfqn.sens.Pfqn_sens_respt} result when the
 * sojourn-time moments were computed. Use it when the per-pair covariances are
 * needed.</p>
 *
 * <p>Exactly one of {@code qlenMva} and {@code qlenLdmx} is non-null on a closed or
 * mixed model; both are null on a purely open model, whose queue-length law is
 * evaluated in closed form and has no lattice recursion behind it. The
 * convenience views {@code Q}, {@code QVar} and {@code X} are populated in every
 * case that defines them; {@code QCov} is null on the purely open branch and on
 * models whose queue-length path does not produce covariance blocks.</p>
 *
 * @see NetworkMomentTable
 * @see NetworkSolver#getMomentTable()
 */
public class NetworkMomentResult {

    /**
     * The closed single-server queue-length result, or null.
     */
    public Ret.pfqnSensMva qlenMva;

    /**
     * The closed-multiserver or mixed queue-length result, or null.
     */
    public Ret.pfqnSensMvaldmx qlenLdmx;

    /**
     * The sojourn-time moment result, or null when response-time moments were not
     * available (open models, or FCFS rates that are not class-independent).
     */
    public Ret.pfqnSensRespt respt;

    /**
     * The per-class queue-length moment result, from {@code Pfqn_sens_mom} with
     * {@code groups = 1..R}, i.e. one class per group. This is the source of the
     * QLenSkew column: scaling L(i,r) alone is Akyildiz-Strelen Theorem 1 with the
     * class subset T = {r}, so it generates the moments of n(i,r) itself rather than
     * of the station total.
     *
     * <p>Null unless moment order 3 was requested AND the model is closed
     * single-server, which is the scope of {@code Pfqn_sens_mom}.</p>
     */
    public Ret.pfqnSensMom qlenmom;

    /**
     * 1 x R throughput of the queue-length path, or null on the purely open branch.
     */
    public Matrix X;

    /**
     * Mq x R mean queue length, one row per queueing station in node order.
     */
    public Matrix Q;

    /**
     * Mq x R queue-length variance, Var[n(i,r)].
     */
    public Matrix QVar;

    /**
     * Per-station covariance blocks; {@code QCov[i]} is R x R with entry (r,s) =
     * Cov[n(i,r),n(i,s)]. Null on the purely open branch.
     */
    public Matrix[] QCov;

    /**
     * Whether a queue-length recursion was run, i.e. whether the model was closed
     * or mixed rather than purely open.
     */
    public boolean hasQlen() {
        return qlenMva != null || qlenLdmx != null;
    }
}
